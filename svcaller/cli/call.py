import logging, os, tempfile, uuid

import click
import pysam
import pybedtools

from svcaller.calling.events import SvType
from svcaller.calling.events import event_filt, clust_filt, call_events, \
    filter_on_shared_termini, chrom_int2str, event_termini_spaced_broadly


def remove_bam_and_bai(bam_filename):
    if os.path.exists(bam_filename):
        os.remove(bam_filename)

    bai_filename = "{}.bai".format(bam_filename)
    if os.path.exists(bai_filename):
        os.remove(bai_filename)


@click.command()
@click.argument('input-bam', type=click.Path(exists=True))
@click.option('--event-type', type=click.Choice([SvType.DEL.value, SvType.INV.value, SvType.TRA.value, SvType.DUP.value]), required=True)
@click.option('--fasta-filename', type=str, required=True)
@click.option('--filter-event-overlap', is_flag = True)
@click.option('--events-gtf', type=str, default = "events.gtf", required=False)
@click.option('--events-bam', type=str, default = "events.bam", required=False)
@click.option('--tmp-dir', type=str, default = "/tmp", required=False)
def run_all_cmd(input_bam, event_type, fasta_filename, events_gtf, events_bam,
                filter_event_overlap, tmp_dir):
    logging.info("Running event calling for event type {} on {}.".format(event_type, input_bam))

    unique_id = uuid.uuid4()
    logging.info("Applying event-specific filter...")
    event_filt_reads_bam_name = "{}/{}_event_filter_{}.bam".format(tmp_dir, unique_id, event_type)
    event_filter_inner(event_type, event_filt_reads_bam_name, input_bam)

    logging.info("Applying read-cluster filter...")
    clust_filt_reads_bam_name = "{}/{}_cluster_filter_{}.bam".format(tmp_dir, unique_id, event_type)
    cluster_filter_inner(clust_filt_reads_bam_name, event_filt_reads_bam_name)

    logging.info("Calling events...")
    with open(events_gtf, 'w') as events_outfile:
        call_events_inner(clust_filt_reads_bam_name, event_type, fasta_filename, events_outfile,
                          events_bam, filter_event_overlap, tmp_dir)

    remove_bam_and_bai(clust_filt_reads_bam_name)
    remove_bam_and_bai(event_filt_reads_bam_name)


@click.command()
@click.option('--event-type', type=click.Choice([SvType.DEL.value, SvType.INV.value, SvType.TRA.value, SvType.DUP.value]), required=True)
@click.option('--output-bam', type=str, default="event_filtered_reads.bam", required=False)
@click.argument('input-bam', type=click.Path(exists=True))
def event_filter_cmd(event_type, output_bam, input_bam):
    logging.info("TRACE: HERE.")
    event_filter_inner(event_type, output_bam, input_bam)


def parse_xa_tok(tok):
    elems = tok.split(",")
    chrom_str = elems[0] # XXX Need to convert ints to X and Y?
    strand = elems[1][0]
    pos = int(elems[1][1:])
    return (chrom_str, pos, strand)


def parse_xa_flag(xa_tag_string):
    toks = xa_tag_string.strip(";").split(";")
    return [parse_xa_tok(tok) for tok in toks]


def positions_consistent_with_no_alternation(read_1, read_2):
    """
    :param read_1: (chrom, pos, strand) tuple
    :param read_2: (chrom, pos, strand) tuple
    :return: True if the read pair can easily be explained without assuming a
    structural alteration, False otherwise.
    """

    read_1_chrom = read_1[0]
    read_1_pos = read_1[1]
    read_1_strand = read_1[2]

    read_2_chrom = read_2[0]
    read_2_pos = read_2[1]
    read_2_strand = read_2[2]

    if read_1_chrom == read_2_chrom:
        if read_1_strand == "+" and read_2_strand == "-":
            if read_1_pos < read_2_pos:
                if read_2_pos - read_1_pos < 1000:
                    return True
                else:
                    return False
            else:
                return False
        elif read_1[2] == "-" and read_2[2] == "+":
            if read_1_pos > read_2_pos:
                if read_1_pos - read_2_pos < 1000:
                    return True
                else:
                    return False
            else:
                return False
    else:
        return False


def no_parsimonious_alternative_mappings(read):
    """
    :param read: A pysam AlignedSegment
    :return: True if this AlignedSegment and it's read pair have an alternative alignment that
    is consistent with the read pairs deriving from a standard genome with no structural alteration
    affecting the read pair.
    """

    if read.has_tag("XA"):
        # Extract strand string for the paired read:
        if not(read.flag & 32):
            paired_read_strand = "+"
        else:
            paired_read_strand = "-"

        # Extract mapping coordinate info for the read pair:
        pair_mapping = (chrom_int2str(read.rnext), read.mpos, paired_read_strand)

        xa_tag_string = read.get_tag("XA")
        alternative_mappings = parse_xa_flag(xa_tag_string)
        for alternative_mapping in alternative_mappings:
            if positions_consistent_with_no_alternation(pair_mapping, alternative_mapping):
                return False

    return True


def event_filter_inner(event_type, output_bam, input_bam):
    logging.info("Running event filter: {}".format(event_type))

    samfile = pysam.AlignmentFile(input_bam, "rb")
    event_filtered_reads = event_filt(samfile, event_type)

    # Impose an additional filter, to throw away reads where the read pair
    # includes an alternative mapping that does not require a structural variant
    # in order to be explained:
    extra_filtered_reads = list(filter(lambda read: no_parsimonious_alternative_mappings(read),
                                       event_filtered_reads))

    with pysam.AlignmentFile(output_bam, "wb", header=samfile.header) as outf:
        for read in extra_filtered_reads:
            outf.write(read)

    pysam.index(str(output_bam))


@click.command()
@click.option('--output-bam', type=str, default="cluster_filtered_reads", required=False)
@click.argument('input-bam', type=click.Path(exists=True))
def cluster_filter_cmd(output_bam, input_bam):
    cluster_filter_inner(output_bam, input_bam)


def cluster_filter_inner(output_bam, input_bam):
    logging.info("Running cluster filter.")

    samfile = pysam.AlignmentFile(input_bam, "rb")
    reads = [r for r in list(samfile)]
    cluster_filtered_reads = clust_filt(reads)
    with pysam.AlignmentFile(output_bam, "wb", header=samfile.header) as outf:
        for read in cluster_filtered_reads:
            outf.write(read)

    pysam.index(str(output_bam))


@click.command()
@click.argument('input-bam', type=click.Path(exists=True))
@click.argument('event-type', type=click.Choice([SvType.DEL.value, SvType.DUP.value, SvType.INV.value, SvType.TRA.value]))
@click.option('--fasta-filename', type=str, required=True)
@click.option('--events-gtf', type=str, default = "events.gtf", required=False)
@click.option('--events-bam', type=str, default = "events.bam", required=False)
@click.option('--filter-event-overlap', is_flag = True)
@click.option('--tmp-dir', type=str, default = "/tmp", required=False)
def call_events_cmd(input_bam, event_type, fasta_filename, events_gtf, events_bam, filter_event_overlap, tmp_dir):
    events_outfile = open(events_gtf, 'w')
    call_events_inner(input_bam, event_type, fasta_filename, events_outfile, events_bam, filter_event_overlap, tmp_dir)


def call_events_inner(filtered_bam, event_type, fasta_filename, events_gff, events_bam, filter_event_overlap, tmp_dir):
    logging.info("Calling events on file {}:".format(filtered_bam))

    samfile = pysam.AlignmentFile(filtered_bam, "rb")
    filtered_reads = [r for r in list(samfile)]

    # Call events:
    pybedtools.set_tempdir(tmp_dir) # Necessary to control temporary folder usage during event calling
    events = list(call_events(filtered_reads, fasta_filename))

    # Filter on discordant read support depth:
    logging.info("Filtering on discordant read support depth...")
    filtered_events = list(filter(lambda event: (len(event._terminus1_reads) >= 3 and len(event._terminus2_reads) >= 3), events))

    # Filter on maximum quality of terminus reads:
    logging.info("Filtering on max quality of terminus reads...")
    filtered_events = list(filter(lambda event: (event.get_t1_mapqual() >= 19 and event.get_t2_mapqual() >= 19), filtered_events))

    logging.info("Filtering on terminus read spacing...")
    filtered_events = list(filter(lambda event: event_termini_spaced_broadly(event), filtered_events))

    if filter_event_overlap:
        # Filter on event terminus sharing (exclude any events that have
        # overlapping termini):
        logging.info("Filtering on terminus event sharing...")
        filtered_events = filter_on_shared_termini(filtered_events)

    # Filter on soft-clipping support:
    logging.info("Filtering on soft-clip support...")
    filtered_events = list(filter(lambda event: event.has_soft_clip_support(), filtered_events))

    # Optionally filter on presence of soft-clipped regions scattered throughout the reads, depending on
    # the event type:
    logging.info("Filtering on scattered soft-clip regions...")
    if event_type == SvType.DUP.value or event_type == SvType.INV.value:
        filtered_events = list(filter(lambda event: not event.has_scattered_soft_clip_regions(), filtered_events))

    # Print them out:
    logging.info("Printing final events...")
    for event in filtered_events:
        print(event.get_gtf(), file=events_gff)

    # Write to a temporary bam file, to facilitate subsequent sorting with pysam. NOTE:
    # Could do sorting in memory since the read count should be low, but it seems less
    # bug-prone to use pysam's sort functionality:
    unique_id = uuid.uuid4()
    tmp_bam_filename = "{}/penultimate_bamfile_{}.bam".format(tmp_dir, unique_id, event_type)
    with pysam.AlignmentFile(tmp_bam_filename, "wb", header=samfile.header) as outf:
        for event in filtered_events:
            for read in event._terminus1_reads + event._terminus2_reads:
                outf.write(read)

    # Sort the intermediate bam file with samtools to produce the final output bam file, then index it:
    pysam.sort("-o", events_bam, tmp_bam_filename)
    pysam.index(str(events_bam))

    remove_bam_and_bai(tmp_bam_filename)

    events_gff.close()
