import logging
import os, pdb

import click
import pysam

from svcaller.cli.calling.calling import event_filt, clust_filt, call_events, filter_on_shared_termini
from svcaller.cli.calling.calling import DEL, INV, TRA, DUP


@click.command()
@click.option('--event-type', type=click.Choice([DEL, INV, TRA, DUP]), required=True)
@click.option('--output-name', type=str, default="SV_CallOutput", required=False)
@click.argument('input-bam', type=click.Path(exists=True))
@click.pass_context
def run_all_cmd(ctx, event_type, output_name, input_bam):
    logging.info("Running event calling: {}".format(event_type))

    logging.info("Applying event-specific filter...")
    event_filt_reads_bam_name = output_name + "_filt1.bam"
    event_filter_cmd(event_type, event_filt_reads_bam_name, input_bam)

    logging.info("Applying read-cluster filter...")
    clust_filt_reads_bam_name = output_name + "_filt2.bam"
    cluster_filter_cmd(clust_filt_reads_bam_name, event_filt_reads_bam_name)

    call_events(clust_filt_reads_bam_name)


@click.command()
@click.option('--event-type', type=click.Choice([DEL, INV, TRA, DUP]), required=True)
@click.option('--output-bam', type=str, default="event_filtered_reads.bam", required=False)
@click.argument('input-bam', type=click.Path(exists=True))
@click.pass_context
def event_filter_cmd(ctx, event_type, output_bam, input_bam):
    event_filter_inner(event_type, output_bam, input_bam)


def event_filter_inner(event_type, output_bam, input_bam):
    logging.info("Running event filter: {}".format(event_type))

    samfile = pysam.AlignmentFile(input_bam, "rb")
    event_filtered_reads = event_filt(samfile, event_type)

    with pysam.AlignmentFile(output_bam, "wb", header=samfile.header) as outf:
        for read in event_filtered_reads:
            outf.write(read)

    pysam.index(str(output_bam))


@click.command()
@click.option('--output-bam', type=str, default="cluster_filtered_reads", required=False)
@click.argument('input-bam', type=click.Path(exists=True))
@click.pass_context
def cluster_filter_cmd(ctx, output_bam, input_bam):
    cluster_filter_inner(output_bam, input_bam)


def cluster_filter_inner(output_bam, input_bam):
    logging.info("Running cluster filter.")

    samfile = pysam.AlignmentFile(input_bam, "rb")
    reads = [r for r in list(samfile)]
    cluster_filtered_reads = clust_filt(reads, samfile)
    with pysam.AlignmentFile(output_bam, "wb", header=samfile.header) as outf:
        for read in cluster_filtered_reads:
            outf.write(read)

    pysam.index(str(output_bam))


@click.command()
@click.argument('input-bam', type=click.Path(exists=True))
@click.option('--fasta-filename', type=str, required=True)
@click.option('--events-gff', type=str, default = "events.gff", required=False)
@click.option('--softclips-gff', type=str, default = "softclippings.gff", required=False)
@click.pass_context
def call_events_cmd(ctx, input_bam, fasta_filename, events_gff, softclips_gff):
    events_outfile = open(events_gff, 'w')
    softclips_outfile = open(softclips_gff, 'w')
    output = call_events_inner(input_bam, fasta_filename, events_outfile, softclips_outfile)
    return output


def call_events_inner(filtered_bam, fasta_filename, events_gff, softclips_gff):
    logging.info("Calling events on file {}:".format(filtered_bam))

    samfile = pysam.AlignmentFile(filtered_bam, "rb")
    filtered_reads = [r for r in list(samfile)]

    # Call events:
    events = list(call_events(filtered_reads, fasta_filename))

    # Filter on discordant read support depth:
    filtered_events = filter(lambda event: (len(event._terminus1_reads) > 5 and len(event._terminus2_reads) > 5), events)

    # Filter on event terminus sharing (exclude any events that have
    # overlapping termini):
    filtered_events = filter_on_shared_termini(filtered_events)

    # Filter on soft-clipping support:
    filtered_events = filter(lambda event: event.has_soft_clip_support(), filtered_events)

    # Print them out:
    for event in filtered_events:
        print >> events_gff, event.to_string()

    for event in filtered_events:
        print >> softclips_gff, event.get_softclip_strings()

    events_gff.close()
    softclips_gff.close()
