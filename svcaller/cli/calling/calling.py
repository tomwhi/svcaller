# coding=utf-8
from datetime import datetime
import functools, logging, sys
import svcaller.cli.calling.softclip as softclip
from collections import defaultdict

DEL = "DEL"
INV = "INV"
TRA = "TRA"
DUP = "DUP"


def calc_cigar_bps(cigar_tup):
    """
    Calculate the effective number of base-pairs specified by this cigar tuple.

    :param cigar_tup: Cigar code and length
    :return: Number of base-pairs spanned
    """

    if cigar_tup[0] == 1:
        # Insertion
        return 0
    elif cigar_tup[0] == 2:
        # Deletion
        return 1 * cigar_tup[1]
    elif cigar_tup[0] == 3:
        # "N": Not defined for DNA
        raise Exception("Unexpected 'N' value in cigar.")
    elif cigar_tup[0] == 4:
        # Soft-clipping
        return 1 * cigar_tup[1]
    elif cigar_tup[0] == 5:
        # Hard-clipping: Not yet sure how to deal with these instances
        raise Exception("Unexpected 'H' value in cigar.")
    elif cigar_tup[0] == 7:
        # Sequence match
        return 1 * cigar_tup[1]
    elif cigar_tup[0] == 7:
        # Sequence mis-match
        return 1 * cigar_tup[1]


def calculate_genomic_span(read):
    """
    Calculate the genomic region spanned by the read, with non-mapping (e.g. soft-clipped)
    regions factored in.

    :param read: A pysam read
    :return: A tuple containing read span start and end integer positions
    """

    cigar = read.cigar

    # Calculate the number of base-pairs to prepend at the start of the alignment
    # (i.e. at the region preceding the first "matched" region):
    idx_of_first_match = min([idx for idx in range(len(cigar)) if cigar[idx][0] == 0])
    cigar_preceding_match = cigar[:idx_of_first_match]
    bps_to_prepend_list = [calc_cigar_bps(tup) for tup in cigar_preceding_match]
    bps_to_prepend = 0
    for bp in bps_to_prepend_list:
        bps_to_prepend += bp

    idx_of_last_match = max([idx for idx in range(len(cigar)) if cigar[idx][0] == 0])
    cigar_after_match = cigar[idx_of_last_match + 1:]
    bps_to_append_list = [calc_cigar_bps(tup) for tup in cigar_after_match]
    bps_to_append = 0
    for bp in bps_to_append_list:
        bps_to_append += bp

    first_match_position = read.pos
    last_match_position = read.reference_end

    return (first_match_position - bps_to_prepend,
            last_match_position + bps_to_append)


def event_termini_spaced_broadly(event):
    """
    Indicates whether the given event contains reads that are mapped to a
    broad region, for both termini. This can be used to filter out a
    particular class of apparent false-positive events (mostly TRA events).

    :param event: An event 
    :return: False if more than 50% of the read-pairs have at least one terminus
    where the start positionis identical, True otherwise.
    """

    terminus1_spans = list(map(lambda read: calculate_genomic_span(read), event._terminus1_reads))
    terminus2_spans = list(map(lambda read: calculate_genomic_span(read), event._terminus2_reads))

    terminus1_starts = [tup[0] for tup in terminus1_spans]
    terminus1_start_counts = {start:terminus1_starts.count(start) for start in terminus1_starts}
    terminus1_ends = [tup[1] for tup in terminus1_spans]
    terminus1_end_counts = {end:terminus1_ends.count(end) for end in terminus1_ends}

    terminus2_starts = [tup[0] for tup in terminus2_spans]
    terminus2_start_counts = {start:terminus2_starts.count(start) for start in terminus2_starts}
    terminus2_ends = [tup[1] for tup in terminus2_spans]
    terminus2_end_counts = {end:terminus2_ends.count(end) for end in terminus2_ends}

    return (max(terminus1_start_counts.values()) <= (len(terminus1_spans)/2.0) and \
            max(terminus1_end_counts.values()) <= (len(terminus1_spans)/2.0) and \
            max(terminus2_start_counts.values()) <= (len(terminus2_spans)/2.0) and \
            max(terminus2_end_counts.values()) <= (len(terminus2_spans)/2.0))


def event_filt(read_iter, event_type, flag_filter=4+8+256+1024+2048):
    filtered_reads = []
    for read in read_iter:
        if not read.flag & flag_filter:
            # This read passes the general bit-wise filter.
            # Apply the event-specific bit-wise flag field and other field
            # filters:
            if event_type == DEL:
                del_filt(filtered_reads, read)
            elif event_type == INV:
                inv_filt(filtered_reads, read)
            elif event_type == TRA:
                tra_filt(filtered_reads, read)
            elif event_type == DUP:
                dup_filt(filtered_reads, read)
            else:
                logging.error("Invalid event_type: " + event_type)

    return filtered_reads


def del_filt(filtered_reads, read):
    # Select reads facing each other, on the same chromosome, with a gap
    # larger than the specified distance threshold:
    if (read.rname == read.rnext):
        if not(read.flag & 16) and (read.flag & 32):
            if read.tlen > 1000:
                filtered_reads.append(read)
        elif (read.flag & 16) and not(read.flag & 32):
            if read.tlen < -1000:
                filtered_reads.append(read)


def dup_filt(filtered_reads, read):
    # Select reads facing away from each other, on the same chromosome. Impose a
    # minimum event size to avoid small artefact calls (of currently unknown origin):
    insert_size_abs = abs(read.tlen)
    if insert_size_abs > 1000:
        if (read.rname == read.rnext):
            if not(read.flag & 16) and (read.flag & 32):
                if read.tlen < 0:
                    filtered_reads.append(read)
            elif (read.flag & 16) and not(read.flag & 32):
                if read.tlen > 0:
                    filtered_reads.append(read)


def inv_filt(filtered_reads, read):
    # Select reads facing in the same direction, on the same chromosome. Impose a
    # minimum event size to avoid small artefact calls (of currently unknown origin):

    insert_size_abs = abs(read.tlen)
    if insert_size_abs > 1000:
        if (read.rname == read.rnext):
            if not(read.flag & 16) and not(read.flag & 32):
                filtered_reads.append(read)
            elif (read.flag & 16) and (read.flag & 32):
                filtered_reads.append(read)


def tra_filt(filtered_reads, read, min_map_qual = 50):
    # First, filter on mapping quality specifically for TRA events, which are
    # prone to mapping artefacts:
    if read.mapping_quality >= min_map_qual and read.get_tag("MQ") >= min_map_qual:
        # Select reads located on distinct chromosomes:
        if (read.rname != read.rnext):
            filtered_reads.append(read)


# FIXME: Quick hack:
def chrom_int2str(chrom_int):
    chrom_str = str(chrom_int + 1)
    if chrom_str == "23":
        chrom_str = "X"
    elif chrom_str == "24":
        chrom_str = "Y"

    return chrom_str


def clust_filt(read_iter, samfile, chrom_of_interest=22):
    '''Filter the input reads to only retain those read-pairs where there are other read-pairs with
    similar coordinates. This is an ugly implementation, needs refactoring.

    chrom_of_interest is a chromosome that is expected to contain
    more reads than the other chromosomes.'''
    idx = 0
    filtered_reads = []
    for read in read_iter:
        half_read_len = int((len(read.seq)/2.0))

        # Determine which of the reads (for this read-pair) will be used to identify
        # nearby read-pairs: We could look at this read, or it's mate-pair.
        base_read_chrom = read.rname
        base_read_pos = read.pos
        other_read_chrom = read.rnext
        other_read_pos = read.pnext

        # PROBLEM: This is a hack making the analysis focused on the X-chromosome.
        # Need to revisit the problem and figure out a more general solution.

        if base_read_chrom == 22:
            base_read_chrom = read.rnext
            base_read_pos = read.pnext
            other_read_chrom = read.rname
            other_read_pos = read.pos

        base_read_chrom_str = chrom_int2str(base_read_chrom)

        try:
            # Only consider reads nearby that are not optical/pcr duplicates and that are not secondary alignments
            # (hence the flag filter):
            reads_nearby_non_unique = [r for r in samfile.fetch(base_read_chrom_str,
                base_read_pos+half_read_len-200,
                base_read_pos+half_read_len+200) if r.flag < 1024]

            reads_nearby_dict = {}
            for curr_read in reads_nearby_non_unique:
                if not curr_read.qname in reads_nearby_dict:
                    reads_nearby_dict[curr_read.qname] = curr_read

            reads_nearby = reads_nearby_dict.values()

            paired_matches = 0
            for read_nearby in reads_nearby:
                read_nearby_pair_chrom = read_nearby.rnext
                if read_nearby_pair_chrom == other_read_chrom and \
                   other_read_pos+half_read_len-1000 < read_nearby.pnext+half_read_len and \
                   other_read_pos+half_read_len+1000 > read_nearby.pnext+half_read_len:
                    paired_matches += 1

            if paired_matches >= 2:
                filtered_reads.append(read)

        except Exception:
            pass

        if idx % 1000 == 0:
            print(idx, file=sys.stderr)
        idx += 1

    # Finally, filter to only retain reads where both reads in the pair are
    # included in this list:
    return paired_filt(filtered_reads)


def paired_filt(reads):
    '''Filter the specified reads, to only retain those where both read-pairs are present
    in the input list.'''

    readname2reads = get_readname2reads(reads)
    return [read for read in reads if len(readname2reads[read.qname]) == 2]


def get_readname2reads(reads):
    readname2reads = {}
    for read in reads:
        if not read.qname in readname2reads:
            readname2reads[read.qname] = []
        readname2reads[read.qname].append(read)

    return readname2reads


def pair_clusters(clusters):
    '''Pair clusters with each other: O(N)'''

    # NOTE: Currently, performing multiple passes over the reads, since that
    # way I can just take the set of clusters as input. Probably won't change
    # this, as this step should not be a time bottleneck anyway.

    # FIXME: A bunch of hacky code to get the data structures required for
    # efficient processing of this data. Could be error-prone / hard to read:
    reads = []
    for cluster in clusters:
        reads = reads + list(cluster.get_reads())

    read2cluster = {}
    for cluster in clusters:
        for read in list(cluster.get_reads()):
            read2cluster[read] = cluster

    readname2reads = get_readname2reads(reads)

    read2mate = {}
    for read in reads:
        all_reads_with_this_name = readname2reads[read.qname]
        other_reads = list(filter(lambda curr_read: curr_read != read, all_reads_with_this_name))
        assert len(other_reads) == 1
        read2mate[read] = other_reads[0]

    # Update each cluster, to indicate a pairing to each other cluster that the mate-pair
    # belongs to, for all mate-pairs of reads in this cluster (hope that makes sense!):
    for cluster in clusters:
        cluster.set_pairings(read2mate, read2cluster)

    # FIXME: This is getting messy:
    return (read2cluster, read2mate)


class SoftClipping:
    def __init__(self, sequence, read, origin_chrom, origin_start, origin_end):
        self._chrom = None
        self._start = None
        self._end = None
        self._seq = sequence
        self._read = read
        self._origin_chrom = origin_chrom
        self._origin_start = origin_start
        self._origin_end = origin_end

    def get_chrom(self):
        return self._chrom

    def get_start(self):
        return self._start

    def get_end(self):
        return self._end

    def is_in(self, chrom, start, end, fasta_filename):
        # Mapping-based approach:
        comment = '''
        if self._chrom != chrom:
            return False
        elif self._end < start:
                return False
        elif self._start > end:
                return False
        else:
            return True'''

        pos_tup = softclip.getRefMatchPos(chrom, start, end, self._seq, fasta_filename)
        if pos_tup != None:
            self._chrom = pos_tup[0]
            self._start = pos_tup[1]
            self._end = pos_tup[2]
            print(pos_tup, file=sys.stderr)

        # Sequence alignment-based approach:
        return pos_tup != None


def get_soft_clippings(read):
    '''Extract all soft-clipping secondary mapping coordinates from this
    read.'''

    cigar_softclip_lengths = []
    dist_to_softclipped = 0
    if read.cigartuples != None:
        cigar_softclip_lengths = list(map(lambda tup: tup[1], list(filter(lambda tup: tup[0] == 4, read.cigartuples))))

        # Calculate displacement of the softclipped region from the start of the *sequence*:
        if len(cigar_softclip_lengths) > 0:
            # This is completely retarded but oh well.
            idxs = range(len(read.cigartuples))
            arr = list(filter(lambda idx: list(map(lambda tup: tup[0] == 4, read.cigartuples))[idx], idxs))
            first_softclip_idx = arr[0]
            before_softclip_tups = read.cigartuples[:first_softclip_idx]
            dist_to_softclipped = sum(list(map(lambda tup: tup[1], before_softclip_tups)))

    softclip_starts = list(map(lambda tup: tuple(tup[1].split(",")[:2]), \
        list(filter(lambda tag: tag[0] == "SA", read.get_tags()))))

    # FIXME: Not sure what the data will be like. Soft-clippings do not
    # always have secondary alignments, so it may not be straightforward
    # or even well-defined to pair up the above two lists.

    # Solution: Currently, only consider reads with *exactly one* soft-clipping, for
    # simplicity:
    soft_clippings = []
    if len(cigar_softclip_lengths) == 1: #len(softclip_starts) == 1 and
        # Record the softclip stated alignment coordinates, if present,
        # although they may not be used:
        chrom = None
        start = None
        end = None
        if len(softclip_starts) == 1:
            chrom = softclip_starts[0][0]
            start = int(softclip_starts[0][1])
            end = start + cigar_softclip_lengths[0]

        # Get the sequence of the (single) soft-clipped section:
        softclipped_seq = read.seq[dist_to_softclipped:dist_to_softclipped+cigar_softclip_lengths[0]]
        softclipped_quals = read.query_qualities[dist_to_softclipped:dist_to_softclipped+cigar_softclip_lengths[0]]
        softclipped_seq_masked = ""
        for idx in range(len(softclipped_seq)):
            if softclipped_quals[idx] < 30:
                softclipped_seq_masked += "N"
            else:
                softclipped_seq_masked += softclipped_seq[idx]

        origin_chrom = chrom_int2str(read.reference_id)
        if dist_to_softclipped == 0:
            origin_start = read.pos - cigar_softclip_lengths[0]
            origin_end = read.pos
        else:
            origin_start = read.pos + dist_to_softclipped
            origin_end = read.pos + dist_to_softclipped + cigar_softclip_lengths[0]

        # FIXME: I really need to make sure the exraction of the softclip sequence/above
        # lengths is correct. It's tricky due to the complicated cigar string structure.
        # Currently, just preventing crashes by not creating an object when the sequence
        # is empty (I'm in a hurry to get some initial results!):
        if len(softclipped_seq_masked) > 0:
            soft_clippings = [SoftClipping(softclipped_seq_masked, read,
                                           origin_chrom, origin_start, origin_end)]

    return soft_clippings


def group_softclip_coords(terminus_reads):
    """
    Group the soft-clipping coordinates according to their start/end points.

    :param soft_clippings: List of terminus reads
    :return: A list of soft-clippings grouped according to their genomic positions.
    """

    soft_clippings = functools.reduce(lambda reads1, reads2: reads1 + reads2,
                                      [get_soft_clippings(read) for read in terminus_reads])

    unique_starts = set([soft_clipping._origin_start for soft_clipping in soft_clippings])
    unique_ends = set([soft_clipping._origin_end for soft_clipping in soft_clippings])

    # NOTE: Not computationally efficient but small numbers => doesn't matter
    start_to_soft_clippings = defaultdict(list)
    for start in unique_starts:
        for soft_clipping in soft_clippings:
            if soft_clipping._origin_start == start:
                start_to_soft_clippings[start].append(soft_clipping)

    end_to_soft_clippings = defaultdict(list)
    for end in unique_ends:
        for soft_clipping in soft_clippings:
            if soft_clipping._origin_end == end:
                end_to_soft_clippings[end].append(soft_clipping)

        end_to_soft_clippings[end].append(soft_clipping)

    # Assign each soft-clipping to a single start or end group (choose largest
    # group):
    start_to_soft_clippings_final = defaultdict(list)
    end_to_soft_clippings_final = defaultdict(list)
    for soft_clipping in soft_clippings:
        start_group = start_to_soft_clippings[soft_clipping._origin_start]
        end_group = end_to_soft_clippings[soft_clipping._origin_end]
        if len(start_group) >= len(end_group):
            start_to_soft_clippings_final[soft_clipping._origin_start].append(soft_clipping)
        else:
            end_to_soft_clippings_final[soft_clipping._origin_end].append(soft_clipping)

    return list(start_to_soft_clippings_final.values()) + list(end_to_soft_clippings_final.values())


def get_softclip_bounds(soft_clippings):
    if len(soft_clippings) > 0:
        # FIXME: Currently, I will just assume that all soft-clippings in the input
        # are overlapping. I should eventually change this to accommodate distinct
        # soft clipping supports.
        chrom = soft_clippings[0].get_chrom()
        start = min(list(map(lambda soft_clipping: soft_clipping.get_start(), soft_clippings)))
        end = max(list(map(lambda soft_clipping: soft_clipping.get_end(), soft_clippings)))

        return [(chrom, start, end)]
    else:
        return []


def test_matched_soft_clipping(chrom, start, end, reads2, fasta_filename):
    '''Examine soft-clipping of all reads in reads2, to see if any of the
    soft-clipped coordinates are consistent with the specified coordinates.'''

    # FIXME: Quick hacky solution:
    # If there are many no-call (phred score == 2) reads in all the reads, then
    # simply return that the read has no useable matching softclipping
    # sequence:
    poor_qual_bp = len(list(filter(lambda qual: qual == 2, functools.reduce(lambda l1, l2: l1 + l2, list(map(lambda read: read.query_qualities, reads2))))))
    total_bp = len(functools.reduce(lambda l1, l2: l1 + l2, list(map(lambda read: read.query_qualities, reads2))))
    if poor_qual_bp/float(total_bp) > 0.2:
        return (defaultdict(list), defaultdict(list))

    # Find all soft-clippings from reads2 that reside in the region spanned
    # by reads in reads1:
    consistent_soft_clippings = []
    for read in reads2:
        curr_soft_clippings = get_soft_clippings(read)
        for soft_clipping in curr_soft_clippings:
            if soft_clipping.is_in(chrom, start, end, fasta_filename):
                consistent_soft_clippings.append(soft_clipping)

    # Derive a set of spatially-distinct, consensus soft-clippings
    # from the above set:
    consensus_soft_clippings = get_softclip_bounds(consistent_soft_clippings)
    return consensus_soft_clippings


class GenomicEvent:
    def __init__(self, terminusA_reads, terminusB_reads):
        # Set the event termini with terminus 1 as the first terminus
        # and terminus 2 as the second as determined by chromosomal
        # numbering/position:
        terminusA_chrom = terminusA_reads[0].rname
        terminusA_pos = min(list(map(lambda read: read.pos, terminusA_reads)))

        terminusB_chrom = terminusB_reads[0].rname
        terminusB_pos = min(list(map(lambda read: read.pos, terminusB_reads)))

        if terminusA_chrom < terminusB_chrom:
            self._terminus1_reads = terminusA_reads
            self._terminus2_reads = terminusB_reads
        elif terminusA_chrom > terminusB_chrom:
            self._terminus1_reads = terminusB_reads
            self._terminus2_reads = terminusA_reads
        else:
            if terminusA_pos < terminusB_pos:
                self._terminus1_reads = terminusA_reads
                self._terminus2_reads = terminusB_reads
            else:
                self._terminus1_reads = terminusB_reads
                self._terminus2_reads = terminusA_reads

        self._matched_softclips_t1 = []
        self._matched_softclips_t2 = []

    def get_t1_mapqual(self):
        mapqual = max(list(map(lambda read: read.mapq, self._terminus1_reads)))
        return mapqual

    def get_t2_mapqual(self):
        return max(list(map(lambda read: read.mapq, self._terminus2_reads)))

    def has_soft_clip_support(self):
        return len(self._matched_softclips_t1) > 0 or len(self._matched_softclips_t2) > 0

    def count_unique_softclip_regions(self):
        self._softclip_groups_t1 = group_softclip_coords(self._terminus1_reads)
        self._softclip_groups_t2 = group_softclip_coords(self._terminus2_reads)

    def has_scattered_soft_clip_regions(self):
        self.count_unique_softclip_regions()

        # Change this once I decide on suitable filtering parameters
        return False

    def get_terminus1_span(self, extension_length=0):
        return self._get_reads_span(self._terminus1_reads, extension_length=extension_length)

    def get_terminus2_span(self, extension_length=0):
        return self._get_reads_span(self._terminus2_reads, extension_length=extension_length)

    def _get_reads_span(self, reads, extension_length=0):
        reads_chrom = reads[0].rname
        reads_strand = "+"
        if reads[0].flag & 16:
            reads_strand = "-"

        # TEMPORARY SANITY CHECK:
        for read in reads:
            curr_strand = "+"
            if read.flag & 16:
                curr_strand = "-"
            assert read.rname == reads_chrom
            assert curr_strand == reads_strand

        reads_start = min(list(map(lambda read: read.pos, reads)))
        reads_end = max(list(map(lambda read: read.pos+read.qlen, reads)))

        # FIXME: Hard-coding an extension value here. Make this a command-line
        # parameter instead:
        start = reads_start - extension_length
        end = reads_end + extension_length
        chrom = chrom_int2str(reads_chrom)
        return (chrom, start, end, reads_strand)

    def test_soft_clipping(self, fasta_filename):
        '''Look at the soft-clipping for reads in terminus 1, and see if they match the
    	position of reads in terminus 2. Then, do the reverse. Record the result'''

        # NOTE: Extend the event termini for the purpose of testing for
        # soft-clipping support:
        terminus1_span = self.get_terminus1_span(extension_length=100)
        terminus2_span = self.get_terminus2_span(extension_length=100)

        # XXX CONTINUE HERE: these are now two dictionaries; use them to store the relevant stats for the events,
        # or in the final filtering stage.
        self._matched_softclips_t1 = test_matched_soft_clipping(\
            terminus1_span[0], terminus1_span[1], terminus1_span[2], \
            self._terminus2_reads, fasta_filename)

        self._matched_softclips_t2 = test_matched_soft_clipping(\
            terminus2_span[0], terminus2_span[1], terminus2_span[2], \
            self._terminus1_reads, fasta_filename)

    def get_t1_depth(self):
        return len(self._terminus1_reads)

    def get_t2_depth(self):
        return len(self._terminus2_reads)

    def get_gtf(self):
        t1_chrom = self._terminus1_reads[0].rname
        t1_first_read_start = min(list(map(lambda read: read.pos, list(self._terminus1_reads))))
        t1_last_read_end = max(list(map(lambda read: read.pos+read.qlen, list(self._terminus1_reads))))

        t2_chrom = self._terminus2_reads[0].rname
        t2_first_read_start = min(list(map(lambda read: read.pos, list(self._terminus2_reads))))
        t2_last_read_end = max(list(map(lambda read: read.pos+read.qlen, list(self._terminus2_reads))))

        num_t1_softclip_groups = len(self._softclip_groups_t1)
        num_t1_softclips_in_groups = functools.reduce(lambda len1, len2: len1 + len2,
                                                      [len(group) for group in self._softclip_groups_t1])
        num_t2_softclip_groups = len(self._softclip_groups_t2)
        num_t2_softclips_in_groups = functools.reduce(lambda len1, len2: len1 + len2,
                                                      [len(group) for group in self._softclip_groups_t2])

        # Temporary code for printing genomic events in gff format...
        t1_gtf = "%s\tSV_event\texon\t%d\t%d\t%d\t.\t.\tgene_id \"%s\"; transcript_id \"%s\";#%d %d\n" % \
          (chrom_int2str(t1_chrom), t1_first_read_start, t1_last_read_end, self.get_t1_depth(),
           self.get_event_name(), self.get_event_name(), num_t1_softclip_groups, num_t1_softclips_in_groups)

        t2_gtf = "%s\tSV_event\texon\t%d\t%d\t%d\t.\t.\tgene_id \"%s\"; transcript_id \"%s\";#%d %d\n" % \
          (chrom_int2str(t2_chrom), t2_first_read_start, t2_last_read_end, self.get_t2_depth(),
           self.get_event_name(), self.get_event_name(), num_t2_softclip_groups, num_t2_softclips_in_groups)

        return t1_gtf + t2_gtf + self.get_softclip_gtf()

    # Hack: Generate an event-name based on the coordinates:
    def get_event_name(self):
        t1_chrom = self._terminus1_reads[0].rname
        t1_first_read_start = min(list(map(lambda read: read.pos, list(self._terminus1_reads))))
        t1_last_read_end = max(list(map(lambda read: read.pos+read.qlen, list(self._terminus1_reads))))

        t2_chrom = self._terminus2_reads[0].rname
        t2_first_read_start = min(list(map(lambda read: read.pos, list(self._terminus2_reads))))
        t2_last_read_end = max(list(map(lambda read: read.pos+read.qlen, list(self._terminus2_reads))))

        return chrom_int2str(t1_chrom) + ":" + str(t1_first_read_start) + "-" + str(t1_last_read_end) + \
          "," + chrom_int2str(t2_chrom) + ":" + str(t2_first_read_start) + "-" + str(t2_last_read_end)

    def get_softclip_gtf(self):
        t1_gtf = ""
        if self._matched_softclips_t1 != []:
            chrom = self._matched_softclips_t1[0][0]
            start = self._matched_softclips_t1[0][1]
            end = self._matched_softclips_t1[0][2]
            t1_gtf = "%s\tSV_event\tCDS\t%d\t%d\t%d\t.\t.\tgene_id \"%s\"; transcript_id \"%s\";\n" % \
            (chrom, start, end, len(self._matched_softclips_t1),
            self.get_event_name(), self.get_event_name())

        t2_gtf = ""
        if self._matched_softclips_t2 != []:
            chrom = self._matched_softclips_t2[0][0]
            start = self._matched_softclips_t2[0][1]
            end = self._matched_softclips_t2[0][2]
            t2_gtf = "%s\tSV_event\tCDS\t%d\t%d\t%d\t.\t.\tgene_id \"%s\"; transcript_id \"%s\";\n" % \
            (chrom, start, end, len(self._matched_softclips_t2),
            self.get_event_name(), self.get_event_name())

        return t1_gtf + t2_gtf


def define_events(clusters, read2cluster, read2mate):
    '''Define a set of putative events by inspecting pairings of
    clusters.'''

    # Keep track of reads that have been assigned to an event:
    reads_already_assigned = set()

    events = set()

    for cluster in clusters:
        for paired_cluster in cluster.get_paired_clusters():
            curr_cluster1_reads = list(filter(lambda read: read2cluster[read2mate[read]] == paired_cluster,
                                         cluster.get_reads()))
            curr_cluster2_reads = list(filter(lambda read: read2cluster[read2mate[read]] == cluster,
                                         paired_cluster.get_reads()))

            clust1_reads_not_assigned = functools.reduce(lambda bool1, bool2: bool1 and bool2,
                [not read in reads_already_assigned for read in curr_cluster1_reads])
            clust2_reads_not_assigned = functools.reduce(lambda bool1, bool2: bool1 and bool2,
                [not read in reads_already_assigned for read in curr_cluster2_reads])

            if clust1_reads_not_assigned and clust2_reads_not_assigned:
               # None of these reads have been assigned to an event. Assign them to a new
               # event:
               curr_event = GenomicEvent(curr_cluster1_reads, curr_cluster2_reads)
               events.add(curr_event)
               reads_already_assigned = reads_already_assigned.union(curr_cluster1_reads).union(curr_cluster2_reads)

    return events


def call_events(filtered_reads, fasta_filename):
    '''Call events on the input reads. The reads must be sorted by chromosome
    and then by position.'''

    # Detect overlapping read clusters...
    clusters = detect_clusters(filtered_reads)

    # Pair up the clusters:
    (read2cluster, read2mate) = pair_clusters(clusters)

    # Define the events, using those pairings:
    putative_events = define_events(clusters, read2cluster, read2mate)

    # Mark each putative event with soft-clipping information of the reads contained in it:
    for event in putative_events:
        event.test_soft_clipping(fasta_filename)

    return putative_events


class ReadCluster:
    def __init__(self):
        self._reads = set()
        self._paired_clusters = set()

    def add_read(self, read):
        self._reads.add(read)

    def get_reads(self):
        return self._reads

    def get_paired_clusters(self):
        return self._paired_clusters

    def set_pairings(self, read2mate, read2cluster):
        '''Record all cluster pairings for this cluster.'''

        for read in self._reads:
            mate = read2mate[read]
            mate_cluster = read2cluster[mate]
            self._paired_clusters.add(mate_cluster)


    def to_string(self):
        read_list = list(self._reads)
        first_read_chrom = read_list[0].rname

        # TEMPORARY SANITY CHECK - implementing quickly now and will write a suite
        # of unit tests later when time abides:
        for read in read_list:
            assert read.rname == first_read_chrom

        first_read_start = min(list(map(lambda read: read.pos, read_list)))

        last_read_end = max(list(map(lambda read: read.pos+read.qlen, read_list)))

        return "%s\t%d\t%d" % (chrom_int2str(first_read_chrom), first_read_start, last_read_end)


def cmp_coords(coords1, coords2):
    if coords1[0] < coords2[0]:
        return -1
    elif coords1[0] > coords2[0]:
        return 1
    else:
        if coords1[1] < coords2[1]:
            return -1
        elif coords1[1] > coords2[1]:
            return 1
        else:
            if coords1[2] < coords2[2]:
                return -1
            elif coords1[2] > coords2[2]:
                return 1
            else:
                return 0


class CoordKey(object):
    def __init__(self, obj, *args):
        self.obj = obj

    def __lt__(self, other):
        return cmp_coords(self.obj, other.obj) < 0

    def __gt__(self, other):
        return cmp_coords(self.obj, other.obj) > 0

    def __eq__(self, other):
        return cmp_coords(self.obj, other.obj) == 0

    def __le__(self, other):
        return cmp_coords(self.obj, other.obj) <= 0

    def __ge__(self, other):
        return cmp_coords(self.obj, other.obj) >= 0

    def __ne__(self, other):
        return cmp_coords(self.obj, other.obj) != 0


def filter_on_shared_termini(events):
    # Get dictionary of events for each span (should almost always be one
    # event per span, as a span is a (chrom, start, end, strand) tuple):
    terminus_span_to_events = {}
    for event in events:
        terminus1_span = event.get_terminus1_span()
        if not terminus1_span in terminus_span_to_events:
            terminus_span_to_events[terminus1_span] = []
        terminus_span_to_events[terminus1_span].append(event)

        terminus2_span = event.get_terminus2_span()
        if not terminus2_span in terminus_span_to_events:
            terminus_span_to_events[terminus2_span] = []
        terminus_span_to_events[terminus2_span].append(event)

    all_termini = list(terminus_span_to_events.keys())
    all_termini.sort(key=CoordKey)

    # Scan over the events in ascending order, and find any sequential
    # overlapping termini. Record all events that have one or more
    # overlapping termini:
    events_with_overlaps = set()
    plus_strand_termini = [tup for tup in all_termini if tup[-1] == "+"]
    minus_strand_termini = [tup for tup in all_termini if tup[-1] == "-"]
    events_with_plus_strand_overlaps = \
        find_terminus_overlaps(plus_strand_termini, terminus_span_to_events)
    events_with_minus_strand_overlaps = \
        find_terminus_overlaps(minus_strand_termini, terminus_span_to_events)

    events_with_overlaps = \
        events_with_plus_strand_overlaps.union(events_with_minus_strand_overlaps)

    filtered_events = list(filter(lambda event: not event in events_with_overlaps,
                             events))

    return filtered_events


def find_terminus_overlaps(termini, terminus_span_to_events):
    events_with_overlaps = set()
    terminus_idx = 0
    prev_chrom = None
    prev_start = -1
    prev_end = -1
    prev_terminus = None
    while terminus_idx < len(termini):
        curr_terminus = termini[terminus_idx]
        if curr_terminus[0] != prev_chrom or curr_terminus[1] > prev_end:
            # New terminus cluster:
            prev_chrom = curr_terminus[0]
            prev_start = curr_terminus[1]
        else:
            # This terminus overlaps the existing cluster of termini:
            for event in terminus_span_to_events[curr_terminus]:
                events_with_overlaps.add(event)
            for event in terminus_span_to_events[prev_terminus]:
                events_with_overlaps.add(event)
        prev_end = curr_terminus[2]
        prev_terminus = curr_terminus
        terminus_idx += 1

    return events_with_overlaps


def detect_clusters(reads):
    clusters = set()

    plus_strand_reads = list(filter(lambda read: not(read.flag & 16), reads))
    minus_strand_reads = list(filter(lambda read: read.flag & 16, reads))
    detect_clusters_single_strand(clusters, plus_strand_reads)
    detect_clusters_single_strand(clusters, minus_strand_reads)

    return clusters


def detect_clusters_single_strand(clusters, reads, prev_read_extn=0):
    prev_read_chrom = None
    prev_read_end_pos = 0
    curr_cluster = None

    for read in reads:
        # Detect if the read is on a different chromosome, or if it does not
        # overlap the previous read
        if read.rname != prev_read_chrom or \
            read.pos > (prev_read_end_pos + prev_read_extn):
            # Start a new cluster:
            curr_cluster = ReadCluster()
            clusters.add(curr_cluster)

        # Add the read to the current cluster:
        curr_cluster.add_read(read)

        prev_read_chrom = read.rname

        # NOTE: Using qlen as the length of the read for the purpose of determining
        # cluster overlaps. I suspect this will perform adequately, but may need to
        # adjust it to instead use matching sequence length instead at some point:
        prev_read_end_pos = read.pos + read.qlen
