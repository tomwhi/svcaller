# coding=utf-8
from datetime import datetime
import pdb, sys
import pysam
import softclip

DEL = "DEL"
INV = "INV"
TRA = "TRA"
DUP = "DUP"


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

#        if read.qname == "HISEQ:226:HK2GCBCXX:2:2216:5917:85768":
#            pdb.set_trace()
#            dummy = 1

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
    # Select reads facing away from each other, on the same chromosome:
    if (read.rname == read.rnext):
        if not(read.flag & 16) and (read.flag & 32):
            if read.tlen < 0:
                filtered_reads.append(read)
        elif (read.flag & 16) and not(read.flag & 32):
            if read.tlen > 0:
                filtered_reads.append(read)


def inv_filt(filtered_reads, read):
    # Select reads facing in the same direction, on the same chromosome:
    if (read.rname == read.rnext):
        if not(read.flag & 16) and not(read.flag & 32):
            filtered_reads.append(read)
        elif (read.flag & 16) and (read.flag & 32):
            filtered_reads.append(read)


def tra_filt(filtered_reads, read):
    # Select reads facing in the same direction, on the same chromosome:
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
                if not reads_nearby_dict.has_key(curr_read.qname):
                    reads_nearby_dict[curr_read.qname] = curr_read

            reads_nearby = reads_nearby_dict.values()

            paired_matches = 0
            for read_nearby in reads_nearby:
                read_nearby_pair_chrom = read_nearby.rnext
                if read_nearby_pair_chrom == other_read_chrom and \
                   other_read_pos+half_read_len-1000 < read_nearby.pnext+half_read_len and \
                   other_read_pos+half_read_len+1000 > read_nearby.pnext+half_read_len:
                    paired_matches += 1

#            if read.qname == "HWI-D00410:258:HTFWFBCXX:1:1210:10002:56665":
#                pdb.set_trace()
#                dummy2 = 1

            if paired_matches >= 2:
                filtered_reads.append(read)

        except Exception, e:
            #print >> sys.stderr, e
            pass

        if idx % 1000 == 0:
            print >> sys.stderr, idx
        idx += 1

    # Finally, filter to only retain reads where both reads in the pair are
    # included in this list:
    return paired_filt(filtered_reads)


def paired_filt(reads):
    '''Filter the specified reads, to only retain those where both read-pairs are present
    in the input list.'''

    readname2reads = get_readname2reads(reads)
    return filter(lambda read: len(readname2reads[read.qname]) == 2, reads)


def get_readname2reads(reads):
    readname2reads = {}
    for read in reads:
        if not readname2reads.has_key(read.qname):
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
        other_reads = filter(lambda curr_read: curr_read != read, all_reads_with_this_name)
        assert len(other_reads) == 1
        read2mate[read] = other_reads[0]

    # Update each cluster, to indicate a pairing to each other cluster that the mate-pair
    # belongs to, for all mate-pairs of reads in this cluster (hope that makes sense!):
    for cluster in clusters:
        cluster.set_pairings(read2mate, read2cluster)

    # FIXME: This is getting messy:
    return (read2cluster, read2mate)


class SoftClipping:
    def __init__(self, sequence):
        self._chrom = None
        self._start = None
        self._end = None
        self._seq = sequence

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
            print >> sys.stderr, pos_tup

        # Sequence alignment-based approach:
        return pos_tup != None


def get_soft_clippings(read):
    '''Extract all soft-clipping secondary mapping coordinates from this
    read.'''

    cigar_softclip_lengths = []
    dist_to_softclipped = 0
    if read.cigartuples != None:
        cigar_softclip_lengths = map(lambda tup: tup[1], filter(lambda tup: tup[0] == 4, read.cigartuples))

        # Calculate displacement of the softclipped region from the start of the *sequence*:
        if len(cigar_softclip_lengths) > 0:
            # This is completely retarded but oh well.
            idxs = range(len(read.cigartuples))
            arr = filter(lambda idx: map(lambda tup: tup[0] == 4, read.cigartuples)[idx], idxs)
            first_softclip_idx = arr[0]
            before_softclip_tups = read.cigartuples[:first_softclip_idx]
            dist_to_softclipped = sum(map(lambda tup: tup[1], before_softclip_tups))

    softclip_starts = map(lambda tup: tuple(tup[1].split(",")[:2]), \
        filter(lambda tag: tag[0] == "SA", read.get_tags()))

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

        soft_clippings = [SoftClipping(softclipped_seq)]

    return soft_clippings


def condense_soft_clippings(soft_clippings):
    if len(soft_clippings) > 0:
        # FIXME: Currently, I will just assume that all soft-clippings in the input
        # are overlapping. I should eventually change this to accommodate distinct
        # soft clipping supports.
        chrom = soft_clippings[0].get_chrom()
        start = min(map(lambda soft_clipping: soft_clipping.get_start(), soft_clippings))
        end = max(map(lambda soft_clipping: soft_clipping.get_end(), soft_clippings))

        return [(chrom, start, end)]
    else:
        return []


def test_matched_soft_clipping(reads1, reads2, fasta_filename):
    '''Examine soft-clipping of all reads in reads2, to see if any of the
    soft-clipped coordinates are consistent with coordinates spanned by
    all reads in reads1.'''

    reads1_chrom = reads1[0].rname

    # TEMPORARY SANITY CHECK:
    for read in reads1:
        assert read.rname == reads1_chrom

    reads1_start = min(map(lambda read: read.pos, reads1))
    reads1_end = max(map(lambda read: read.pos+read.qlen, reads1))

    # Find all soft-clippings from reads2 that reside in the region spanned
    # by reads in reads1:
    consistent_soft_clippings = []
    for read in reads2:
        curr_soft_clippings = get_soft_clippings(read)
        for soft_clipping in curr_soft_clippings:
            if soft_clipping.is_in(chrom_int2str(reads1_chrom), reads1_start, reads1_end, fasta_filename):
                consistent_soft_clippings.append(soft_clipping)

    # Derive a set of spatially-distinct, consensus soft-clippings
    # from the above set:
    consensus_soft_clippings = condense_soft_clippings(consistent_soft_clippings)
    return consensus_soft_clippings


class GenomicEvent:
    def __init__(self, terminus1_reads, terminus2_reads):
        self._terminus1_reads = terminus1_reads
        self._terminus2_reads = terminus2_reads

        self._matched_softclips_t1 = None
        self._matched_softclips_t2 = None

    def has_soft_clip_support(self):
        return len(self._matched_softclips_t1) > 0 or len(self._matched_softclips_t1) > 0

    def test_soft_clipping(self, fasta_filename):
    	'''Look at the soft-clipping for reads in terminus 1, and see if they match the
    	position of reads in terminus 2. Then, do the reverse. Record the result'''

        self._matched_softclips_t1 = test_matched_soft_clipping(\
            self._terminus1_reads, self._terminus2_reads, fasta_filename)

        self._matched_softclips_t2 = test_matched_soft_clipping(\
            self._terminus2_reads, self._terminus1_reads, fasta_filename)

    def to_string(self):
        t1_chrom = self._terminus1_reads[0].rname
        t1_first_read_start = min(map(lambda read: read.pos, list(self._terminus1_reads)))
        t1_last_read_end = max(map(lambda read: read.pos+read.qlen, list(self._terminus1_reads)))

        t2_chrom = self._terminus2_reads[0].rname
        t2_first_read_start = min(map(lambda read: read.pos, list(self._terminus2_reads)))
        t2_last_read_end = max(map(lambda read: read.pos+read.qlen, list(self._terminus2_reads)))

        # Temporary code for printing genomic events in gff format...
        # Hack: Generate an event-name based on the coordinates:
        event_name = chrom_int2str(t1_chrom) + "_" + str(t1_first_read_start) + "_" + str(t2_last_read_end)

        return "chr%s\ta\tb\t%d\t%d\t500\t+\t.\t%s\n" % (chrom_int2str(t1_chrom), t1_first_read_start, t1_last_read_end, event_name) + \
            "chr%s\ta\tb\t%d\t%d\t500\t+\t.\t%s" % (chrom_int2str(t2_chrom), t2_first_read_start, t2_last_read_end, event_name)

    def get_softclip_strings(self):
        t1_chrom = self._terminus1_reads[0].rname
        t1_first_read_start = min(map(lambda read: read.pos, list(self._terminus1_reads)))
        t1_last_read_end = max(map(lambda read: read.pos+read.qlen, list(self._terminus1_reads)))

        t2_chrom = self._terminus2_reads[0].rname
        t2_first_read_start = min(map(lambda read: read.pos, list(self._terminus2_reads)))
        t2_last_read_end = max(map(lambda read: read.pos+read.qlen, list(self._terminus2_reads)))

        event_name = chrom_int2str(t1_chrom) + "_" + str(t1_first_read_start) + "_" + str(t2_last_read_end)

        out_string = ""
        if self._matched_softclips_t1 != None:
            chrom = self._matched_softclips_t1[0][0]
            start = self._matched_softclips_t1[0][1]
            end = self._matched_softclips_t1[0][2]
            out_string += "chr%s\ta\tb\t%d\t%d\t500\t+\t.\t%s\n" % (chrom, start, end, event_name)

        if self._matched_soft_clips_t2 != None:
            chrom = self._matched_softclips_t2[0][0]
            start = self._matched_softclips_t2[0][1]
            end = self._matched_softclips_t2[0][2]
            out_string += "chr%s\ta\tb\t%d\t%d\t500\t+\t.\t%s\n" % (chrom, start, end, event_name)

        return out_string


def define_events(clusters, read2cluster, read2mate):
    '''Define a set of putative events by inspecting pairings of
    clusters.'''

    # Keep track of reads that have been assigned to an event:
    reads_already_assigned = set()

    events = set()

    for cluster in clusters:
        for paired_cluster in cluster.get_paired_clusters():
            curr_cluster1_reads = filter(lambda read: read2cluster[read2mate[read]] == paired_cluster,
                                         cluster.get_reads())
            curr_cluster2_reads = filter(lambda read: read2cluster[read2mate[read]] == cluster,
                                         paired_cluster.get_reads())

            if reduce(lambda bool1, bool2: bool1 and bool2,
                      [not read in reads_already_assigned for read in curr_cluster1_reads]) and \
               reduce(lambda bool1, bool2: bool1 and bool2,
                      [not read in reads_already_assigned for read in curr_cluster2_reads]):
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

    # TMP: PRINTING EVENTS TO A GFFFILE:
#    header = '''browser position chrX:66761874-66952461
#browser hide all
#track name=genomic_events\tdescription="Genomic_events"\tvisibility=2'''
#    print header
#
#    for event in list(putative_events):
#        print event.to_string()


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

        first_read_start = min(map(lambda read: read.pos, read_list))

        last_read_end = max(map(lambda read: read.pos+read.qlen, read_list))

        return "%s\t%d\t%d" % (chrom_int2str(first_read_chrom), first_read_start, last_read_end)


def detect_clusters(reads):
    clusters = set()

    plus_strand_reads = filter(lambda read: not(read.flag & 16), reads)
    minus_strand_reads = filter(lambda read: read.flag & 16, reads)
    detect_clusters_single_strand(clusters, plus_strand_reads)
    detect_clusters_single_strand(clusters, minus_strand_reads)

    return clusters


def detect_clusters_single_strand(clusters, reads):
    prev_read_chrom = None
    prev_read_end_pos = 0
    curr_cluster = None

    for read in reads:
        # Detect if the read is on a different chromosome, or if it does not
        # overlap the previous read
        if read.rname != prev_read_chrom or read.pos > prev_read_end_pos:
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
































