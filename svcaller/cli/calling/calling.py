# coding=utf-8
from datetime import datetime
import pdb, sys
import pysam

DEL = "DEL"
INV = "INV"
TRA = "TRA"
DUP = "DUP"

def call_event(input_bam, output_name, event_type):
    samfile = pysam.AlignmentFile(input_bam, "rb")
    if event_type == DEL:
        event_filtered_reads = del_filter(samfile)
    elif event_type == INV:
        event_filtered_reads = inv_filter(samfile)
    elif event_type == TRA:
        event_filtered_reads = tra_filter(samfile)
    elif event_type == DUP:
        event_filtered_reads = dup_filter(samfile)
    else:
        logging.error("Invalid event_type: " + event_type)

    cluster_filtered_reads = cluster_filter(event_filtered_reads, samfile)
    pdb.set_trace()
    with pysam.AlignmentFile(output_name + ".bam", "wb", header=samfile.header) as outf:
        for read in cluster_filtered_reads:
            outf.write(read)
    #XXX = call_events(second_filtered_bam)
    #XXX = check_softclipping(second_filtered_bam)
    #generate_outputs(XXX, XXX)


def del_filter(samfile):
    filtered_reads = []
    for read in samfile:
        # Select reads facing each other, on the same chromosome, with a gap
        # larger than the specified distance threshold:
        if (read.rname == read.rnext):
            if not(read.flag & 16) and (read.flag & 32):
                if read.tlen > 1000:
                    filtered_reads.append(read)
            elif (read.flag & 16) and not(read.flag & 32):
                if read.tlen < -1000:
                    filtered_reads.append(read)

    return filtered_reads


def cluster_filter(input_reads, samfile, chrom_of_interest=22):
    '''Filter the input reads to only retain those read-pairs where there are other read-pairs with
    similar coordinates. This is an ugly implementation, needs refactoring.

    chrom_of_interest is a chromosome that is expected to contain
    more reads than the other chromosomes.'''
    idx = 0
    filtered_reads = []
    for read in input_reads:
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

        base_read_chrom_str = str(base_read_chrom + 1)
        if base_read_chrom_str == "23":
            base_read_chrom_str = "X"
        elif base_read_chrom_str == "24":
            base_read_chrom_str = "Y"

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
            print >> sys.stderr, e
            pass

        if idx % 10 == 0:
            print >> sys.stderr, idx
        idx += 1

    return filtered_reads


comment = '''
import pdb, sys
import pysam

samfile = pysam.AlignmentFile(sys.argv[1], "rb")

allReads = list(samfile)
filteredReads = []


outFile = open("dummy.sam", 'w')
outSam = pysam.AlignmentFile(outFile, "w", template=samfile)
out2 = open(sys.argv[2] + ".sam", 'w')
for read in filteredReads:
  print >> out2, read.tostring(outSam)

out2.close()
'''