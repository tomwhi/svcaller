# coding=utf-8
from datetime import datetime
import pdb, sys
import pysam

DEL = "DEL"
INV = "INV"
TRA = "TRA"
DUP = "DUP"


#def call_event(input_bam, output_name, event_type):
#
#    cluster_filtered_reads = cluster_filter(event_filtered_reads, samfile)
#    pdb.set_trace()
#    #XXX = call_events(second_filtered_bam)
#    #XXX = check_softclipping(second_filtered_bam)
#    #generate_outputs(XXX, XXX)


def event_filt(read_iter, event_type, flag_filter=256+1024+2048):
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
            #print >> sys.stderr, e
            pass

        if idx % 1000 == 0:
            print >> sys.stderr, idx
        idx += 1

    return filtered_reads