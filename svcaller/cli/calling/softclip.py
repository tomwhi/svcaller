# coding=utf-8
import pdb
import pybedtools
from past import autotranslate
autotranslate(['swalign'])
import swalign

def revcomp(seq):
    return "".join([{"A":"T","C":"G","G":"C","T":"A","N":"N"}[lett] for lett in [seq[i] for i in [len(seq) - idx - 1 for idx in range(len(seq))]]])


def getMatchSize(tup):
    if tup[1] == "M":
        return tup[0]
    else:
        return 0


def getAlign(ref_seq, query_seq):
    match = 1
    mismatch = -5
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring)  # you can also choose gap penalties, etc...
    alignment = sw.align(ref_seq, query_seq)
    return alignment


def processAlign(alignment):
    matchSizes = list(map(lambda tup: getMatchSize(tup), alignment.cigar))
    # FIXME: Another nasty quick fix:
    if len(alignment.cigar) == 0:
        return (0, -1)
    try:
        maxMatch = max(matchSizes)
    except Exception:
        pdb.set_trace()
        dummy = 1
    matchingIdxs = []
    for idx in range(len(matchSizes)):
        if matchSizes[idx] == maxMatch:
            matchingIdxs.append(idx)
    firstLargestMatchIdx = matchingIdxs[0]
    displacementFromRefStart = alignment.r_pos
    for idx in range(firstLargestMatchIdx):
        tup = alignment.cigar[idx]
        if tup[1] == "M":
            displacementFromRefStart += tup[0]
        elif tup[1] == "D":
            displacementFromRefStart += tup[0]
    return (maxMatch, displacementFromRefStart)


def getRefMatch(eventChrom, eventStart, eventEnd, softclippedSeq, fasta_filename):
    regs = [pybedtools.Interval(eventChrom, eventStart, eventEnd)]
    tmp_bedtool = pybedtools.BedTool(regs)
    refseq = tmp_bedtool.seq((eventChrom, eventStart, eventEnd), fasta_filename)
    fwalign = getAlign(refseq, softclippedSeq)
    fw_match = processAlign(fwalign)

    refseq_revcomp = revcomp(refseq)
    revalign = getAlign(refseq_revcomp, softclippedSeq)
    rev_match = processAlign(revalign)

    if fw_match[0] > rev_match[0]:
        return (fw_match[0], fw_match[1], "+")
    else:
        return (rev_match[0], rev_match[1], "-")


def getRefMatchPos(eventChrom, eventStart, eventEnd, softclippedSeq, fasta_filename):
#    if eventEnd == 66944558 and softclippedSeq == "ATAAATAAATAAATAAATAAATACGTACATACATACACACATACATACA":
#        pdb.set_trace()
#        dummy = 1
    (maxMatch, displacementFromRefStart, strand) = getRefMatch(eventChrom, eventStart, eventEnd, softclippedSeq, fasta_filename)
    if strand == "+":
        matchRefPosStart = eventStart + displacementFromRefStart
        matchRefPosEnd = matchRefPosStart + maxMatch
    else:
        matchRefPosEnd = eventEnd - displacementFromRefStart
        matchRefPosStart = eventEnd - displacementFromRefStart - maxMatch
    if maxMatch > 15: # Hard-coded match threshold:
        return (eventChrom, matchRefPosStart+1, matchRefPosEnd)
    else:
        return None
