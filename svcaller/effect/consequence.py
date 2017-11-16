from pybedtools import BedTool
from enum import Enum


class SvEffect(Enum):
    NO_OVERLAP = "NO_OVERLAP"
    OVERLAP_WITH_EFFECT = "OVERLAP_WITH_EFFECT"
    OVERLAP_UNKNOWN_EFFECT = "OVERLAP_UNKNOWN_EFFECT"
    GENE_FUSION = "GENE_FUSION"


class GeneClass(Enum):
    TUMOUR_SUPPRESSOR = "TUMOUR_SUPPRESSOR"
    AR = "AR"
    FUSION_CANDIDATE = "FUSION_CANDIDATE"


def predict_effects(svs_filename, ts_filename, ar_filename, fusion_filename):
    """
    Predict the consequence of the specified structural variants on the specified
    tumour suppressors, androgen receptor, and gene fusion candidate.

    :param svs_filename: Location of bed file specifying the structural variant coordinates.
    :param ts_filename: Location of bed file specifying the tumour supressor gene region coords.
    :param ar_filename: Location of bed file specifying Androgren Receptor gene region coords.
    :param fusion_filename: Location of bed file specifying two broad gene fusion region.
    :param effects_file: Output filehandle for writing JSON data stating predicted effects
    :return: 
    """

    gene_classes = [enum.value for enum in list(GeneClass)]

    svs_bed = BedTool(svs_filename)

    gene_region_beds = [BedTool(filename) for filename in
                        [ts_filename, ar_filename, fusion_filename]]

    gene_class_to_gene_region_bed = dict(zip(gene_classes, gene_region_beds))

    gene_class_to_results = {}
    for gene_class, gene_region_bed in gene_class_to_gene_region_bed:
        gene_class_to_results[gene_class] = \
            predict_effects_for_class(svs_bed, gene_class, gene_region_bed)


def predict_effects_for_class(svs_bed, gene_class, gene_region_bed):
    """
    Predict the effects of svs for all genes specified in gene_region_bed, given
    the specified gene_class.

    :param svs_bed: BedTool object providing structural variant coords
    :param gene_class: GeneClass enumeration value
    :param gene_region_bed: BedTool object providing coords of functional regions for genes
    :return: Dictionary linking gene to predicted SvEffect enumeration value
    """

    # XXX CONTINUE HERE: PROCESS THE GENE_REGION_BED TO EXTRACT UNIQUE GENES
    for gene in gene_region_bed:
    return {}

































