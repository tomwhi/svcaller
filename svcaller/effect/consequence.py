from enum import Enum
import pandas as pd
import sys
from svcaller.calling.events import SvType


class SvEffect(Enum):
    NO_OVERLAP = "NO_OVERLAP"
    OVERLAP_WITH_EFFECT = "OVERLAP_WITH_EFFECT"
    OVERLAP_UNKNOWN_EFFECT = "OVERLAP_UNKNOWN_EFFECT"


class GeneClass(Enum):
    TUMOUR_SUPRESSOR = "TUMOUR_SUPRESSOR"
    AR = "AR"
    FUSION_CANDIDATE = "FUSION_CANDIDATE"


class InvalidSvEffectException(Exception):
    pass


def parse_bed_to_dict(bed_file):
    df = pd.read_table(bed_file, header=None, sep="\t")
    return {name: table for name, table in df.groupby(df[3])}


def region1_overlaps_region2(region1, region2):
    return (region1[0] == region2[0]) and (region1[2] > region2[1]) and (region1[1] < region2[2])


def sv_in_regions(sv, regions):
    return (sv.iloc[0] == regions.iloc[0,0]) and \
           (sv.iloc[2] > min(regions.iloc[:,1])) and \
           (sv.iloc[1] < max(regions.iloc[:,2]))


def predict_del_effect(sv, _, functional_regions):
    # Ignore gene_class; just determine effect based on overlap with any of
    # the specified functional_regions.

    assert sv[3] == SvType.DEL.value

    curr_prediction = SvEffect.NO_OVERLAP
    if sv_in_regions(sv, functional_regions):
        curr_prediction = SvEffect.OVERLAP_UNKNOWN_EFFECT

    region_overlap_values = [region1_overlaps_region2(sv, region)
                             for (_, region) in functional_regions.iterrows()]
    if any(region_overlap_values):
        curr_prediction = SvEffect.OVERLAP_WITH_EFFECT

    return curr_prediction


def predict_inv_effect(sv, gene_class, functional_regions):
    return None


def predict_dup_effect(sv, gene_class, functional_regions):
    return None


def predict_tra_effect(sv, gene_class, functional_regions):
    return None


def collapse_sv_predictions(sv_effects):
    if not all([effect in SvEffect for effect in sv_effects]):
        raise InvalidSvEffectException("Input list {} includes invalid SvEffect value:".format(sv_effects))

    curr_max_effect = SvEffect.NO_OVERLAP
    for effect in sv_effects:
        if effect == SvEffect.OVERLAP_UNKNOWN_EFFECT:
            if curr_max_effect == SvEffect.NO_OVERLAP:
                curr_max_effect = SvEffect.OVERLAP_UNKNOWN_EFFECT
        elif effect == SvEffect.OVERLAP_WITH_EFFECT:
            curr_max_effect = SvEffect.OVERLAP_WITH_EFFECT

    return curr_max_effect


def predict_svs_gene_effect(svs, gene_class, gene_regions):
    # Data structure for emulating switch statement:
    svtype_to_function = {
        SvType.DEL: predict_del_effect,
        SvType.INV: predict_inv_effect,
        SvType.DUP: predict_dup_effect,
        SvType.TRA: predict_tra_effect,
    }

    # Currently calculate each sv effect individually and collapse down
    # to a single prediction:
    all_svs_effects = []
    for sv_type in [type_ for type_ in iter(SvType)]:
        print("TRACE: sv_type: {}".format(sv_type), file = sys.stderr)

        # Retrieve the SVs of the specified type:
        svs_of_type = svs.get(sv_type.value, pd.DataFrame({}))

        # Retrieve the relevant predictor function:
        predictor_function = svtype_to_function[sv_type]

        # Apply that function to each SV:
        for _, sv_row in svs_of_type.iterrows():
            print("TRACE: sv_row:", file=sys.stderr)
            print(sv_row, file=sys.stderr)
            all_svs_effects.append(predictor_function(sv_row, gene_class, gene_regions))

    return collapse_sv_predictions(all_svs_effects)


def filter_svtype_to_table(svtype_to_table, gene_regions):
    def sv_in_regions_tmp(sv):
        return sv_in_regions(sv, gene_regions)

    return {svtype: svs_table[svs_table.apply(sv_in_regions_tmp, axis=1)]
            for (svtype, svs_table) in svtype_to_table.items()}


def predict_svs_effects_for_class(svtype_to_table, gene_class, gene_to_table):
    gene_to_effect = {}
    for gene in gene_to_table:
        gene_regions = gene_to_table[gene]

        svtype_to_table_overlapping = filter_svtype_to_table(svtype_to_table, gene_regions)
        gene_svs_effect = predict_svs_gene_effect(svtype_to_table_overlapping, gene_class, gene_regions)
        gene_to_effect[gene] = gene_svs_effect

    return gene_to_effect


def predict_effects(svs_file, ts_file, ar_file, fusion_file):
    """
    Predict the consequence of the specified structural variants on the specified
    tumour suppressors, androgen receptor, and gene fusion candidate.

    :param svs_file: An open bed file specifying the structural variant coordinates.
    :param ts_file: An open bed file specifying the tumour supressor gene region coords.
    :param ar_file: An open bed file specifying Androgren Receptor gene region coords.
    :param fusion_file: An open bed file specifying two broad gene fusion region.

    :return: A dictionary with gene class as key and results dictionary as value
    """

    svtype_to_table = parse_bed_to_dict(svs_file)

    gene_classes = [GeneClass.TUMOUR_SUPRESSOR, GeneClass.AR, GeneClass.FUSION_CANDIDATE]
    gene_to_bed_tables = [parse_bed_to_dict(curr_file) for curr_file in
                          [ts_file, ar_file, fusion_file]]
    gene_class_to_gene_region_bed = dict(zip(gene_classes, gene_to_bed_tables))

    gene_class_to_results = {}
    for gene_class, gene_region_bed in gene_class_to_gene_region_bed.items():
        gene_class_to_results[gene_class] = \
            predict_svs_effects_for_class(svtype_to_table, gene_class, gene_region_bed)

    return gene_class_to_results
