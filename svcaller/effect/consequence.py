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


def region1_overlaps_regions2(region1, regions2):
    return any([region1_overlaps_region2(region1, region2)
                for (_, region2) in regions2.iterrows()])


def predict_effect_overlap(sv, functional_regions):
    """
    Predict the effect of the given structural variant on the
    functional regions, based solely on overlap of the SV
    with any of those regions.

    :param sv: Series object representing a structural variant
    :param functional_regions: Data frame representing the functional regions
    :return: SvEffect enumeration value
    """
    curr_prediction = SvEffect.NO_OVERLAP
    if sv_in_regions(sv, functional_regions):
        curr_prediction = SvEffect.OVERLAP_UNKNOWN_EFFECT

    if region1_overlaps_regions2(sv, functional_regions):
        curr_prediction = SvEffect.OVERLAP_WITH_EFFECT

    return curr_prediction


# FIXME: Representing AR TSS and last non-LBD exon as the first two
# gene_regions rows. This seems potentially inelegant and bug-prone, and
# affects the functions below. Doing this in order to be able to re-use
# much of the code for given SV type effects across the two classes (tumour
# suppressor and AR). However, there is probably a cleaner solution.
def get_ar_regions(gene_regions):
    assert gene_regions.iloc[0,1] == gene_regions.iloc[0,2]
    return (gene_regions.iloc[0,:], gene_regions.iloc[1,:],
            gene_regions.iloc[2:,:])


def predict_del_effect(sv, gene_class, gene_regions):
    assert sv[3] == SvType.DEL.value

    functional_regions = gene_regions
    if gene_class == GeneClass.AR:
        (_, last_non_lbd_region, functional_regions) = get_ar_regions(gene_regions)

    prediction = predict_effect_overlap(sv, functional_regions)

    # The deletion must occur after the end of the last non-LBD region
    # for the deletion to have an effect in the case of AR:
    if gene_class == GeneClass.AR and prediction == SvEffect.OVERLAP_WITH_EFFECT:
        assert sv[0] == gene_regions.iloc[0,0]
        if sv[1] < last_non_lbd_region[2]:
            prediction = SvEffect.OVERLAP_UNKNOWN_EFFECT

    return prediction


def predict_inv_effect(sv, gene_class, gene_regions):
    assert sv[3] == SvType.INV.value

    functional_regions = gene_regions
    if gene_class == GeneClass.AR:
        (_, last_non_lbd_region, functional_regions) = get_ar_regions(gene_regions)

    prediction = predict_effect_overlap(sv, functional_regions)

    if prediction == SvEffect.OVERLAP_WITH_EFFECT:
        assert sv[0] == gene_regions.iloc[0,0]
        if gene_class == GeneClass.TUMOUR_SUPRESSOR:
            # If the inversion spans the *entire* gene region then it is counted
            # as an unknown significance overlap:
            if sv[1] < functional_regions.iloc[0,1] and sv[2] > functional_regions.iloc[-1,2]:
                prediction = SvEffect.OVERLAP_UNKNOWN_EFFECT
        if gene_class == GeneClass.AR:
            # If the inversion includes some of the last non-lbd region, change call
            # to unknown significance:
            if sv[1] < last_non_lbd_region[2]:
                prediction = SvEffect.OVERLAP_UNKNOWN_EFFECT

    return prediction


def predict_dup_effect(sv, gene_class, gene_regions):
    assert sv[3] == SvType.DUP.value

    functional_regions = gene_regions
    if gene_class == GeneClass.AR:
        (tss, last_non_lbd_region, functional_regions) = get_ar_regions(gene_regions)

    prediction = predict_effect_overlap(sv, functional_regions)

    # FIXME: It feels like there must be a more elegant, bug-proof, readable
    # way to do this (i.e. dealing with the various DUP overlap scenarios):
    if prediction == SvEffect.OVERLAP_WITH_EFFECT:
        assert sv[0] == gene_regions.iloc[0,0]
        if gene_class == GeneClass.TUMOUR_SUPRESSOR:
            # If the duplication overhangs either end of the overall functional
            # region, then predict unknown effect:
            if (sv[1] < functional_regions.iloc[0,1]) or \
               (sv[2] > functional_regions.iloc[-1,2]):
                prediction = SvEffect.OVERLAP_UNKNOWN_EFFECT
        if gene_class == GeneClass.AR:
            # If the duplication overlaps the last non-lbd region or the end
            # of the AR region, then change call to unknown significance:
            if (sv[1] < last_non_lbd_region[2]) or \
               (sv[2] > functional_regions.iloc[-1, 2]):
                prediction = SvEffect.OVERLAP_UNKNOWN_EFFECT

    # Deal with separate scenario for AR, in which a duplication that includes the
    # TSS and last non-LBD region - but excludes part of the LB functional region -
    # could result in an mRNA lacking the LBD, depending on the size of the
    # duplicated region preceding the TSS:
    if gene_class == GeneClass.AR:
        assert sv[0] == tss[0]
        if (sv[1] < tss[1]) and (sv[2] > last_non_lbd_region[2]) and \
           (sv[2] < functional_regions.iloc[-1, 2]):
            prediction = SvEffect.OVERLAP_WITH_EFFECT

    return prediction


def predict_tra_effect(sv, gene_class, gene_regions):
    assert sv[3] == SvType.TRA.value

    functional_regions = gene_regions
    if gene_class == GeneClass.AR:
        (_, last_non_lbd_region, functional_regions) = get_ar_regions(gene_regions)

    assert sv[0] == functional_regions.iloc[0,0]

    prediction = SvEffect.NO_OVERLAP

    tra_point_position = sv[2]
    if sv[5] == "-":
        tra_point_position = sv[1]

    # Fall back to calling unknown significance overlap if the TRA occurs
    # somewhere within the functional region:
    if (tra_point_position > functional_regions.iloc[0,1]) and \
       (tra_point_position < functional_regions.iloc[-1,2]):
        prediction = SvEffect.OVERLAP_UNKNOWN_EFFECT

    # For determining if the SV has a known effect, only consider TRA events
    # pointing in the same direction as the gene, irrespective of gene class:
    if sv[5] == gene_regions.iloc[0, 5]:
        if gene_class == GeneClass.TUMOUR_SUPRESSOR:
            # A TRA anywhere in the functional region pointing the same
            # direction as the gene transcription will be called as having
            # an effect in the case of tumour supressors:
            if (tra_point_position > functional_regions.iloc[0,1]) and \
               (tra_point_position < functional_regions.iloc[-1,2]):
                prediction = SvEffect.OVERLAP_WITH_EFFECT
        elif gene_class == GeneClass.AR:
            # A TRA pointing the same direction as the gene transcription
            # and occuring after the last non-LBD region but before the
            # end of the LBD functional region will be called as having
            # an effect in the case of AR:
            if (tra_point_position > last_non_lbd_region[2]) and \
               (tra_point_position < functional_regions.iloc[-1,2]):
                prediction = SvEffect.OVERLAP_WITH_EFFECT

    return prediction


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
        # Retrieve the SVs of the specified type:
        svs_of_type = svs.get(sv_type.value, pd.DataFrame({}))

        # Retrieve the relevant predictor function:
        predictor_function = svtype_to_function[sv_type]

        # Apply that function to each SV:
        for _, sv_row in svs_of_type.iterrows():
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
