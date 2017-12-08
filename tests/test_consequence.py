import unittest

from svcaller.effect.consequence import *
from unittest.mock import patch, mock_open
from collections import OrderedDict


class TestConsequenceFunctions(unittest.TestCase):
    def setUp(self):
        self.test_svs_content = '''1\t1000\t1100\tDEL\t0\t+
1\t3000\t3100\tDEL\t0\t+
2\t1000\t1100\tTRA\t0\t+
'''

        self.eg_svs_del_df = pd.DataFrame(OrderedDict({
            "chrom": [1, 1, 1, 1],
            "start": [1000, 3000, 1500, 1000],
            "end": [1100, 3100, 1600, 3000],
            "type": ["DEL", "DEL", "DEL", "DEL"],
            "score": [0, 0, 0, 0],
            "strand": ["+", "+", "+", "+"]
        }))

        self.eg_svs_inv_df = pd.DataFrame(OrderedDict({
            "chrom": [1, 1, 1, 1],
            "start": [1000, 1500, 1000, 1900],
            "end": [1100, 1600, 3000, 2200],
            "type": ["INV", "INV", "INV", "INV"],
            "score": [0, 0, 0, 0],
            "strand": ["+", "+", "+", "+"]
        }))

        self.eg_svs_dup_df = pd.DataFrame(OrderedDict({
            "chrom": [1, 1, 1, 1, 1],
            "start": [1000, 1500, 1900, 400, 400],
            "end": [1100, 1600, 2200, 900, 1300],
            "type": ["DUP", "DUP", "DUP", "DUP", "DUP"],
            "score": [0, 0, 0, 0, 0],
            "strand": ["+", "+", "+", "+", "+"]
        }))

        self.eg_svs_tra_df = pd.DataFrame(OrderedDict({
            "chrom": [1, 1],
            "start": [500, 1500],
            "end": [600, 1600],
            "type": ["TRA", "TRA"],
            "score": [0, 0],
            "strand": ["+", "+"],
        }))

        self.eg_regions_df_ts1 = pd.DataFrame(OrderedDict({
            "chrom": [1, 1],
            "start": [800, 2000],
            "end": [1200, 2100],
            "gene": ["TS1", "TS1"],
            "score": [0, 0],
            "strand": ["+", "+"],
        }))

        self.eg_regions_df_ts2 = pd.DataFrame(OrderedDict({
            "chrom": [1, 1],
            "start": [800, 2000],
            "end": [1200, 2100],
            "gene": ["TS2", "TS2"],
            "score": [0, 0],
            "strand": ["-", "-"],
        }))

        self.eg_regions_df_ar = pd.DataFrame(OrderedDict({
            "chrom": [1, 1, 1, 1],
            "start": [500, 800, 1700, 2000],
            "end": [500, 1200, 1800, 2100],
            "gene": ["AR", "AR", "AR", "AR"],
            "score": [0, 0, 0, 0],
            "strand": ["+", "+", "+", "+"],
        }))


        self.eg_svtype_to_table = {
            "DEL": pd.DataFrame(OrderedDict({
                "chrom": [1],
                "start": [1000],
                "end": [1100],
                "type": ["DEL"],
                "score": [0],
                "strand": ["+"]
            })),
            "DUP": pd.DataFrame(OrderedDict({
                "chrom": [1],
                "start": [1000],
                "end": [1100],
                "type": ["DUP"],
                "score": [0],
                "strand": ["+"]
            })),
            "INV": pd.DataFrame(OrderedDict({
                "chrom": [1],
                "start": [1000],
                "end": [1100],
                "type": ["INV"],
                "score": [0],
                "strand": ["+"]
            })),
            "TRA": pd.DataFrame(OrderedDict({
                "chrom": [1],
                "start": [1000],
                "end": [1100],
                "type": ["TRA"],
                "score": [0],
                "strand": ["+"]
            })),
        }

        self.eg_bed_content_1 = '''1\t800\t1200\tTS1
1\t2000\t2100\tTS1
2\t2000\t2100\tTS2
'''

        self.eg_bed_content_fusion = '''1\t10000\t20000\tFUSION1
1\t30000\t40000\tFUSION2
'''

    def test_parse_bed_to_dict(self):
        open_name = '%s.open' % __name__
        with patch(open_name, mock_open(read_data=self.test_svs_content), create=True):
            with open("dummy_filename.bed") as test_svs:
                type2df = parse_bed_to_dict(test_svs)
                self.assertIn("DEL", type2df)
                self.assertIn("TRA", type2df)
                self.assertEquals(type2df["DEL"].iloc[0,1], 1000)

    def test_region1_overlaps_regions2_1(self):
        self.assertTrue(region1_overlaps_region2(('1', 1000, 2000), ('1', 500, 1500)))
        self.assertTrue(region1_overlaps_region2(('1', 1000, 2000), ('1', 500, 2500)))
        self.assertFalse(region1_overlaps_region2(('1', 1000, 2000), ('1', 2500, 3000)))
        self.assertFalse(region1_overlaps_region2(('1', 1000, 2000), ('1', 500, 800)))

    def test_sv_in_regions_true(self):
        self.assertTrue(sv_in_regions(
            self.eg_svs_del_df.iloc[0, :],
            self.eg_regions_df_ts1
        ))

    def test_sv_in_regions_false(self):
        self.assertFalse(sv_in_regions(
            self.eg_svs_del_df.iloc[1, :],
            self.eg_regions_df_ts1
        ))

    def test_collapse_sv_predictions_1(self):
        self.assertEquals(collapse_sv_predictions([SvEffect.NO_EFFECT]), SvEffect.NO_EFFECT)

    def test_collapse_sv_predictions_2(self):
        self.assertEquals(
            collapse_sv_predictions([SvEffect.NO_EFFECT, SvEffect.OVERLAP_UNKNOWN_EFFECT, SvEffect.NO_EFFECT]),
            SvEffect.OVERLAP_UNKNOWN_EFFECT
        )

    def test_collapse_sv_predictions_3(self):
        self.assertEquals(
            collapse_sv_predictions([SvEffect.NO_EFFECT, SvEffect.OVERLAP_UNKNOWN_EFFECT, SvEffect.OVERLAP_WITH_EFFECT]),
            SvEffect.OVERLAP_WITH_EFFECT
        )

    def test_predict_del_effect_with_effect(self):
        self.assertEquals(
            predict_del_effect(self.eg_svs_del_df.iloc[0, :], GeneClass.TUMOUR_SUPRESSOR, self.eg_regions_df_ts1),
            SvEffect.OVERLAP_WITH_EFFECT
        )

    def test_predict_del_effect_no_effect(self):
        self.assertEquals(
            predict_del_effect(self.eg_svs_del_df.iloc[1, :], GeneClass.TUMOUR_SUPRESSOR, self.eg_regions_df_ts1),
            SvEffect.NO_EFFECT
        )

    def test_predict_del_effect_unknown_effect(self):
        self.assertEquals(
            predict_del_effect(self.eg_svs_del_df.iloc[2, :], GeneClass.TUMOUR_SUPRESSOR, self.eg_regions_df_ts1),
            SvEffect.OVERLAP_UNKNOWN_EFFECT
        )

    def test_predict_del_effect_overlap_first_exon_ts(self):
        self.assertEquals(
            predict_del_effect(self.eg_svs_del_df.iloc[0, :], GeneClass.TUMOUR_SUPRESSOR, self.eg_regions_df_ts1),
            SvEffect.OVERLAP_WITH_EFFECT
        )

    def test_predict_del_effect_overlap_non_lbd_ar(self):
        self.assertEquals(
            predict_del_effect(self.eg_svs_del_df.iloc[0, :], GeneClass.AR, self.eg_regions_df_ar),
            SvEffect.NO_EFFECT
        )

    def test_predict_del_effect_overlap_non_lbd_to_first_lbd_ar(self):
        self.assertEquals(
            predict_del_effect(self.eg_svs_del_df.iloc[3, :], GeneClass.AR, self.eg_regions_df_ar),
            SvEffect.OVERLAP_UNKNOWN_EFFECT
        )

    def test_predict_inv_effect_first_region(self):
        self.assertEquals(
            predict_inv_effect(self.eg_svs_inv_df.iloc[0, :], GeneClass.TUMOUR_SUPRESSOR, self.eg_regions_df_ar),
            SvEffect.OVERLAP_WITH_EFFECT
        )

    def test_predict_inv_effect_in_the_middle(self):
        self.assertEquals(
            predict_inv_effect(self.eg_svs_inv_df.iloc[1, :], GeneClass.TUMOUR_SUPRESSOR, self.eg_regions_df_ar),
            SvEffect.OVERLAP_UNKNOWN_EFFECT
        )

    def test_predict_inv_effect_first_and_last_overlap(self):
        self.assertEquals(
            predict_inv_effect(self.eg_svs_inv_df.iloc[2, :], GeneClass.TUMOUR_SUPRESSOR, self.eg_regions_df_ar),
            SvEffect.OVERLAP_WITH_EFFECT
        )

    def test_predict_inv_effect_first_region_ar(self):
        self.assertEquals(
            predict_inv_effect(self.eg_svs_inv_df.iloc[0, :], GeneClass.AR, self.eg_regions_df_ar),
            SvEffect.NO_EFFECT
        )

    def test_predict_inv_effect_in_the_middle_ar(self):
        test_sv = self.eg_svs_inv_df.iloc[1, :]
        self.assertEquals(
            predict_inv_effect(test_sv, GeneClass.AR, self.eg_regions_df_ar),
            SvEffect.NO_EFFECT
        )

    def test_predict_inv_effect_first_and_last_overlap_ar(self):
        self.assertEquals(
            predict_inv_effect(self.eg_svs_inv_df.iloc[2, :], GeneClass.AR, self.eg_regions_df_ar),
            SvEffect.OVERLAP_UNKNOWN_EFFECT
        )

    def test_predict_inv_effect_last_overlap_ar(self):
        self.assertEquals(
            predict_inv_effect(self.eg_svs_inv_df.iloc[3, :], GeneClass.AR, self.eg_regions_df_ar),
            SvEffect.OVERLAP_WITH_EFFECT
        )

    def test_predict_dup_effect_first_region(self):
        self.assertEquals(
            predict_dup_effect(self.eg_svs_dup_df.iloc[0, :], GeneClass.TUMOUR_SUPRESSOR, self.eg_regions_df_ts1),
            SvEffect.OVERLAP_WITH_EFFECT
        )

    def test_predict_dup_effect_in_the_middle(self):
        self.assertEquals(
            predict_dup_effect(self.eg_svs_dup_df.iloc[1, :], GeneClass.TUMOUR_SUPRESSOR, self.eg_regions_df_ts1),
            SvEffect.OVERLAP_UNKNOWN_EFFECT
        )

    def test_predict_dup_effect_first_and_last_overlap(self):
        self.assertEquals(
            predict_dup_effect(self.eg_svs_dup_df.iloc[2, :], GeneClass.TUMOUR_SUPRESSOR, self.eg_regions_df_ts1),
            SvEffect.OVERLAP_UNKNOWN_EFFECT
        )

    def test_predict_dup_effect_3(self):
        self.assertEquals(
            predict_dup_effect(self.eg_svs_dup_df.iloc[3, :], GeneClass.TUMOUR_SUPRESSOR, self.eg_regions_df_ts1),
            SvEffect.OVERLAP_UNKNOWN_EFFECT
        )

    def test_predict_dup_effect_4(self):
        self.assertEquals(
            predict_dup_effect(self.eg_svs_dup_df.iloc[4, :], GeneClass.TUMOUR_SUPRESSOR, self.eg_regions_df_ts1),
            SvEffect.OVERLAP_UNKNOWN_EFFECT
        )

    def test_predict_dup_ar_0(self):
        test_sv = self.eg_svs_dup_df.iloc[0, :]
        self.assertEquals(
            predict_dup_effect(test_sv, GeneClass.AR, self.eg_regions_df_ar),
            SvEffect.NO_EFFECT
        )

    def test_predict_dup_ar_1(self):
        test_sv = self.eg_svs_dup_df.iloc[1, :]
        self.assertEquals(
            predict_dup_effect(test_sv, GeneClass.AR, self.eg_regions_df_ar),
            SvEffect.NO_EFFECT
        )

    def test_predict_dup_ar_2(self):
        test_sv = self.eg_svs_dup_df.iloc[2, :]
        self.assertEquals(
            predict_dup_effect(test_sv, GeneClass.AR, self.eg_regions_df_ar),
            SvEffect.OVERLAP_UNKNOWN_EFFECT
        )

    def test_predict_dup_ar_3(self):
        test_sv = self.eg_svs_dup_df.iloc[3, :]
        self.assertEquals(
            predict_dup_effect(test_sv, GeneClass.AR, self.eg_regions_df_ar),
            SvEffect.NO_EFFECT
        )

    def test_predict_dup_ar_4(self):
        test_sv = self.eg_svs_dup_df.iloc[4, :]
        self.assertEquals(
            predict_dup_effect(test_sv, GeneClass.AR, self.eg_regions_df_ar),
            SvEffect.OVERLAP_WITH_EFFECT
        )

    def test_predict_tra_ts1_0(self):
        test_sv = self.eg_svs_tra_df.iloc[0, :]
        self.assertEquals(
            predict_tra_effect(test_sv, GeneClass.TUMOUR_SUPRESSOR, self.eg_regions_df_ts1),
            SvEffect.NO_EFFECT
        )

    def test_predict_tra_ts1_1(self):
        test_sv = self.eg_svs_tra_df.iloc[1, :]
        self.assertEquals(
            predict_tra_effect(test_sv, GeneClass.TUMOUR_SUPRESSOR, self.eg_regions_df_ts1),
            SvEffect.OVERLAP_WITH_EFFECT
        )

    def test_predict_tra_ts2_0(self):
        test_sv = self.eg_svs_tra_df.iloc[0, :]
        self.assertEquals(
            predict_tra_effect(test_sv, GeneClass.TUMOUR_SUPRESSOR, self.eg_regions_df_ts2),
            SvEffect.NO_EFFECT
        )

    def test_predict_tra_ts2_1(self):
        test_sv = self.eg_svs_tra_df.iloc[1, :]
        self.assertEquals(
            predict_tra_effect(test_sv, GeneClass.TUMOUR_SUPRESSOR, self.eg_regions_df_ts2),
            SvEffect.OVERLAP_UNKNOWN_EFFECT
        )

    def test_predict_tra_ar_0(self):
        test_sv = self.eg_svs_tra_df.iloc[0, :]
        self.assertEquals(
            predict_tra_effect(test_sv, GeneClass.AR, self.eg_regions_df_ar),
            SvEffect.NO_EFFECT
        )

    def test_predict_tra_ar_1(self):
        test_sv = self.eg_svs_tra_df.iloc[1, :]
        self.assertEquals(
            predict_tra_effect(test_sv, GeneClass.AR, self.eg_regions_df_ar),
            SvEffect.OVERLAP_WITH_EFFECT
        )

    def test_filter_svtype_to_table(self):
        svtype_to_table = {
            "DEL": pd.DataFrame(OrderedDict({
            "chrom": [1, 1, 1],
            "start": [1000, 3000, 1500],
            "end": [1100, 3100, 1600],
            "type": ["DEL", "DEL", "DEL"],
            "score": [0, 0, 0],
            "strand": ["+", "+", "+"]
        }))}

        filtered = filter_svtype_to_table(svtype_to_table, self.eg_regions_df_ts1)
        self.assertIn("DEL", filtered)
        self.assertEquals(2, len(list(filtered.values())[0]))

    @patch('svcaller.effect.consequence.predict_del_effect')
    @patch('svcaller.effect.consequence.predict_dup_effect')
    @patch('svcaller.effect.consequence.predict_inv_effect')
    @patch('svcaller.effect.consequence.predict_tra_effect')
    def test_predict_svs_gene_effect_1(self, mock_predict_tra_effect, mock_predict_inv_effect,
                                     mock_predict_dup_effect, mock_predict_del_effect):
        mock_predict_del_effect.return_value = SvEffect.NO_EFFECT
        mock_predict_dup_effect.return_value = SvEffect.NO_EFFECT
        mock_predict_inv_effect.return_value = SvEffect.NO_EFFECT
        mock_predict_tra_effect.return_value = SvEffect.NO_EFFECT

        self.assertEquals(
            predict_svs_gene_effect(self.eg_svtype_to_table, GeneClass.TUMOUR_SUPRESSOR, "Dummy"),
            SvEffect.NO_EFFECT
        )

    @patch('svcaller.effect.consequence.predict_del_effect')
    @patch('svcaller.effect.consequence.predict_dup_effect')
    @patch('svcaller.effect.consequence.predict_inv_effect')
    @patch('svcaller.effect.consequence.predict_tra_effect')
    def test_predict_svs_gene_effect_2(self, mock_predict_tra_effect, mock_predict_inv_effect,
                                     mock_predict_dup_effect, mock_predict_del_effect):
        mock_predict_del_effect.return_value = SvEffect.OVERLAP_UNKNOWN_EFFECT
        mock_predict_dup_effect.return_value = SvEffect.NO_EFFECT
        mock_predict_inv_effect.return_value = SvEffect.NO_EFFECT
        mock_predict_tra_effect.return_value = SvEffect.NO_EFFECT

        self.assertEquals(
            predict_svs_gene_effect(self.eg_svtype_to_table, GeneClass.TUMOUR_SUPRESSOR, "Dummy"),
            SvEffect.OVERLAP_UNKNOWN_EFFECT
        )

    @patch('svcaller.effect.consequence.predict_del_effect')
    @patch('svcaller.effect.consequence.predict_dup_effect')
    @patch('svcaller.effect.consequence.predict_inv_effect')
    @patch('svcaller.effect.consequence.predict_tra_effect')
    def test_predict_svs_gene_effect_3(self, mock_predict_tra_effect, mock_predict_inv_effect,
                                     mock_predict_dup_effect, mock_predict_del_effect):
        mock_predict_del_effect.return_value = SvEffect.NO_EFFECT
        mock_predict_dup_effect.return_value = SvEffect.OVERLAP_WITH_EFFECT
        mock_predict_inv_effect.return_value = SvEffect.OVERLAP_UNKNOWN_EFFECT
        mock_predict_tra_effect.return_value = SvEffect.NO_EFFECT

        self.assertEquals(
            predict_svs_gene_effect(self.eg_svtype_to_table, GeneClass.TUMOUR_SUPRESSOR, "Dummy"),
            SvEffect.OVERLAP_WITH_EFFECT
        )

    @patch('svcaller.effect.consequence.predict_svs_gene_effect')
    @patch('svcaller.effect.consequence.filter_svtype_to_table')
    def test_predict_svs_effects_for_class(self, mock_filter_svtype_to_table, mock_predict_svs_gene_effect):
        mock_filter_svtype_to_table.return_value = {}
        mock_predict_svs_gene_effect.return_value = SvEffect.OVERLAP_UNKNOWN_EFFECT

        gene_to_table = {
            "TS1": pd.DataFrame(OrderedDict({
                "chrom": [1, 1],
                "start": [800, 2000],
                "end": [1200, 2100],
                "gene": ["TS1", "TS1"],
            }))
        }

        gene_to_effect = predict_svs_effects_for_class({}, None, gene_to_table)

        self.assertIn("TS1", gene_to_effect)

    @patch('svcaller.effect.consequence.predict_svs_effects_for_class')
    def test_predict_effects(self, mock_predict_svs_effects_for_class):
        mock_predict_svs_effects_for_class.return_value = {}
        open_name = '%s.open' % __name__
        with patch(open_name, mock_open(read_data=self.test_svs_content), create=True) as open1, \
            patch(open_name, mock_open(read_data=self.eg_bed_content_1), create=True) as open2, \
            patch(open_name, mock_open(read_data=self.eg_bed_content_1), create=True) as open3, \
            patch(open_name, mock_open(read_data=self.eg_bed_content_fusion), create=True) as open4:
            with open1("dummy_svs.bed") as test_svs_file, \
                open2("dummy_ts.bed") as test_ts_file, \
                open3("dummy_ar.bed") as test_ar_file, \
                open4("dummy_fusion.bed") as test_fusion_file:
                gene_class_to_results = \
                    predict_effects(test_svs_file, test_ts_file, test_ar_file, test_fusion_file)
                self.assertIn(GeneClass.TUMOUR_SUPRESSOR.value, gene_class_to_results)
                self.assertIn(GeneClass.AR.value, gene_class_to_results)
                self.assertIn(GeneClass.FUSION_CANDIDATE.value, gene_class_to_results)