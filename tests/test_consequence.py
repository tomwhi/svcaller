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

        self.test_svs_df = pd.DataFrame(OrderedDict({
            "chrom": [1, 1, 2],
            "start": [1000, 3000, 1000],
            "end": [1100, 3100, 1100],
            "type": ["DEL", "DEL", "TRA"],
            "score": [0, 0, 0],
            "strand": ["+", "+", "+"]
        }))

        self.test_regions_df_ts1 = pd.DataFrame(OrderedDict({
            "chrom": [1, 1],
            "start": [800, 2000],
            "end": [1200, 2100],
            "gene": ["TS1", "TS1"],
        }))

        self.test_bed_content_1 = '''1\t800\t1200\tTS1
1\t2000\t2100\tTS1
2\t2000\t2100\tTS2
'''

        self.test_bed_content_fusion = '''1\t10000\t20000\tFUSION1
1\t30000\t40000\tFUSION2
'''

    def test_extract_groups(self):
        open_name = '%s.open' % __name__
        with patch(open_name, mock_open(read_data=self.test_svs_content), create=True):
            with open("dummy_filename.bed") as test_svs:
                type2df = extract_groups(test_svs)
                self.assertIn("DEL", type2df)
                self.assertIn("TRA", type2df)
                self.assertEquals(type2df["DEL"].iloc[0,1], 1000)

    def test_region1_overlaps_regions2_1(self):
        self.assertTrue(region1_overlaps_region2(('1', 1000, 2000), ('1', 500, 1500)))
        self.assertTrue(region1_overlaps_region2(('1', 1000, 2000), ('1', 500, 2500)))
        self.assertFalse(region1_overlaps_region2(('1', 1000, 2000), ('1', 2500, 3000)))
        self.assertFalse(region1_overlaps_region2(('1', 1000, 2000), ('1', 500, 800)))

    @patch('svcaller.effect.consequence.predict_svs_effects_for_class')
    def test_predict_effects(self, mock_predict_svs_effects_for_class):
        mock_predict_svs_effects_for_class.return_value = {}
        open_name = '%s.open' % __name__
        with patch(open_name, mock_open(read_data=self.test_svs_content), create=True) as open1, \
            patch(open_name, mock_open(read_data=self.test_bed_content_1), create=True) as open2, \
            patch(open_name, mock_open(read_data=self.test_bed_content_1), create=True) as open3, \
            patch(open_name, mock_open(read_data=self.test_bed_content_fusion), create=True) as open4:
            with open1("dummy_svs.bed") as test_svs_file, \
                open2("dummy_ts.bed") as test_ts_file, \
                open3("dummy_ar.bed") as test_ar_file, \
                open4("dummy_fusion.bed") as test_fusion_file:
                gene_class_to_results = \
                    predict_effects(test_svs_file, test_ts_file, test_ar_file, test_fusion_file)
                self.assertIn(GeneClass.TUMOUR_SUPRESSOR, gene_class_to_results)
                self.assertIn(GeneClass.AR, gene_class_to_results)
                self.assertIn(GeneClass.FUSION_CANDIDATE, gene_class_to_results)

    def test_sv_in_regions_true(self):
        self.assertTrue(sv_in_regions(
            self.test_svs_df.iloc[0,:],
            self.test_regions_df_ts1
        ))

    def test_sv_in_regions_false(self):
        self.assertFalse(sv_in_regions(
            self.test_svs_df.iloc[1,:],
            self.test_regions_df_ts1
        ))

    def test_collapse_sv_predictions_1(self):
        self.assertEquals(collapse_sv_predictions([SvEffect.NO_OVERLAP]), SvEffect.NO_OVERLAP)

    def test_collapse_sv_predictions_2(self):
        self.assertEquals(
            collapse_sv_predictions([SvEffect.NO_OVERLAP, SvEffect.OVERLAP_UNKNOWN_EFFECT, SvEffect.NO_OVERLAP]),
            SvEffect.OVERLAP_UNKNOWN_EFFECT
        )

    def test_collapse_sv_predictions_3(self):
        self.assertEquals(
            collapse_sv_predictions([SvEffect.NO_OVERLAP, SvEffect.OVERLAP_UNKNOWN_EFFECT, SvEffect.OVERLAP_WITH_EFFECT]),
            SvEffect.OVERLAP_WITH_EFFECT
        )