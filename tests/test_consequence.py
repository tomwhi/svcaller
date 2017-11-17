import unittest

from svcaller.effect.consequence import *
from unittest.mock import patch, mock_open


class TestConsequenceFunctions(unittest.TestCase):
    def setUp(self):
        self.test_svs_content = '''1\t1000\t1100\tDEL\t0\t+
1\t3000\t3100\tDEL\t0\t+
2\t1000\t1100\tTRA\t0\t+
'''

        self.test_svs_dict = pd.DataFrame({
    "chrom": [1, 1, 2],
    "start": [1000, 3000, 1000],
    "end": [1100, 3100, 1100],
    "type": ["DEL", "DEL", "TRA"],
    "score": [0, 0, 0],
    "strand": ["+", "+", "+"]
})

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

    def test_predict_effects(self):
        self.assertEquals
        