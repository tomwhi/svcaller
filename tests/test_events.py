import unittest

from svcaller.calling.events import *
from unittest.mock import patch, mock_open


class TestConsequenceFunctions(unittest.TestCase):
    def setUp(self):
        self.test_sv_gtf_content_1 = '''
17	SV_event	exon	7184103	7184551	50	+	.	gene_id "17:7184103-7184551,17:7576474-7576873"; transcript_id "17:7184103-7184551,17:7576474-7576873";
18	SV_event	exon	7576474	7576873	50	-	.	gene_id "17:7184103-7184551,17:7576474-7576873"; transcript_id "17:7184103-7184551,17:7576474-7576873";
17	SV_event	CDS	7184534	7184551	1	+	.	gene_id "17:7184103-7184551,17:7576474-7576873"; transcript_id "17:7184103-7184551,17:7576474-7576873";
17	SV_event	CDS	7576476	7576516	1	-	.	gene_id "17:7184103-7184551,17:7576474-7576873"; transcript_id "17:7184103-7184551,17:7576474-7576873";

21	SV_event	exon	39874777	39875065	46	+	.	gene_id "21:39874777-39875065,21:42871970-42872366"; transcript_id "21:39874777-39875065,21:42871970-42872366";
20	SV_event	exon	42871970	42872366	46	-	.	gene_id "21:39874777-39875065,21:42871970-42872366"; transcript_id "21:39874777-39875065,21:42871970-42872366";
21	SV_event	CDS	39875015	39875064	1	+	.	gene_id "21:39874777-39875065,21:42871970-42872366"; transcript_id "21:39874777-39875065,21:42871970-42872366";
21	SV_event	CDS	42871971	42872014	1	-	.	gene_id "21:39874777-39875065,21:42871970-42872366"; transcript_id "21:39874777-39875065,21:42871970-42872366";

'''

        self.test_sv_gtf_content_2 = '''
17	SV_event	exon	7184103	7184551	50	+	.	gene_id "17:7184103-7184551,17:7576474-7576873"; transcript_id "17:7184103-7184551,17:7576474-7576873";
17	SV_event	exon	7576474	7576873	50	-	.	gene_id "17:7184103-7184551,17:7576474-7576873"; transcript_id "17:7184103-7184551,17:7576474-7576873";
17	SV_event	CDS	7184534	7184551	1	+	.	gene_id "17:7184103-7184551,17:7576474-7576873"; transcript_id "17:7184103-7184551,17:7576474-7576873";
17	SV_event	CDS	7576476	7576516	1	-	.	gene_id "17:7184103-7184551,17:7576474-7576873"; transcript_id "17:7184103-7184551,17:7576474-7576873";

21	SV_event	exon	39874777	39875065	46	+	.	gene_id "21:39874777-39875065,21:42871970-42872366"; transcript_id "21:39874777-39875065,21:42871970-42872366";
21	SV_event	exon	42871970	42872366	46	-	.	gene_id "21:39874777-39875065,21:42871970-42872366"; transcript_id "21:39874777-39875065,21:42871970-42872366";
21	SV_event	CDS	39875015	39875064	1	+	.	gene_id "21:39874777-39875065,21:42871970-42872366"; transcript_id "21:39874777-39875065,21:42871970-42872366";
21	SV_event	CDS	42871971	42872014	1	-	.	gene_id "21:39874777-39875065,21:42871970-42872366"; transcript_id "21:39874777-39875065,21:42871970-42872366";

'''

        self.test_sv_gtf_content_3 = '''
    17	SV_event	exon	7184103	7184551	50	+	.	gene_id "17:7184103-7184551,17:7576474-7576873"; transcript_id "17:7184103-7184551,17:7576474-7576873";
    17	SV_event	exon	7576474	7576873	50	+	.	gene_id "17:7184103-7184551,17:7576474-7576873"; transcript_id "17:7184103-7184551,17:7576474-7576873";
    17	SV_event	CDS	7184534	7184551	1	+	.	gene_id "17:7184103-7184551,17:7576474-7576873"; transcript_id "17:7184103-7184551,17:7576474-7576873";
    17	SV_event	CDS	7576476	7576516	1	+	.	gene_id "17:7184103-7184551,17:7576474-7576873"; transcript_id "17:7184103-7184551,17:7576474-7576873";

    21	SV_event	exon	39874777	39875065	46	+	.	gene_id "21:39874777-39875065,21:42871970-42872366"; transcript_id "21:39874777-39875065,21:42871970-42872366";
    21	SV_event	exon	42871970	42872366	46	+	.	gene_id "21:39874777-39875065,21:42871970-42872366"; transcript_id "21:39874777-39875065,21:42871970-42872366";
    21	SV_event	CDS	39875015	39875064	1	+	.	gene_id "21:39874777-39875065,21:42871970-42872366"; transcript_id "21:39874777-39875065,21:42871970-42872366";
    21	SV_event	CDS	42871971	42872014	1	+	.	gene_id "21:39874777-39875065,21:42871970-42872366"; transcript_id "21:39874777-39875065,21:42871970-42872366";

    '''

    def test_read_sv_gtf_1(self):
        open_name = '%s.open' % __name__
        with patch(open_name, mock_open(read_data=self.test_sv_gtf_content_2), create=True):
            with open("dummy_filename.gtf") as test_gtf_file:
                svs_bed_table = read_sv_gtf(test_gtf_file, SvType.DEL)
                self.assertEquals(len(svs_bed_table), 2)
                self.assertEquals(svs_bed_table.iloc[0, 1], 7184103)
                self.assertEquals(svs_bed_table.iloc[0, 2], 7576873)
                self.assertEquals(svs_bed_table.iloc[0, 5], None)


    def test_read_sv_gtf_2(self):
        open_name = '%s.open' % __name__
        with patch(open_name, mock_open(read_data=self.test_sv_gtf_content_1), create=True):
            with open("dummy_filename.gtf") as test_gtf_file:
                svs_bed_table = read_sv_gtf(test_gtf_file, SvType.TRA)
                self.assertEquals(len(svs_bed_table), 4)
                self.assertEquals(svs_bed_table.iloc[0, 1], 7184103)
                self.assertEquals(svs_bed_table.iloc[0, 2], 7184551)
                self.assertEquals(svs_bed_table.iloc[0, 5], "+")
                self.assertEquals(svs_bed_table.iloc[1, 1], 7576474)
                self.assertEquals(svs_bed_table.iloc[1, 2], 7576873)
                self.assertEquals(svs_bed_table.iloc[1, 5], "-")


    def test_read_sv_gtf_3(self):
        open_name = '%s.open' % __name__
        with patch(open_name, mock_open(read_data=self.test_sv_gtf_content_3), create=True):
            with open("dummy_filename.gtf") as test_gtf_file:
                svs_bed_table = read_sv_gtf(test_gtf_file, SvType.INV)
                self.assertEquals(len(svs_bed_table), 2)
                self.assertEquals(svs_bed_table.iloc[0, 1], 7184103)
                self.assertEquals(svs_bed_table.iloc[0, 2], 7576873)
                self.assertEquals(svs_bed_table.iloc[0, 5], "+")