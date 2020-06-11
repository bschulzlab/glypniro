import unittest
from glypnirO.common import GlypnirOComponent, GlypnirO, load_fasta
import pandas as pd

o_glycan_file = r"C:\Users\localadmin\PycharmProjects\glypnirO\Olink_10_20ppm_TN_CSF_062617_04_A_R2.raw_Byonic.xlsx"
n_glycan_file = r"C:\Users\localadmin\PycharmProjects\glypnirO\Nlink_10_20ppm_TN_CSF_062617_04_A_R2.raw_Byonic.xlsx"
area_file = r"C:\Users\localadmin\PycharmProjects\glypnirO\TN_CSF_062617_04_MSnSpectrumInfo_A_R2.xlsx"
fasta_file = r"C:\Users\localadmin\PycharmProjects\glypnirO\Schulz_CP_Homo_sapiens_proteomeUP000005640reveiwed_Download_20180420_proteins_20303.fasta"

job = [
    {
        "filename": r"C:\Users\localadmin\PycharmProjects\glypnirO\data\70_Olink_10_20ppm_TN_CSF_062617_02_H_R1.raw_Byonic.xlsx",
        "area_filename": r"C:\Users\localadmin\PycharmProjects\glypnirO\data\TN_CSF_062617_02_MSnSpectrumInfo.txt",
        "replicate_id": r"R1",
        "condition_id": r"H"
    },
    {
        "filename": r"C:\Users\localadmin\PycharmProjects\glypnirO\data\70_Olink_10_20ppm_TN_CSF_062617_03_A_R1.raw_Byonic.xlsx",
        "area_filename": r"C:\Users\localadmin\PycharmProjects\glypnirO\data\TN_CSF_062617_03_MSnSpectrumInfo.txt",
        "replicate_id": r"R1",
        "condition_id": r"A"
    },
    {
        "filename": r"C:\Users\localadmin\PycharmProjects\glypnirO\data\70_Olink_10_20ppm_TN_CSF_062617_04_A_R2.raw_Byonic.xlsx",
        "area_filename": r"C:\Users\localadmin\PycharmProjects\glypnirO\data\TN_CSF_062617_04_MSnSpectrumInfo.txt",
        "replicate_id": r"R2",
        "condition_id": r"A"
    },
    {
        "filename": r"C:\Users\localadmin\PycharmProjects\glypnirO\data\70_Olink_10_20ppm_TN_CSF_062617_05_H_R2.raw_Byonic.xlsx",
        "area_filename": r"C:\Users\localadmin\PycharmProjects\glypnirO\data\TN_CSF_062617_05_H_R2_MSnSpectrumInfo.txt",
        "replicate_id": r"R2",
        "condition_id": r"H"
    },
    {
        "filename": r"C:\Users\localadmin\PycharmProjects\glypnirO\data\70_Olink_10_20ppm_TN_CSF_062617_06_A_R3.raw_Byonic.xlsx",
        "area_filename": r"C:\Users\localadmin\PycharmProjects\glypnirO\data\TN_CSF_062617_06_A_R3_MSnSpectrumInfo.txt",
        "replicate_id": r"R3",
        "condition_id": r"A"
    },
    {
        "filename": r"C:\Users\localadmin\PycharmProjects\glypnirO\data\70_Olink_10_20ppm_TN_CSF_062617_07_A_R4.raw_Byonic.xlsx",
        "area_filename": r"C:\Users\localadmin\PycharmProjects\glypnirO\data\TN_CSF_062617_07_A_R4_MSnSpectrumInfo.txt",
        "replicate_id": r"R4",
        "condition_id": r"A"
    }
]

job_N = [
    {
        "filename": r"C:\Users\localadmin\PycharmProjects\glypnirO\data\Nlink_Ragged_10_20ppm_TN_CSF_062617_02_H_R1.raw_Byonic.xlsx",
        "area_filename": r"C:\Users\localadmin\PycharmProjects\glypnirO\data\TN_CSF_062617_02_MSnSpectrumInfo.txt",
        "replicate_id": r"R1",
        "condition_id": r"H"
    },
    {
        "filename": r"C:\Users\localadmin\PycharmProjects\glypnirO\data\Nlink_Ragged_10_20ppm_TN_CSF_062617_03_A_R1.raw_Byonic.xlsx",
        "area_filename": r"C:\Users\localadmin\PycharmProjects\glypnirO\data\TN_CSF_062617_03_MSnSpectrumInfo.txt",
        "replicate_id": r"R1",
        "condition_id": r"A"
    }
]

class CommonTest(unittest.TestCase):
    def test_load_fasta(self):
        a = load_fasta(fasta_file)
        self.assertEqual(len(a), 305)


class GlypnirOComponentCase(unittest.TestCase):
    def test_init(self):
        a = GlypnirOComponent(o_glycan_file, area_file, replicate_id="test_replicate1", condition_id="test_sample1", protein_name="sp|P01024|CO3_HUMAN Complement C3 OS=Homo sapiens OX=9606 GN=C3 PE=1 SV=2")

    def test_analyze(self):
        a = GlypnirOComponent(o_glycan_file, area_file, replicate_id="test_replicate1", condition_id="test_sample1", protein_name="sp|P02649|APOE_HUMAN Apolipoprotein E OS=Homo sapiens OX=9606 GN=APOE PE=1 SV=1")
        fasta_library = load_fasta(fasta_file)
        a.process("[S|T]", fasta_library, analysis="O-glycan")
        a.data.to_csv("test.csv")
        r = a.analyze()

    def test_n_analyze(self):
        a = GlypnirOComponent(n_glycan_file, area_file, replicate_id="test_replicate1", condition_id="test_sample1",
                              protein_name="sp|P02649|APOE_HUMAN Apolipoprotein E OS=Homo sapiens OX=9606 GN=APOE "
                                           "PE=1 SV=1")
        fasta_library = load_fasta(fasta_file)
        a.process("(?=(N[^PX][ST]))", fasta_library, analysis="N-glycan")


class GlynirOCase(unittest.TestCase):
    def test_analyze(self):
        a = GlypnirO(fasta_file)
        a.add_batch_component(job, 0, "sp|P02649|APOE_HUMAN Apolipoprotein E OS=Homo sapiens OX=9606 GN=APOE PE=1 SV=1")
        a.process_components("[S|T]", "O-glycan")
        a.analyze_components().to_csv("test_analyze.csv")

    def test_analyze_all(self):
        a = GlypnirO(fasta_file)
        a.add_batch_component(job, 0)
        a.process_components("[S|T]", "O-glycan")
        result = a.analyze_components("O-glycan")
        with pd.ExcelWriter("test_analyze.xlsx") as writer:
            result.to_excel(writer)


class GlynirNCase(unittest.TestCase):
    def test_analyze(self):
        a = GlypnirO(fasta_file)
        a.add_batch_component(job_N, 0)
        a.process_components("(?=(N[^PX][ST]))", "N-glycan")
        result = a.analyze_components("N-glycan")
        with pd.ExcelWriter("test_analyze_N.xlsx") as writer:
            result.to_excel(writer)


if __name__ == '__main__':
    unittest.main()
