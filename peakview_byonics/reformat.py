import argparse
import pandas as pd
from sequal.sequence import Sequence
from sequal.resources import proton
import os
parser = argparse.ArgumentParser(description="Utility script for automated formatting and generation of GlypNirO input from Byonic and PeakView files.")
parser.add_argument("-b",
                    "-byonic",
                    type=str,
                    help="Filepath to Byonic output xlsx file."
                         "replicate_id, filename, and area_filename", dest="b")
parser.add_argument("-p",
                    "-peakview",
                    type=str,
                    help="Filepath to PeakView peptide output in xlsx format", dest="p")

parser.add_argument("-o",
                    "-output",
                    type=str,
                    help="Filepath to output", dest="o")
#byonic_file = r"C:\Users\localadmin\PycharmProjects\glypniro2\Yeast_Byonic.unique_EDK.xlsx"
#peakview_file = r"C:\Users\localadmin\PycharmProjects\glypniro2\20200904_Ed_beers_beerdb_SWATH+yeast+barley_glycans_IonLib_SWATH - Peptides_Recalc.xlsx"
#dia_library_file = r"C:\Users\localadmin\PycharmProjects\glypniro2\20190605_Combined_Byonic_Grape_Yeast_ProteinPilot_Pau4_2_No_LOC100247730_Added_KPYK2.txt"
#output = r"description_yeast.xlsx"

if __name__ == "__main__":
    args = parser.parse_args()

    byonic = pd.read_excel(args.b)
    peakview = pd.read_excel(args.p)
    result = []
    for i, _ in byonic.groupby(["Peptide\n< ProteinMetrics Confidential >", "Glycans\nNHFAGNa", "Protein Name", "Starting\nposition", "Calc.\nmass (M+H)", "z"]):
        seq = Sequence(i[0])
        flank_removed = str(seq)[2:-2]
        stripped_seq = seq.to_stripped_string()[2:-2]
        calculated_mz = (i[4] - proton + proton*i[5])/i[5]
        result.append(list(i[:-1]) + [flank_removed, stripped_seq, round(calculated_mz, 3)])

    df = pd.DataFrame(result, columns=["Peptide", "Glycans\nNHFAGNa", "Protein Name", "Starting\nposition", "Calc.\nmass (M+H)", "Flank Removed", "Stripped Sequence", "mz"])

    description = []
    stripped_seq_dict = {}
    scan = 0
    for c in peakview.columns[5:]:
        area = []
        byonic_glycan = []
        condition, replicate = c.rsplit("_", 1)
        filepath = os.path.join(os.path.dirname(os.path.abspath(__file__)), c)
        description.append([condition, replicate, filepath+".xlsx", filepath+".txt"])
        columns = list(peakview.columns[:5]) + [c]
        df_result = peakview[columns]

        for i, r in df_result.iterrows():
            scan += 1
            key = r["Peptide"]
            data = {}
            if key not in stripped_seq_dict:
                seq = Sequence("-.{}.-".format(r["Peptide"]))
                for aa in range(len(seq)):
                    if seq[aa].mods:
                        try:
                            m = float(seq[aa].mods[0].value)
                            round_3 = round(m, 3)
                            print(m, round_3)
                            if round_3 >= 0:
                                before, after = str(round_3).split(".")
                                if len(after) < 3:
                                    for d in range(3-len(after)):
                                        after += "0"
                                seq[aa].mods[0].value = "+{}.{}".format(before, after)
                            else:
                                seq[aa].mods[0].value = str(round_3)

                        except ValueError:
                            pass

                str_seq = str(seq)[2:-2]

                stripped_seq = seq.to_stripped_string()[2:-2]

                sub_df = df[df["Stripped Sequence"] == stripped_seq]
                if not sub_df.empty:

                    sub_sub_df = sub_df[sub_df["Flank Removed"] == str_seq]

                    if not sub_sub_df.empty:

                        data = {
                            "Peptide": "{}".format(sub_sub_df["Peptide"].values[0]),
                            "Glycans\nNHFAGNa": sub_sub_df["Glycans\nNHFAGNa"].values[0],
                            "Protein Name": sub_sub_df["Protein Name"].values[0],
                            "Starting\nposition": sub_sub_df["Starting\nposition"].values[0],
                            "Calc.\nmass (M+H)":  sub_sub_df["Calc.\nmass (M+H)"].values[0]}
                        stripped_seq_dict[key] = data

                    else:
                        prec_mz = round(r["Precursor MZ"], 3)
                        prec_sub_df = sub_df[sub_df["mz"] == prec_mz]
                        if not prec_sub_df.empty:
                            print(prec_sub_df)
                            data = {
                                "Peptide": "{}".format(prec_sub_df["Peptide"].values[0]),
                                "Glycans\nNHFAGNa": prec_sub_df["Glycans\nNHFAGNa"].values[0],
                                "Protein Name": prec_sub_df["Protein Name"].values[0],
                                "Starting\nposition": prec_sub_df["Starting\nposition"].values[0],
                                "Calc.\nmass (M+H)": prec_sub_df["Calc.\nmass (M+H)"].values[0]}
                            stripped_seq_dict[key] = data
                        else:
                            data = {
                            "Peptide": "{}{}{}".format(sub_df["Peptide"].values[0][:2], r["Peptide"], sub_df["Peptide"].values[0][-2:]),
                            "Glycans\nNHFAGNa": "",
                            "Protein Name": sub_df["Protein Name"].values[0],
                            "Starting\nposition": sub_df["Starting\nposition"].values[0],
                            "Calc.\nmass (M+H)": r["Precursor MZ"]}
                            stripped_seq_dict[key] = data
                else:
                    data = {"Peptide": "{}".format(str(seq)), "Glycans\nNHFAGNa": "", "Protein Name": r["Protein"], "Starting\nposition": "", "Calc.\nmass (M+H)": r["Precursor MZ"]}
                    stripped_seq_dict[key] = data
            else:
                data = stripped_seq_dict[key]

            byonic_glycan.append([data["Protein Name"], data["Glycans\nNHFAGNa"], data["Peptide"], 1, "scan={}".format(str(scan)), data["Starting\nposition"], data["Calc.\nmass (M+H)"]])
            area.append([scan, r[c]])
        byonic_df = pd.DataFrame(byonic_glycan, columns=["Protein Name", "Glycans\nNHFAGNa", "Peptide\n< ProteinMetrics Confidential >", "Score", "Scan #", "Starting\nposition", "Calc.\nmass (M+H)"])
        with pd.ExcelWriter(filepath+".xlsx") as writer:
            byonic_df.to_excel(writer, sheet_name="Spectra", index=False)
        area_df = pd.DataFrame(area, columns=["First Scan", "Area"])
        area_df.to_csv(filepath+".txt", sep="\t", index=False)
    description = pd.DataFrame(description, columns=["condition_id", "replicate_id", "filename", "area_filename"])
    with pd.ExcelWriter(args.o) as writer:
        description.to_excel(writer, index=False)
