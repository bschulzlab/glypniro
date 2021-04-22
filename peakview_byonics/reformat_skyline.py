import re

import pandas as pd
import argparse
from sequal.sequence import Sequence
import pathlib

glycan_column = "Glycans\nNHFAGNa"
peptide_column = "Peptide\n< ProteinMetrics Confidential >"
scan_number_regex = re.compile("scan=(\d+)")

parser = argparse.ArgumentParser(description="Utility script for automated formatting and generation of GlypNirO input from a specific Skyline output.")
parser.add_argument("-b",
                    "-byonic",
                    type=str,
                    help="Filepath to Byonic output xlsx file."
                         "replicate_id, filename, and area_filename", dest="b")
parser.add_argument("-s",
                    "-skyline",
                    type=str,
                    help="Filepath to Skyline peptide output in xlsx format", dest="s")

parser.add_argument("-o",
                    "-output",
                    type=str,
                    help="Filepath to output folder", dest="o")

replicate_regex = re.compile("(\w+)(\d+)$")

if __name__ == "__main__":
    work = []
    args = parser.parse_args()
    df = pd.read_csv(args.s)
    byonic = pd.read_excel(args.b)
    split_data = df["Peptide"].str.split("_", expand=True)
    df[peptide_column] = split_data[0]
    df[glycan_column] = split_data[1]
    df["z"] = split_data[2]
    byonic.drop([glycan_column, "z"], axis=1, inplace=True)
    df = df.merge(byonic, on=peptide_column)
    current_condition = ""
    current_replicate = ""
    matches = df["Replicate"].str.extractall(replicate_regex).reset_index()
    #scan_number = df["Scan #"].str.extractall(scan_number_regex).reset_index()
    df["Condition id"] = matches[0]
    df["Replicate id"] = matches[1]
    df["Area"] = df["Normalized Area"]
    output_folder = args.o
    pathlib.Path(output_folder).mkdir(parents=True, exist_ok=True)
    for g, d in df.groupby(["Condition id", "Replicate id"]):
        pd_file = str(pathlib.Path(output_folder).joinpath("_".join(g) + ".txt"))
        byonic_file = str(pathlib.Path(output_folder).joinpath("_".join(g) + ".xlsx"))
        scan_pd = []
        scan_byonic = []
        for i in range(d.shape[0]):
            scan_pd.append(str(i))
            scan_byonic.append("scan="+str(i))
        with pd.ExcelWriter(byonic_file) as writer:
            d["Scan #"] = pd.Series(scan_byonic, index=d.index)
            d["First Scan"] = pd.Series(scan_pd, index=d.index)
            d[["Protein Name", "Glycans\nNHFAGNa", "Peptide\n< ProteinMetrics Confidential >", "Score", "Scan #", "Starting\nposition", "Calc.\nmass (M+H)"]].to_excel(writer, sheet_name="Spectra", index=False)
            d[["First Scan", "Area"]].to_csv(pd_file, sep="\t", index=False)
            work.append([g[0], g[1], byonic_file, pd_file])

    work = pd.DataFrame(work, columns=["condition_id", "replicate_id", "filename", "area_filename"])
    with pd.ExcelWriter(str(pathlib.Path(output_folder).joinpath("work.xlsx"))) as writer:
        work.to_excel(writer, index=False)