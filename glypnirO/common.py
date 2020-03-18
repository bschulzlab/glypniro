from copy import deepcopy
import numpy as np
import pandas as pd
import re

from sequal.sequence import Sequence

sequence_column_name = "Peptide\n< ProteinMetrics Confidential >"
glycans_column_name = "Glycans\nNHFAGNa"
modifications_column_name = "Modification Type(s)"
observed_mz = "Calc.\nmass (M+H)"
protein_column_name = "Protein Name"
rt = "Scan Time"
selected_aa = {"N", "S", "T"}

regex_glycan_number_pattern = "\d+"
glycan_number_regex = re.compile(regex_glycan_number_pattern)
regex_pattern = "\.[\[\]\w\.\+\-]*\."
sequence_regex = re.compile(regex_pattern)


def get_mod_value(amino_acid):
    if amino_acid.mods:
        if amino_acid.mods[0].value.startswith("+"):
            return float(amino_acid.mods[0].value[1:])
        else:
            return -float(amino_acid.mods[0].value[1:])
    else:
        return 0


def load_fasta(fasta_file_path, selected=None, selected_prefix=""):
    with open(fasta_file_path, "rt") as fasta_file:
        result = {}
        current_seq = ""
        for line in fasta_file:
            line = line.strip()

            if line.startswith(">"):
                if selected:
                    if selected_prefix + line[1:] in selected:
                        result[line[1:]] = ""
                        current_seq = line[1:]
                else:
                    result[line[1:]] = ""
                    current_seq = line[1:]
            else:
                result[current_seq] += line
        return result


class Result:
    def __init__(self, df):
        self.df = df
        self.empty = df.empty

    def calculate_proportion(self):
        df = self.df.copy()
        for i in df.columns:
            df[i] = df[i]/df[i].sum()
        return df

    def to_summary(self, df=None, name=""):
        if df is None:
            df = self.df
        temp = df.unstack()
        temp.index.names = ["Peptides", "Glycans"]
        temp.name = name
        return temp



class GlypnirOComponent:
    def __init__(self, filename, area_filename, replicate_id, condition_id, protein_name, minimum_area=0):
        if type(filename) == pd.DataFrame:
            data = filename.copy()
        else:
            data = pd.read_excel(filename, sheet_name="Spectra")
        if type(area_filename) == pd.DataFrame:
            file_with_area = area_filename
        else:
            if area_filename.endswith("xlsx"):
                file_with_area = pd.read_excel(area_filename)
            else:
                file_with_area = pd.read_csv(area_filename, sep="\t")
        data["Scan number"] = pd.to_numeric(data["Scan #"].str.extract("scan=(\d+)", expand=False))
        data = pd.merge(data, file_with_area, left_on="Scan number", right_on="First Scan")
        self.protein_name = protein_name
        self.data = data.sort_values(by=['Area'], ascending=False)
        self.replicate_id = replicate_id
        self.condition_id = condition_id
        self.data = data[data["Area"].notnull()]
        self.data = data[(data["Area"] >= minimum_area) & (data["Protein Name"] == ">"+protein_name)]
        if len(self.data.index) > 0:
            self.empty = False
        else:
            self.empty = True
        self.row_to_glycans = {}
        self.glycan_to_row = {}
        self.analysis = ""

    def process(self, motif, fasta_library=None, expand_window=0, analysis="N-glycan"):
        self.analysis = analysis
        entries_number = len(self.data.index)
        if analysis == "N-glycan":
            self.data["total_number_of_asn"] = pd.Series([0]*entries_number, index=self.data.index, dtype=int)
            self.data["total_number_of_n-linked_sequon"] = pd.Series([0]*entries_number, index=self.data.index, dtype=int)
            self.data["total_number_of_hexnac"] = pd.Series([0]*entries_number, index=self.data.index, dtype=int)
            self.data["total_number_of_deamidation"] = pd.Series([0]*entries_number, index=self.data.index, dtype=int)
            self.data["total_number_of_modded_asn"] = pd.Series([0]*entries_number, index=self.data.index, dtype=int)
            self.data["total_number_of_unmodded_asn"] = pd.Series([0] * entries_number, index=self.data.index, dtype=int)
        elif analysis == "O-glycan":
            self.data["total_number_of_hex"] = pd.Series([0]*entries_number, index=self.data.index, dtype=int)
            self.data["total_number_of_modded_ser_thr"] = pd.Series([0]*entries_number, index=self.data.index, dtype=int)
            self.data["total_number_of_unmodded_ser_or_thr"] = pd.Series([0]*entries_number, index=self.data.index, dtype=int)
            self.data["o_glycosylation_status"] = pd.Series([False]*entries_number, index=self.data.index, dtype=bool)
        for i, r in self.data.iterrows():

            search = sequence_regex.search(r[sequence_column_name])
            seq = Sequence(search.group(0))
            stripped_seq = seq.to_stripped_string()
            modifications = {}
            if pd.notnull(r[modifications_column_name]):
                for mod in r[modifications_column_name].split(","):

                    number = 1
                    if "*" in mod:
                        m = mod.split("*")
                        minimod = Sequence(m[0].strip())
                        number = int(m[1].strip())
                    else:
                        minimod = Sequence(mod.strip())
                    if minimod[0].mods[0].value not in modifications:
                        modifications[minimod[0].mods[0].value] = {}
                    modifications[minimod[0].mods[0].value][minimod[0].value] = {"mod": deepcopy(minimod[0].mods[0]),
                                                                                 "number": number}

                    if minimod[0].value == "N":
                        if analysis == "N-glycan":
                            if minimod[0].mods[0].value == 1:
                                self.data.at[i, "total_number_of_deamidation"] += number
                            self.data.at[i, "total_number_of_modded_asn"] += number
                    elif minimod[0].value in "ST":
                        if analysis == "O-glycan":
                            self.data.at[i, "total_number_of_modded_ser_thr"] += number
            glycans = []
            if pd.notnull(r[glycans_column_name]):
                glycans = r[glycans_column_name].split(",")
            if search:

                self.data.at[i, "stripped_seq"] = stripped_seq.rstrip(".").lstrip(".")

                origin_seq = fasta_library[self.protein_name].index(self.data.at[i, "stripped_seq"])
                self.data.at[i, "origin_start"] = origin_seq + 1

                motifs = [match for match in seq.find_with_regex(motif, ignore=seq.gaps())]
                if expand_window:
                    expanded_window = Sequence(fasta_library[self.protein_name][origin_seq: origin_seq + len(self.data.at[i, "stripped_seq"]) + expand_window])
                    expanded_window_motifs = [match for match in expanded_window.find_with_regex(motif, ignore=seq.gaps())]
                    origin_map = [i.start + origin_seq - 1 for i in expanded_window_motifs]
                    if len(expanded_window_motifs) > len(motifs):
                        self.data.at[i, "expanded_motif"] = str(expanded_window[expanded_window_motifs[-1]])
                    self.data.at[i, "expanded_aa"] = str(expanded_window[-expand_window:])

                else:
                    origin_map = [i.start + origin_seq - 1 for i in motifs]

                if analysis == "N-glycan":
                    self.data.at[i, "total_number_of_asn"] = seq.count("N", 0, len(seq))
                    if expand_window:
                        self.data.at[i, "total_number_of_n-linked_sequon"] = len(expanded_window_motifs)
                    else:
                        self.data.at[i, "total_number_of_n-linked_sequon"] = len(motifs)
                    self.data.at[i, "total_number_of_unmodded_asn"] = self.data.at[i, "total_number_of_asn"] - self.data.at[i, "total_number_of_modded_asn"]
                elif analysis == "O-glycan":
                    self.data.at[i, "total_number_of_ser_thr"] = seq.count("S", 0, len(seq)) + seq.count("T", 0, len(seq))
                    self.data.at[i, "total_number_of_unmodded_ser_or_thr"] = self.data.at[i, "total_number_of_modded_ser_thr"] - self.data.at[i, "total_number_of_modded_ser_thr"]


                current_glycan = 0
                max_glycans = len(glycans)
                if max_glycans:
                    self.row_to_glycans[i] = np.sort(glycans)
                    for g in glycans:
                        self.glycan_to_row[g] = i
                for aa in range(1, len(seq) - 1):
                    if seq[aa].mods:
                        round_mod_value = round(float(seq[aa].mods[0].value))
                        mod_value = seq[aa].mods[0].value[0] + str(round_mod_value)
                        if mod_value in modifications:
                            if seq[aa].value in "ST" and analysis == "O-glycan":
                                if round_mod_value == 80:
                                    continue

                            if seq[aa].value in modifications[mod_value]:
                                if seq[aa].value == "N" and round_mod_value == 1:
                                    seq[aa].extra = "Deamidated"

                                if modifications[mod_value][seq[aa].value]['number'] > 0:
                                    modifications[mod_value][seq[aa].value]['number'] -= 1
                                    seq[aa].mods[0].mass = float(seq[aa].mods[0].value)

                                    if max_glycans and current_glycan != max_glycans:

                                        seq[aa].mods[0].value = glycans[current_glycan]
                                        seq[aa].extra = "Glycosylated"

                                        if seq[aa].value == "N":
                                            if analysis == "N-glycan":
                                                if "hexnac" in glycans[current_glycan].lower():
                                                    self.data.at[i, "total_number_of_hexnac"] += 1

                                        elif seq[aa].value in "ST":
                                            if analysis == "O-glycan":
                                                self.data.at[i, "total_number_of_hex"] += 1

                                        current_glycan += 1
                                        #if current_glycan == max_glycans:
                                            #break

                glycosylation_count = 1
                for n in origin_map:
                    position = "{}_position".format(str(glycosylation_count))
                    self.data.at[i, position] = seq[n-origin_seq+1].value + str(
                        n + 1)

                    if seq[n-origin_seq+1].extra == "Glycosylated":
                        self.data.at[i, position + "_match"] = "H"
                    else:
                        self.data.at[i, position + "_match"] = "U"

                    if analysis == "N-glycan":
                        if self.data.at[i, "total_number_of_n-linked_sequon"] != self.data.at[i, "total_number_of_hexnac"]:
                            if seq[n-origin_seq+1].extra == "Deamidated":
                                if self.data.at[i, "total_number_of_hexnac"] > 0:
                                    self.data.at[i, position + "_match"] = "D/H"
                                    if self.data.at[i, "total_number_of_unmodded_asn"] > 0:
                                        self.data.at[i, position + "_match"] = "D/H/U"
                                else:
                                    self.data.at[i, position + "_match"] = "D"
                            else:
                                if self.data.at[i, "total_number_of_hexnac"] > 0:
                                    if self.data.at[i, "total_number_of_deamidation"] == 0:
                                        self.data.at[i, position + "_match"] = "H"
                                    else:
                                        self.data.at[i, position + "_match"] ="D/H"
                                    if self.data.at[i, "total_number_of_unmodded_asn"] > 0:
                                        self.data.at[i, position + "_match"] = "D/H/U"
                            if not seq[n-origin_seq+1].extra:
                                if self.data.at[i, "total_number_of_hexnac"] > 0 and self.data.at[i, "total_number_of_deamidation"]> 0:
                                    self.data.at[i, position + "_match"] = "D/H"
                                    if self.data.at[i, "total_number_of_unmodded_asn"] > 0:
                                        self.data.at[i, position + "_match"] = "D/H/U"
                                elif self.data.at[i, "total_number_of_hexnac"] > 0:
                                    self.data.at[i, position + "_match"] = "H"
                                    if self.data.at[i, "total_number_of_unmodded_asn"] > 0:
                                        self.data.at[i, position + "_match"] = "D/H/U"
                                else:
                                    self.data.at[i, position + "_match"] = "U"
                    glycosylation_count += 1
                if analysis == "O-glycan":
                    if max_glycans:
                        self.data.at[i, "o_glycosylation_status"] = True

    def analyze(self, max_sites = 1):
        temp = self.data.sort_values(["Area", "Score"], ascending=False)
        if self.analysis == "N-glycan":
            temp = temp[temp["total_number_of_n-linked_sequon"] <= max_sites]

        result = {}
        result_total = {}
        temp[glycans_column_name] = temp[glycans_column_name].fillna("None")
        for i, g in temp.groupby([sequence_column_name, "z", glycans_column_name]):
            unique_row = g.loc[g["Area"].idxmax()]
            if unique_row["stripped_seq"] not in result:
                result[unique_row["stripped_seq"]] = {}
                result_total[unique_row["stripped_seq"]] = 0
            if i[2] == "None":
                if "U" not in result[unique_row["stripped_seq"]]:
                    result[unique_row["stripped_seq"]]["U"] = 0
                result[unique_row["stripped_seq"]]["U"] += unique_row["Area"]

            else:
                for gly in self.row_to_glycans[unique_row.name]:
                    if gly not in result[unique_row["stripped_seq"]]:
                        result[unique_row["stripped_seq"]][gly] = 0
                    result[unique_row["stripped_seq"]][gly] += unique_row["Area"]

        df = pd.DataFrame(result)
        return Result(df)


class GlypnirO:
    def __init__(self, fasta_filename):
        self.fasta = load_fasta(fasta_filename)
        self.components = None

    def add_component(self, filename, area_filename, replicate_id, sample_id):
        component = GlypnirOComponent(filename, area_filename, replicate_id, sample_id)

    def add_batch_component(self, component_list, minimum_area, protein=None):
        self.components = pd.DataFrame(component_list)

        if protein is not None:
            self.components["Protein"] = pd.Series([protein]*len(self.components.index), index=self.components.index)
            for i, r in self.components.iterrows():
                comp = GlypnirOComponent(r["filename"], r["area_filename"], r["replicate_id"], condition_id=r["condition_id"], protein_name=protein, minimum_area=minimum_area)
                self.components.at[i, "component"] = comp
                print("{} - {}, {} peptides has been successfully loaded".format(r["condition_id"], r["replicate_id"], str(len(comp.data.index))))

        else:
            components = []
            for i, r in self.components.iterrows():
                data = pd.read_excel(r["filename"], sheet_name="Spectra")
                if r["area_filename"].endswith("xlsx"):
                    file_with_area = pd.read_excel(r["area_filename"])
                else:
                    file_with_area = pd.read_csv(r["area_filename"], sep="\t")
                for index, g in data.groupby([protein_column_name]):

                    u = index[1:]
                    if not u.startswith("Reverse") and not u.endswith("(Common contaminant protein)"):
                        comp = GlypnirOComponent(g, file_with_area, r["replicate_id"],
                                                 condition_id=r["condition_id"], protein_name=u,
                                                 minimum_area=minimum_area)
                        if not comp.empty:
                            components.append({"filename": r["filename"], "area_filename": r["area_filename"], "condition_id": r["condition_id"], "replicate_id": r["replicate_id"], "Protein": u, "component": comp})
                print(
                    "{} - {} peptides has been successfully loaded".format(r["condition_id"],
                                                                                      r["replicate_id"]))
            self.components = pd.DataFrame(components, columns=list(self.components.columns) + ["component", "Protein"])



    def process_components(self, motif, analysis):
        for i, r in self.components.iterrows():
            # print("Processing {} - {} {} for {}".format(r["condition_id"], r["replicate_id"], r["Protein"], analysis))
            r["component"].process(motif, self.fasta, analysis=analysis)

    def analyze_components(self):
        # template = self.components[["Protein", "condition_id", "replicate_id"]].sort_values(["Protein", "condition_id", "replicate_id"])
        # template["label"] = pd.Series(["Raw"]*len(template.index), index=template.index)
        # template_proportion = template.copy()
        # template_proportion["label"] = pd.Series(["Proportion"]*len(template_proportion.index), index=template_proportion.index)
        #
        # index = pd.MultiIndex.from_frame(pd.concat([template, template_proportion], ignore_index=True)[["Protein", "label", "condition_id", "replicate_id"]])
        result = []
        for i, r in self.components.iterrows():
            # print(r["Protein"], r["condition_id"], r["replicate_id"])
            analysis = r["component"].analyze()
            if not analysis.empty:

                a = analysis.to_summary(name="Raw")
                pro = analysis.calculate_proportion()
                b = analysis.to_summary(pro, "Proportion")
                temp_df = pd.concat([a, b], axis=1)
                temp_df.columns.name = "Label"
                temp_df = temp_df.stack()
                lc = [temp_df]
                for c in ["Protein", "condition_id", "replicate_id"]:
                    lc.append(pd.Series([r[c]]*len(temp_df.index), index=temp_df.index, name=c))
                temp_df = pd.concat(lc, axis=1)
                temp_df = temp_df.reset_index()
                result.append(temp_df)

        result = pd.concat(result)
        result = result.reset_index(drop=True)
        groups = result.groupby(by=["Protein", "Peptides"])
        result = groups.filter(lambda x: (len(x["Glycans"].unique()) > 1) or ("U" not in x["Glycans"].unique()))
        result = result.set_index(["Label", "condition_id", "replicate_id", "Protein", "Peptides", "Glycans"])
        result = result.unstack(["Label", "condition_id", "replicate_id"])
        result = result.sort_index(level=["Label", "condition_id", "replicate_id"], axis=1)
        result = result.sort_index(level=["Peptides"])
        result.columns = result.columns.droplevel()
        # result = result.stack("Protein")
        # result = result.swaplevel("Protein", "Peptides")
        # result = result.swaplevel("Glycans", "Peptides")
        return result

