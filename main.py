import argparse
import pandas as pd
from glypnirO.common import GlypnirO, process_tmt_pd_byonic, GlypnirOComponent, filter_with_U

parser = argparse.ArgumentParser(description="Automated workflow for processing and combining Byonic and PD output")

parser.add_argument("-i",
                    "-input_file",
                    type=str,
                    help="Filepath to experiment description file where each row has 4 columns, condition_id, "
                         "replicate_id, filename, and area_filename", dest="i")
parser.add_argument("-o",
                    "-output",
                    type=str,
                    help="Filepath to output", dest="o")
parser.add_argument("-s",
                    "--score_cutoff",
                    type=float,
                    help="Filter the data by specific Byonic cut-off quality score",
                    default=200, dest="s")

parser.add_argument("-t",
                    "--trust_byonic",
                    help="Enable site specific glycan position parsing from Byonic glycan assignment",
                    action="store_true", dest="trust_byonic")

parser.add_argument("-g",
                    "--get_uniprot",
                    help="Required internet connection. Enable parsing of UniProt ID from protein name and request "
                         "the original protein name from the UniProt databas.",
                    action="store_true", dest="get_uniprot")

parser.add_argument("-p",
                    "--parse_uniprot",
                    help="Attempt to parse UniProt accession ID using regular expression and use them as master id",
                    action="store_true", dest="parse_uniprot")

parser.add_argument("-d",
                    "--debug",
                    help="In conjunction to the final output, the script would also create debug files that contain the "
                         "unique PSM selected for calculation of the data in the final output",
                    action="store_true", dest="debug")

parser.add_argument("-m",
                    "--mode",
                    type=int,
                    help="Deciding the operating mode of GlypNirO\n"
                         "1 - normal pd + byonic dda"
                         "2 - pd + byonic tmt",
                    default=1,
                    dest="m")

if __name__ == "__main__":

    args = parser.parse_args()
    print("Input: " + args.i)
    print("Output: " + args.o)
    ex = GlypnirO(trust_byonic=args.trust_byonic, get_uniprot=args.get_uniprot, debug=args.debug,
                  parse_uniprot=args.parse_uniprot)
    if args.m == 1:
        for i, r in ex.add_batch_component(args.i, args.s):
            print(r)
        ex.process_components()
        result = ex.analyze_components()
        with pd.ExcelWriter(args.o) as writer:
            print("Writing Occupancy data to excel sheets.")
            if not result["Occupancy"].empty:
                result["Occupancy"].to_excel(writer, sheet_name="Unglyco_and_Glycoforms_Prop")
            if not result["Occupancy_Without_Proportion_U"].empty:
                result["Occupancy_Without_Proportion_U"].to_excel(writer, sheet_name="Unglyco_and_Glycoforms_Sep")
            if not result["Occupancy_With_U"].empty:
                result["Occupancy_With_U"].to_excel(writer, sheet_name="Unglycosylated")
            if not result["Glycoforms"].empty:
                result["Glycoforms"].to_excel(writer, sheet_name="Glycoforms")
        if ex.debug:
            for u in ex.unique_dict:
                with pd.ExcelWriter(args.o + "_" + u + ".xlsx") as writer:
                    df = pd.DataFrame(ex.unique_dict[u])
                    df.to_excel(writer, index=False)
    elif args.m == 2:
        data = pd.read_csv(args.i, sep="\t", encoding="utf-8")
        data, sample_info = process_tmt_pd_byonic(data)

        ex.uniprot_parsed_data = data[["Master Protein Accessions", "Protein Descriptions"]].rename(columns={"Master Protein Accessions": "Entry", "Protein Descriptions": "Protein names"})
        component = GlypnirOComponent(data, mode=args.m, trust_byonic=args.trust_byonic)
        component.process(mode=args.m, tmt_info=sample_info, protein_column="Master Protein Accessions", sequence_column="Annotated Sequence", glycans_column="Glycan composition", starting_position_column="Position in Protein")

        result = component.analyze(debug=args.debug, protein_column="Master Protein Accessions", glycans_column="Glycan composition", starting_position_column="Position in Protein", mode=args.m, tmt_info=sample_info, observed_mz_column="m/z [Da]")
        if args.debug:
            pd.DataFrame(component.unique_rows).to_csv(args.o + "_debug.txt", sep="\t", index=False)
        output1 = []
        output_without_u = []
        output_occupancy_no_calculation_u = []

        if not result.empty:
            for s in result.separate_result():
                print(s.df)
                a = s.to_summary(name="Raw", trust_byonic=args.trust_byonic)

                pro = s.calculate_proportion()
                b = s.to_summary(pro, "Proportion", trust_byonic=args.trust_byonic)
                condition = "None"
                for cond in sample_info:
                    if s.df.columns[-1] in sample_info[cond]:
                        condition = cond
                temp_df = ex._summary(a, b, add_protein=False, condition=condition, replicate=s.df.columns[-1])
                # print(temp_df)
                output1.append(temp_df)

                # proportion for glycoforms here are calculated without unglycosylated form.
                a_without_u = s.to_summary(name="Raw", trust_byonic=args.trust_byonic, occupancy=False)
                pro_without_u = s.calculate_proportion(occupancy=False)
                b_without_u = s.to_summary(pro_without_u, "Proportion", trust_byonic=args.trust_byonic,
                                                         occupancy=False)

                temp_df_without_u = ex._summary(a_without_u, b_without_u, add_protein=False, condition=condition, replicate=s.df.columns[-1])
                output_without_u.append(temp_df_without_u)

                temp_df_no_calculation_u = ex._summary(a, b_without_u, add_protein=False, condition=condition, replicate=s.df.columns[-1])
                output_occupancy_no_calculation_u.append(temp_df_no_calculation_u)

        result_occupancy = ex._summary_format(output1)
        result_occupancy_with_u = ex._summary_format(output1, filter_with_U, True)
        result_glycoform = ex._summary_format(output_without_u)
        tempdf_index_reset_result_occupancy_with_u = result_occupancy_with_u.reset_index()
        tempdf_index_reset_result_glycoform = result_glycoform.reset_index()
        result_occupancy_glycoform_sep = pd.concat(
            [tempdf_index_reset_result_glycoform, tempdf_index_reset_result_occupancy_with_u])
        # format the output with the correct column name for site specific or peptide level analysis
        if args.trust_byonic:
            result_occupancy_glycoform_sep = result_occupancy_glycoform_sep.set_index(["Protein", "Protein names",
                                                                                       # "Isoform",
                                                                                       "Glycosylated positions in "
                                                                                       "peptide",
                                                                                       "Glycans"])
            result_occupancy_glycoform_sep = result_occupancy_glycoform_sep.sort_index(
                level=["Protein", "Protein names",
                       # "Isoform",
                       "Glycosylated positions in peptide"])
        else:

            result_occupancy_glycoform_sep = result_occupancy_glycoform_sep.set_index(
                ["Protein", "Protein names",
                 # "Isoform",
                 "Position peptide N-terminus", "Peptides", "Glycans"])
            result_occupancy_glycoform_sep = result_occupancy_glycoform_sep.sort_index(
                level=["Protein", "Protein names",
                       # "Isoform",
                       "Position peptide N-terminus", "Peptides"])
        with pd.ExcelWriter(args.o) as writer:
            print("Writing Occupancy data to excel sheets.")
            if not result_occupancy.empty:
                result_occupancy.to_excel(writer, sheet_name="Unglyco_and_Glycoforms_Prop")
            if not result_occupancy_glycoform_sep.empty:
                result_occupancy_glycoform_sep.to_excel(writer, sheet_name="Unglyco_and_Glycoforms_Sep")
            if not result_occupancy_with_u.empty:
                result_occupancy_with_u.to_excel(writer, sheet_name="Unglycosylated")
            if not result_glycoform.empty:
                result_glycoform.to_excel(writer, sheet_name="Glycoforms")
    print("Completed.")

