import pandas as pd
import sys
import subprocess

def calculate_rpkm(counts_file, bed_file_HepG2_Ctrl,bed_file_K562_Ctrl,bed_file_HepG2_rep1,bed_file_K562_rep1,lengths_file):
    def get_bed_line_count(bed_file):
        try:
            line_count = int(subprocess.check_output(f"wc -l < {bed_file}", shell=True))
            return line_count
        except subprocess.CalledProcessError as e:
            print(f"Error counting lines in {bed_file}: {e}")
            sys.exit(1)
            
    total_reads_HepG2_ctrl = get_bed_line_count(bed_file_HepG2_Ctrl)
    total_reads_K562_ctrl = get_bed_line_count(bed_file_K562_Ctrl)
    total_reads_HepG2_rep1 = get_bed_line_count(bed_file_HepG2_rep1)
    total_reads_K562_rep1 = get_bed_line_count(bed_file_K562_rep1)
#    total_reads_HepG2_rep2 = get_bed_line_count(bed_file_HepG2_rep2)
#    total_reads_K562_rep2 = get_bed_line_count(bed_file_K562_rep2)
    lengths_df = pd.read_csv(lengths_file, sep='\t', header=None, names=['Gene Name', 'length'])
    



    


    counts_df = pd.read_csv(counts_file, sep="\t", header=None)
    counts_df.columns = ["Gene Name","HepG2_ctrl", "K562_ctrl",  "HepG2_rep1", "K562_rep1"]
    counts_df = counts_df.merge(lengths_df, on='Gene Name', how='left')
    counts_df.columns = ["Gene Name","HepG2_ctrl", "K562_ctrl",  "HepG2_rep1", "K562_rep1","length"]
    counts_df[["HepG2_ctrl", "K562_ctrl",  "HepG2_rep1", "K562_rep1","length"]] = counts_df[["HepG2_ctrl", "K562_ctrl",  "HepG2_rep1", "K562_rep1","length"]].apply(pd.to_numeric, errors='coerce')

    counts_df["RPKM_HepG2_ctrl"] = (counts_df["HepG2_ctrl"] / total_reads_HepG2_ctrl /counts_df["length"]) * 1e9
    counts_df["RPKM_K562_ctrl"] = (counts_df["K562_ctrl"] / total_reads_K562_ctrl /counts_df["length"]) * 1e9
    counts_df["RPKM_HepG2_rep1"] = (counts_df["HepG2_rep1"] / total_reads_HepG2_rep1 /counts_df["length"]) * 1e9
    counts_df["RPKM_K562_rep1"] = (counts_df["K562_rep1"] / total_reads_K562_rep1 /counts_df["length"]) * 1e9
#    counts_df["RPKM_HepG2_rep2"] = (counts_df["HepG2_rep2"] / total_reads_HepG2_rep2 /counts_df["length"]) * 1e9
#    counts_df["RPKM_K562_rep2"] = (counts_df["K562_rep2"] / total_reads_K562_rep2 /counts_df["length"]) * 1e9

#    counts_df["Avg_Rep1"] = (counts_df["RPKM_K562_rep1"] + counts_df["RPKM_HepG2_rep2"] + counts_df["RPKM_K562_rep2"]) / 3
#    nan_condition = counts_df["Avg_Rep1"].notna() 
#    counts_df = counts_df[nan_condition]
#    filter_condition_low_rpkm = (counts_df["Avg_Rep1"] >= 100)
#    counts_df = counts_df[filter_condition_low_rpkm]
#    counts_df = counts_df.drop(columns=["Avg_Rep1"])
    
    output_file = sys.argv[7]
    counts_df.to_csv(output_file, sep="\t", index=False)

    print(f"Calculations completed. Results saved in {output_file}")

if __name__ == "__main__":
    if len(sys.argv) !=8:
        print("Usage: python calculate_rpm_fc.py <counts_file> <bed_file_siCtrl> <bed_file_siRBP>")
        sys.exit(1)

    counts_file = sys.argv[1]
    bed_file_HepG2_Ctrl = sys.argv[2]
    bed_file_K562_Ctrl = sys.argv[3]    
    bed_file_HepG2_rep1 = sys.argv[4]
    bed_file_K562_rep1 = sys.argv[5]
#    bed_file_HepG2_rep2 = sys.argv[6]
#    bed_file_K562_rep2 = sys.argv[7]

    lengths_file = sys.argv[6]
    calculate_rpkm(counts_file, bed_file_HepG2_Ctrl,bed_file_K562_Ctrl,bed_file_HepG2_rep1,bed_file_K562_rep1,lengths_file)
