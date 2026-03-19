import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap



def merge_er_files(pattern_suffixes):
    rbp_files = {}
    for suffix in pattern_suffixes:
        for f in glob.glob(f"*{suffix}"):
            rbp_name = f.split('_')[0]
            if rbp_name not in rbp_files:
                rbp_files[rbp_name] = []
            rbp_files[rbp_name].append(f)
            
    if not rbp_files:
        print(f" No files: {pattern_suffixes}")
        return pd.DataFrame()
        
    merged_df = None
    
    for rbp, flist in rbp_files.items():
        rbp_df_list = []
        for f in flist:
            df = pd.read_csv(f, sep="\t")
            
            if 'log2FC_Mean' in df.columns:
                val_col = 'log2FC_Mean'
            elif 'log2FC_rep1' in df.columns:
                val_col = 'log2FC_rep1'
            else:
                continue
                
            df_sub = df[['Gene Name', val_col]].rename(columns={val_col: rbp})
            rbp_df_list.append(df_sub)
            
        rbp_combined = pd.concat(rbp_df_list, ignore_index=True).drop_duplicates(subset=['Gene Name'])
        
        if merged_df is None:
            merged_df = rbp_combined
        else:
            merged_df = pd.merge(merged_df, rbp_combined, on="Gene Name", how="outer")
            
    if merged_df is None or merged_df.empty:
         return pd.DataFrame()

    merged_df = merged_df.fillna(0)
    merged_df = merged_df.set_index("Gene Name")
    
    merged_df.rename(index={'U8': 'SNORD118'}, inplace=True)
    
    target_miRNAs = ['MIR21', 'MIR584', 'MIR644', 'MIR7705', 'MIR877', 'MIR1304', 'AL513534.2', 'MIR1296', 'MIR3691', 'MIR624']
    target_snoRNAs = [
        "SNORD96", "SCARNA14", "SNORD118", "SNORD3A", "SNORD3C", 
        "SNORD89", "SNORD101", "SNORD11B", "SNORD62", "SNORD11", 
        "SNORD20", "SNORD25", "SNORD58", "SNORD79", "SCARNA6", 
        "SNORD16", "SNORA75", "SNORD59"
    ]

    target_RNAs = target_miRNAs + target_snoRNAs
    merged_df = merged_df.reindex(index=target_RNAs).fillna(0) 
    
    target_rbps = [
        "SMNDC1", "TRA2A", "SNU13", "U1A", "LARP7", "NOLC1", "PTBP1"
    ]
    merged_df = merged_df.reindex(columns=target_rbps).fillna(0)
    
    return merged_df, len(target_miRNAs)

def draw_simple_heatmap(df, out_pdf, title_text, split_index=None):
    if df.empty:
        print(f"Skip{out_pdf}")
        return

    colors = ['#3182bd', '#ffffff', '#de2d26'] 
    cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)
    

    plt.figure(figsize=(5.5, 9))
    
    ax = sns.heatmap(
        df, 
        cmap=cmap, 
        vmin=-2.0, vmax=4.0,            
        center=0,                       
        linewidths=0.5, linecolor='white', 
        cbar_kws={'label': 'Log2FC', 'shrink': 0.8} 
    )
    

    if split_index:
        ax.axhline(split_index, color='black', linewidth=1.5, linestyle='--')
        

        mirna_center = split_index / 2
        snorna_center = split_index + (len(df) - split_index) / 2
        

        ax.text(-0.35, mirna_center, 'miRNA', 
                transform=ax.get_yaxis_transform(), 
                va='center', ha='right', fontsize=12, fontweight='bold', rotation=90, clip_on=False)
                
        ax.text(-0.35, snorna_center, 'snoRNA', 
                transform=ax.get_yaxis_transform(), 
                va='center', ha='right', fontsize=12, fontweight='bold', rotation=90, clip_on=False)

    plt.title(title_text, fontsize=14, fontweight='bold', pad=20)
    
    plt.xticks(rotation=45, ha="right", fontsize=11)
    plt.yticks(rotation=0, fontsize=11)
    
    ax.set_ylabel('')
    ax.set_xlabel('')
    
    plt.savefig(out_pdf, format='pdf', bbox_inches='tight', dpi=600)
    print(f"SAVED: {out_pdf}")
    
    plt.close()

if __name__ == "__main__":
    
    print("Processing Mature RNA Data...")
    df_sno_mi, split_idx = merge_er_files(["_miRNA_ER.txt","_snoRNA_ER.txt"])
    draw_simple_heatmap(
        df_sno_mi, 
        "Figure3A_Target_RNA.png", 
        "RBP eCLIP Enrichment on snRNA-interacting snoRNAs", 
        split_index=split_idx
    )

    print("\nProcessing  Data...")
    df_pre_mi, split_idx = merge_er_files(["_miRNA_5'3'_ER.txt","_snoRNA_5'3'_ER.txt" ])
    draw_simple_heatmap(
        df_pre_mi, 
        "Figure3B_Target_RNA5'3'.png", 
        "RBP eCLIP Enrichment on snRNA-interacting snoRNA 5' and 3' ends", 
        split_index=split_idx
    )