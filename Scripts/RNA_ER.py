import glob
import pandas as pd

def find_depleted_mirnas(min_rbps=5):

    files = glob.glob("*_miRNA_ER.txt")
    
    if not files:
        print("No *_miRNA_ER.txt, Check.")
        return pd.DataFrame()
        
    merged_df = None
    

    for f in files:
        rbp_name = f.split('_')[0]
        df = pd.read_csv(f, sep="\t")
        

        if 'log2FC_Mean' in df.columns:
            val_col = 'log2FC_Mean'
        elif 'log2FC_rep1' in df.columns:
            val_col = 'log2FC_rep1'
        else:
            print(f"No lines, skip.")
            continue

        df_sub = df[['Gene Name', val_col]].rename(columns={val_col: rbp_name})
        

        df_sub = df_sub.groupby('Gene Name', as_index=False).max()
        

        if merged_df is None:
            merged_df = df_sub
        else:
            merged_df = pd.merge(merged_df, df_sub, on="Gene Name", how="outer")
            
    if merged_df is None or merged_df.empty:
        print("Empty.")
        return pd.DataFrame()
        

    merged_df = merged_df.fillna(0)
    merged_df = merged_df.set_index("Gene Name")
    

    rbp_cols = merged_df.columns.tolist()
    total_rbps = len(rbp_cols)
    

    merged_df['Depleted_Count'] = (merged_df[rbp_cols] < 0).sum(axis=1)
    

    filtered_df = merged_df[merged_df['Depleted_Count'] >= min_rbps].copy()
    filtered_df['Mean_Log2FC'] = filtered_df[rbp_cols].mean(axis=1)
    filtered_df = filtered_df.sort_values(by=['Depleted_Count', 'Mean_Log2FC'], ascending=[False, True])
    

    if filtered_df.empty:
        print("No such miRNA.")
    else:
        print(f"Find top 20:\n")
        print(filtered_df.head(20).to_string())
    print("-" * 60)
    
    return filtered_df

if __name__ == "__main__":
    best_negative_controls = find_depleted_mirnas(min_rbps=5)
    
    if not best_negative_controls.empty:
        top_10_mirnas = best_negative_controls.head(10).index.tolist()

        print(f'target_miRNAs = {top_10_mirnas}')