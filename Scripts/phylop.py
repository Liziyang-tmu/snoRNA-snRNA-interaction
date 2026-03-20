import pyBigWig
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

def extract_phylop_scores(bed_file, bw_file):

    bw = pyBigWig.open(bw_file)
    
    data = []
    
    with open(bed_file, 'r') as f:
        for line in f:
            if line.strip() == "": continue
            parts = line.strip().split()
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            
 
            length = end - start
            
 
            up_start, up_end = start - length, start
            down_start, down_end = end, end + length
            
            try:

                score_interaction = bw.stats(chrom, start, end, type="mean")[0]
                score_up = bw.stats(chrom, up_start, up_end, type="mean")[0]
                score_down = bw.stats(chrom, down_start, down_end, type="mean")[0]
                

                if all(v is not None for v in [score_interaction, score_up, score_down]):

                    score_flank = (score_up + score_down) / 2.0
                    
                    data.append({
                        'Region': 'Flanking',
                        'PhyloP_Score': score_flank,
                        'Pair_ID': f"{chrom}:{start}-{end}"
                    })
                    data.append({
                        'Region': 'Interaction',
                        'PhyloP_Score': score_interaction,
                        'Pair_ID': f"{chrom}:{start}-{end}"
                    })
            except RuntimeError:

                continue
                
    bw.close()
    return pd.DataFrame(data)


bed_file = "interaction.bed"      
bw_file = "hg38.phyloP100way.bw"        


df = extract_phylop_scores(bed_file, bw_file)


interaction_scores = df[df['Region'] == 'Interaction']['PhyloP_Score'].values
flank_scores = df[df['Region'] == 'Flanking']['PhyloP_Score'].values

stat, p_val = wilcoxon(interaction_scores, flank_scores, alternative='greater')
print(f"Paired Wilcoxon P-value (Interaction > Flanking): {p_val}")


sns.set_theme(style="ticks", font_scale=1.2)
plt.figure(figsize=(4.5, 6))


ax = sns.boxplot(x='Region', y='PhyloP_Score', data=df, 
                 palette=['#CCCCCC', '#E64B35'], width=0.5, showfliers=False)


for pair_id in df['Pair_ID'].unique():
    pair_data = df[df['Pair_ID'] == pair_id]
    plt.plot([0, 1], 
             [pair_data[pair_data['Region']=='Flanking']['PhyloP_Score'].values[0],
              pair_data[pair_data['Region']=='Interaction']['PhyloP_Score'].values[0]], 
             color='gray', alpha=0.3, linewidth=1, zorder=1)


sns.stripplot(x='Region', y='PhyloP_Score', data=df, 
              color='black', alpha=0.6, size=5, zorder=2, ax=ax)


y_max = df['PhyloP_Score'].max()
y_range = y_max - df['PhyloP_Score'].min()
y_annot = y_max + y_range * 0.05

plt.plot([0, 0, 1, 1], [y_annot, y_annot+0.05, y_annot+0.05, y_annot], lw=1.5, c='k')

if p_val < 0.001:
    p_text = "***"
elif 0.001< p_val < 0.01:
    p_text = "**"
elif 0.01 < p_val < 0.05:
    p_text = "*"
else:
    p_text = "ns"
plt.text(0.5, y_annot+0.06, p_text, ha='center', va='bottom', color='k', fontsize=14)


plt.ylabel("PhyloP Score (100-way vertebrates)")
plt.xlabel("")
plt.title("Conservation of Interaction Sites", pad=15)
plt.ylim(df['PhyloP_Score'].min() - 0.5, y_annot + y_range * 0.15)
sns.despine()

plt.tight_layout()
plt.savefig('PhyloP_Conservation_Paired.pdf', bbox_inches='tight',dpi=600)
plt.show()