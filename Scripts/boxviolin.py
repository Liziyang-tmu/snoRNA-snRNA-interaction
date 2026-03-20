import pandas as pd
import io
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon


df = pd.read_csv('NMFE.csv') 


stat_wilcox, p_val_wilcox = wilcoxon(df['Interactions'], df['Random'], alternative='less')


df_melted = df.melt(id_vars='DGs', 
                    value_vars=['Interactions', 'Random'], 
                    var_name='Group', 
                    value_name='Normalized_MFE')


sns.set_theme(style="ticks", font_scale=1.2) 
plt.figure(figsize=(5, 6))
my_palette = ['#E64B35', '#4DBBD5']




ax = sns.violinplot(
    x='Group', y='Normalized_MFE', data=df_melted, 
    palette=my_palette, 
    inner=None,        
    width=0.6,           
    linewidth=1.5,
    alpha=0.6          
)


sns.boxplot(
    x='Group', y='Normalized_MFE', data=df_melted,
    palette=my_palette,  
    width=0.2,          
    boxprops=dict(edgecolor='black', linewidth=1.5), 
    whiskerprops=dict(color='black', linewidth=1.5),
    capprops=dict(color='black', linewidth=1.5),
    medianprops=dict(color='black', linewidth=2),
    showfliers=False,    
    ax=ax                
)


sns.stripplot(
    x='Group', y='Normalized_MFE', data=df_melted,
    color='black',       
    alpha=0.5,           
    jitter=True,         
    size=5,              
    zorder=3,            
    ax=ax
)
# --------------------------------------------


y_max = df_melted['Normalized_MFE'].max()
y_range = df_melted['Normalized_MFE'].max() - df_melted['Normalized_MFE'].min()
y_annot = y_max + y_range * 0.08 

plt.plot([0, 0, 1, 1], [y_annot, y_annot+0.02, y_annot+0.02, y_annot], lw=1.5, c='k')


if p_val_wilcox < 0.0001:
    p_text = "****"
elif 0.001 < p_val_wilcox < 0.01:
    p_text = f"**"
elif 0.01 < p_val_wilcox < 0.05:
    p_text = f"*"
elif 0.0001 < p_val_wilcox < 0.001:
    p_text = f"***"
else:
    p_text = f"ns"
    
plt.text(0.5, y_annot+0.03, p_text, ha='center', va='bottom', color='k', fontsize=12)


plt.ylabel("Normalized Minimum Free Energy (NMFE)")
plt.xlabel("") 
plt.title("Random vs Verified Interactions", pad=15)
plt.ylim(df_melted['Normalized_MFE'].min() - 0.1, y_annot + 0.15) 
sns.despine() 

plt.tight_layout()
plt.savefig('NMFE_Violin_Plot.pdf', bbox_inches='tight', dpi=600) 
plt.show()

print(f"P-value: {p_val_wilcox}")