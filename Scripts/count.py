import sys
import pandas as pd

def count_genes(files):
    gene_counts = {}


    for i, file in enumerate(files):
        df = pd.read_csv(file, sep="\t", header=None, low_memory=False)  
        genes_in_file = df[9].tolist() 
        

        for gene in genes_in_file:
            if gene not in gene_counts:
                gene_counts[gene] = [0] * len(files) 
            gene_counts[gene][i] += 1 
    

    result_df = pd.DataFrame.from_dict(gene_counts, orient="index")
    result_df.index.name = 'Gene Name'
    result_df.columns = [f'File{i+1} Count' for i in range(len(files))]
    

    result_df.to_csv(sys.argv[5], sep="\t", header=True)


files = [sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]]

count_genes(files)