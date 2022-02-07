## from https://github.com/chaodi51/snakemake-workflows-rna-seq-star-deseq2/blob/master/scripts/count-matrix.py
import pandas as pd

def get_count_column(strand):
    if strand == "no":
        return 1 #non stranded protocol
    elif strand == "yes":
        return 2 #3rd column in STAR output {sample}.ReadsPerGene.out.tab
    elif strand == "reverse":
        return 3 #4th column, usually for Illumina truseq
    else:
        raise ValueError("'strand' column should be empty or has the value of 'none', 'yes' or 'reverse'!!!")

counts = [pd.read_table(count_file, index_col=0, usecols=[0, get_count_column(strand)], 
          header=None, skiprows=4) 
          for count_file, strand in zip(snakemake.input, snakemake.params.strand_list)]

for df, sample in zip(counts, snakemake.params.sample_list):
    df.columns = [sample]

matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
# collapse technical replicates
# matrix = matrix.groupby(matrix.columns, axis=1).sum()
matrix.to_csv(snakemake.output[0], sep="\t")
