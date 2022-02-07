## return a list of library strand, as 'none', 'yes' or 'reverse' denoted as in HTSeq -s option
def get_strandness(sample_table):
    if "strand" in sample_table.columns:
        return sample_table["strand"].tolist()

def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]

# def get_deseq2_threads(wildcards=None):
#     # https://twitter.com/mikelove/status/918770188568363008
#     few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
#     return 1 if len(sample_list) < 100 or few_coeffs else 6   

rule count_matrix:
    input:
        expand("STAR_align/{sample}.ReadsPerGene.out.tab", sample=SAMPLES)
    output:
        "../results/tables/all_readCount.tsv"
    params:
        sample_list = sample_table["sample"].tolist(),
        strand_list = get_strandness(sample_table),
        mem = '2G',
        jobName = "count_matrix" 
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"

rule deseq2_allsample_plot:
    input:
        count_table="../results/tables/all_readCount.tsv"
    output:
        heatmap_exp = report("../results/diffexp/heatmap_expression_matrix.pdf", "../report/heatmap.rst"),
        heatmap_dist = report("../results/diffexp/heatmap_sample_distance.pdf", "../report/heatmap_dist.rst"),
        pca_plot = report("../results/diffexp/pca-plot.pdf", "../report/pca.rst")
    params:
        sample_table = config["samples"],
        mem = '2G',
        jobName = 'deseq2_allsample_plot'
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/deseq2_allsample_plot.log"
    threads: 1
    script:
        "../scripts/deseq2_allsample_plot.R"
    

rule deseq2_diffexp:
    input:
        count_table="../results/tables/all_readCount.tsv"
    output:
        table = "../results/diffexp/{contrast}.diffexp.tsv",
        ma_plot = "../results/diffexp/{contrast}.ma-plot.pdf"
    params:
        sample_table = config["samples"],
        mem = '2G',
        jobName = 'deseq2_diffexp.{contrast}'
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    threads: 1
    script:
        "../scripts/deseq2_diffexp.R"

