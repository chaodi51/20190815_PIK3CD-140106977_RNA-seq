
## download fastq files from GEO
# 1, get SRA numbers from SRA Run Selector and save in "sra.id"
# 2, get corresponding sample ids with GSMXXX on the GSEXXXX webpage, save in "gsm.id", make sure the order is correct
# 3, put 1 and 2 to a single file called "sample_table.tsv"
rule download_fastq:
    output:
        "raw_fq/{sample}_1.fastq", 
        "raw_fq/{sample}_2.fastq"
    params:
        outdir = raw_fq,
        mem = '6G',
        jobName = "download_fastq.{sample}",
        sraID = lambda wildcards: sample_sra_pair[wildcards.sample]
    conda: "../envs/rnaseq_env.yaml"
    shell:
        "fastq-dump -A {wildcards.sample} --gzip -O {params.outdir} --split-files {params.sraID}"
