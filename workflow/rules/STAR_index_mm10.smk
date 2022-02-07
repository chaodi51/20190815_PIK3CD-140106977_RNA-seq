shell.prefix("source ~/.bash_profile; ")

import os

configfile: "../config/config_star_index_mm10.yaml"

# snakemake -s rules/STAR_index.smk -c "qsub -l h_vmem={params.mem} -l mem_free={params.mem} -pe smp {threads} -V -cwd -e qsub/{params.jobName}.e -o qsub/{params.jobName}.o" -j -p

## STAR indexing the genome
rule star_index:
    input:
        fa = config["genome"]["sequence"], # reference FASTA file
        gtf = config["genome"]["annotation"] # GTF file
    output:
        directory("/mnt/isilon/dbhi_bfx/dic/genomes/mm10/STAR_mm10_index") # the index folder
    threads: 6 # set the maximum number of available cores
    params: 
        mem = "10G",
        jobName = "star_index"
    shell:
        """
        conda activate snakemake
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} \\
        --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang 100
        """

