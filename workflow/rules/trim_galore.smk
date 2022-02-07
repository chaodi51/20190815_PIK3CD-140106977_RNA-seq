rule trim_galore:
    input: 
        get_fastq ## call the function in main Snakefile
    output:
        fq1 = "trimmed_fq/{sample}_1_trimmed.fq.gz",
        fq1_rep = "trimmed_fq/{sample}_1.fq.gz_trimming_report.txt",
        fq2 = "trimmed_fq/{sample}_2_trimmed.fq.gz",
        fq2_rep = "trimmed_fq/{sample}_2.fq.gz_trimming_report.txt"
    log: "logs/trimmed/{sample}.log"
    threads: 4
    params:
        outdir = trimmed_fq,
        mem = '6G',
        jobName = "trim_galore.{sample}" 
    conda: "../envs/rnaseq_env.yaml"
    shell:
        "trim_galore --cores {threads} --gzip --fastqc --paired -o {params.outdir} {input} &> {log} && "
        "mv trimmed_fq/{wildcards.sample}_1_val_1.fq.gz {output.fq1} &>> {log} && "
        "mv trimmed_fq/{wildcards.sample}_2_val_2.fq.gz {output.fq2} &>> {log}"  