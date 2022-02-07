# run STAR mapping with only uniquely mapped reads
rule star_map:
    input:
        fq1 = "trimmed_fq/{sample}_1_trimmed.fq.gz",
        fq2 = "trimmed_fq/{sample}_2_trimmed.fq.gz"
    output:
        "STAR_align/{sample}.bam",
        "STAR_align/{sample}.Log.final.out"
    params:
        genome_dir = STAR_index,
        annotation = gtf_file,
        mem = '10G',
        jobName = "star_map.{sample}" 
    threads: 4
    conda: "../envs/rnaseq_env.yaml"
    shell:
        "STAR --runThreadN {threads} --genomeDir {params.genome_dir} "
        "--outFileNamePrefix STAR_align/{wildcards.sample}. --outSAMtype BAM SortedByCoordinate "
        "--outSAMmapqUnique 255 --outFilterMultimapNmax 1 "
        "--quantMode GeneCounts --sjdbGTFfile {params.annotation} "
        "--readFilesIn {input} --outSJfilterReads Unique --readFilesCommand gunzip -c && "
        "mv STAR_align/{wildcards.sample}.Aligned.sortedByCoord.out.bam STAR_align/{wildcards.sample}.bam"

rule samtools_index:
    input:
        "STAR_align/{sample}.bam"
    output:
        "STAR_align/{sample}.bam.bai"
    params: 
      mem = '10G',
      jobName = "samtools_index.{sample}"
    conda: "../envs/rnaseq_env.yaml"
    shell:
        "samtools index {input} {output}"