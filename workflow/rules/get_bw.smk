## convert the bam file to bigwig format (signal normalized to RPM) for viewing the data on Genome Browse

rule get_bw:
    input: 
        bam = "STAR_align/{sample}.bam",
        rpm_factors = 'report/rpm_factor.txt'
    output: 
        "bw_rpm/{sample}.bw"
    threads: 1
    params:
        mem = '10G',
        # rpmFactor = lambda wildcards: sample_factor_dict[wildcards.sample],
        jobName = "get_bw.{sample}"
    conda: "../envs/rnaseq_env.yaml"
    shell:
        "factor=`cat {input.rpm_factors} |  grep {wildcards.sample} | cut -f2` && "
        "genomeCoverageBed -split -bg -ibam {input.bam} -scale $factor > bw_rpm/{wildcards.sample}.bg && "
        "bedtools sort -i bw_rpm/{wildcards.sample}.bg > bw_rpm/{wildcards.sample}.sort.bg && "
        "bedGraphToBigWig bw_rpm/{wildcards.sample}.sort.bg STAR_index/chrNameLength.txt {output} && "
        "rm bw_rpm/{wildcards.sample}.bg bw_rpm/{wildcards.sample}.sort.bg"

