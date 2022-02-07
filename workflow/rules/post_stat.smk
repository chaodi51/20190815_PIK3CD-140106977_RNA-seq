rule post_stat:
    input:
        expand("STAR_align/{sample}.Log.final.out", sample=SAMPLES)
    output:
        mappedReads = "report/mapped_reads.txt",
        rpmFactor = "report/rpm_factor.txt"
    params:
        mem = '1G',
        jobName = "post_stat"
    shell:
        '''
        rm -f {output};
        for i in {input}; do
            sampleName="$(basename $i .Log.final.out)";
            cat $i | grep 'Number of input reads' | awk '{{print $6}}' > foo1;
            cat $i | grep 'Uniquely mapped reads number' | awk '{{print $6}}' > foo2;
            cat $i | grep 'Number of reads mapped to multiple loci' | awk '{{print $9}}' > foo3;
            cat $i | grep "Uniquely mapped reads %" | awk '{{print $6}}' > foo4;
            cat $i | grep "% of reads mapped to multiple loci"|awk '{{print $9}}' > foo5;
            paste foo1 foo2 foo3 foo4 foo5 | awk '{{print "'$sampleName'\t"$1"\t"$2+$3"\t"$4+$5}}' >> {output.mappedReads}
        done
        cat {output.mappedReads} | awk '{{print $1"\t"1000000/$2}}' > {output.rpmFactor}
        sed -i '1isample\ttotal_reads\tmapped_reads\t%mapped' {output.mappedReads};
        rm -f foo*
        '''
