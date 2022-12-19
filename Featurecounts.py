
(SAMPLE,) = glob_wildcards("aligned/{sample}.sorted.bam")

rule all:
    input:
    
        expand("aligned/{sample}.sorted.bam.bai", sample=SAMPLE),
        expand("samtools_stats/{sample}.stats.txt", sample=SAMPLE),
        expand("samtools_stats/{sample}.flagstat.txt", sample=SAMPLE),
        #rawCounts
        "raw_Counts",
    



rule samtools_index:
    input:
        "aligned/{sample}.sorted.bam",
    output:
        "aligned/{sample}.sorted.bam.bai",
    shell:
        "samtools index {input} {output} "

rule samtools_stats:
    input:
        bam="aligned/{sample}.sorted.bam",
    output:
        "samtools_stats/{sample}.stats.txt",
    shell:
       "samtools stats {input.bam} > {output} "

rule samtools_flagstat:
    input:
        "aligned/{sample}.sorted.bam",
    output:
        "samtools_stats/{sample}.flagstat.txt",
    shell:
        "samtools flagstat {input} > {output} "
        

rule featureCounts:
    input:
        samples=expand("aligned/{sample}.bam", sample=SAMPLE), 
        annotation = rules.get_genome_gtf.output.gtf
    output:
        "raw_Counts"
    threads:
        16
    shell:
        "featureCounts -T {threads} -a {input.annotation} -o {output} {input.samples}"
