
(SAMPLE,) = glob_wildcards("aligned/{sample}_hisat2_sorted.bam")

rule all:
    input:
        expand("aligned/{sample}_hisat2_sorted.bam.bai", sample=SAMPLE),
        expand("samtools_stats/{sample}.stats.txt", sample=SAMPLE),
        expand("samtools_stats/{sample}.flagstat.txt", sample=SAMPLE),
        
    
rule samtools_index:
    input:
        "aligned/{sample}_hisat2_sorted.bam",
        "hisat2/genome_hisat2.idx"
    output:
        "aligned/{sample}_hisat2_sorted.bam.bai",
        "hisat2/genome_hisat2.idx"
    shell:
        "samtools idx {input} {output}"

rule samtools_stats:
    input:
        bam="aligned/{sample}_hisat2_sorted.bam",
    output:
        "samtools_stats/{sample}.stats.txt",
    shell:
       "samtools stats {input.bam} > {output} "

rule samtools_flagstat:
    input:
        "aligned/{sample}_hisat2_sorted.bam",
    output:
        "samtools_stats/{sample}.flagstat.txt",
    shell:
        "samtools flagstat {input} > {output} "
