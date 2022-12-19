(SAMPLE,) = glob_wildcards("aligned/{sample}.sorted.bam")

rule all:
    input: 
        expand ("aligned/{sample}_hisat2_sorted.bam.bai"),
    

rule samtools_index:
    input:
        "aligned/{sample}_hisat2_sorted.bam",       
    output:
        "aligned/{sample}sorted.bam.bai",
    params:
        threads=20,
        idx="hisat2/genome_hisat2.idx"
    shell:
        "samtools index {input} {output}"
