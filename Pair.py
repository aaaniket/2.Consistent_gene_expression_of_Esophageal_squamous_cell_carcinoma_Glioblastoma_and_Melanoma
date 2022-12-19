
r1="Reads/{sample}_1.fastq.gz"
r2="Reads/{sample}_2.fastq.gz"

(SAMPLES,)=glob_wildcards("Reads/{sample}_1.fastq.gz")


rule all:
    input:
        expand("hisat2/{sample}.sam",sample=SAMPLES),
       

rule hisat2_Alignment:
    input:
        idxdone="hisat2/makeidx.done",
        trim1="Reads/{sample}_1.fastq.gz",
        trim2="Reads/{sample}_2.fastq.gz",
    output:
        "hisat2/{sample}.sam"
    params:
        idx="hisat2/genome_hisat2.idx",
        threads=5
    shell:
        """
        hisat2 -p {params.threads} -x {params.idx}  -1 {input.trim1} -2 {input.trim2} -S {output}
        """
        
