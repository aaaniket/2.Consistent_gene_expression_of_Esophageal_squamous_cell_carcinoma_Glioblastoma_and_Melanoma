(SAMPLE,)=glob_wildcards("trimmedReads/{sample}.fastq")

rule all:
    input:
        expand("hisat2/{sample}.sam",sample=SAMPLE),

rule hisat2_Alignment:
    input:
        idxdone="hisat2/makeidx.done",
        trim1="trimmedReads/{sample}.fastq",
    output:
        "hisat2/{sample}.sam"
    params:
        idx="hisat2/genome_hisat2.idx",
        threads=20
    shell:
        """
        hisat2 -p {params.threads} -x {params.idx}  -U {input.trim1} -S {output}
        """