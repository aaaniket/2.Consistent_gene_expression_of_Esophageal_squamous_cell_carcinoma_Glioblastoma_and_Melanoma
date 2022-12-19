(SAMPLE,)=glob_wildcards("hisat2/{sample}.sam")

rule all:
    input:
        expand("hisat2/{sample}_hisat2_sorted.bam",sample=SAMPLE),
        
                
rule convert_samtobam:
    input:
        "hisat2/{sample}.sam"
    output:
        "hisat2/{sample}_hisat2_sorted.bam"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} -b -o {output} {input}
        """
        
