
rule all:
    input:
        #get fa and gtf files
        "genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
        "genome/Homo_sapiens.GRCh38.106.gtf.gz",
        #Get annotation GTF

rule get_genome_gtf:
    "Downloading Genome annotation file from Ensemble, Homo sapiens primary assembly (GRCh38)"
    output:
        gtf = "genome/Homo_sapiens.GRCh38.106.gtf.gz"
    shell:
        "cd genome"
        " && wget ftp://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz"
        " && gunzip -k Homo_sapiens.GRCh38.106.gtf.gz 