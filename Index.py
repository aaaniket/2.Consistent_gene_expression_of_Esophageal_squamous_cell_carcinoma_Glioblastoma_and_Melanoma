
rule hisat2_indexing:
    input:
        ref="genome.fa"
    output:
        touch("hisat2/makeidx.done")
    params:
        threads=20,
        idx="hisat2/genome_hisat2.idx"
    shell:
        """
        hisat2-build -p {params.threads} {input.ref} {params.idx}
        """
            