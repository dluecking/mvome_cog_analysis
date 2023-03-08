
configfile: "config.yml"



################################################################################
# all
rule all:
    input:
        expand("reads/{sample}_true.mvome.reads1.subsample.fq.gz", sample = config["samples"]),
        expand("reads/{sample}_true.mvome.reads2.subsample.fq.gz", sample = config["samples"])

################################################################################
# single rules
rule subsample_reads:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        r1="reads/{sample}_true.mvome.reads1.fq.gz",
        r2="reads/{sample}_true.mvome.reads2.fq.gz"
    output:
        r1="reads/{sample}_true.mvome.reads1.subsample.fq",
        r2="reads/{sample}_true.mvome.reads2.subsample.fq"
    threads: 8
    params: SUBSAMPLE_SIZE=config["SUBSAMPLE_SIZE"]
    log: "log/subsample_reads.{sample}.log"
    shell:
        """
        seqtk sample -s seed=100 {input.r1} {params.SUBSAMPLE_SIZE} \
        > {output.r1}
        
        seqtk sample -s seed=100 {input.r2} {params.SUBSAMPLE_SIZE} \
        > {output.r2}
        """

rule zip_subsamples:
    input:
        r1="reads/{sample}_true.mvome.reads1.subsample.fq",
        r2="reads/{sample}_true.mvome.reads2.subsample.fq"
    output:
        r1="reads/{sample}_true.mvome.reads1.subsample.fq.gz",
        r2="reads/{sample}_true.mvome.reads2.subsample.fq.gz"
    threads: 8
    shell:
        """
        gzip {input.r1}
        """


rule concatenate_MAGs_by_type:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        top20_df=config["top20_df"]
    output:
        "MAGs_sorted/sorting_done.log"
    threads: 8
    log: "log/concatenate_MAGs_by_type.log"
    shell:
        "Rscript scripts/concatenate_MAGs_by_type.R && touch {output}"
    



