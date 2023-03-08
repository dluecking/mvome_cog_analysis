"""
command to run:
snakemake -pr -j 16 --use-conda --rerun-incomplete --conda-prefix ~/envs/ \
--cluster  "sbatch --time={resources.time} --mem={resources.mem} --cpus={resources.cpus}"
"""
 
configfile: "config.yml"
localrules: subsample_reads, zip_subsamples, concatenate_MAGs_by_type

################################################################################
# all
rule all:
    input:
        "plots/cog_barplot.png",
        "plots/cog_barplot.svg"

################################################################################
# single rules

rule plotting:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        "intermediate/collected_data.tsv"
    output:
        png="plots/cog_barplot.png",
        svg="plots/cog_barplot.svg"
    shell:
        """
        Rscript scripts/plotting.R
        """
rule collect_blast_out_data:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        "intermediate/donefiles/blast_microbial.done",
        "intermediate/donefiles/blast_peDNA.done"
    output:
        "intermediate/collected_data.tsv" 
    shell:
        """
        Rscript scripts/collect_data.R
        """

rule blast_microbial:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        "============================================="
    output:
        "intermediate/donefiles/blast_microbial.done"
    shell:
        """
        echo alala
        """

rule blast_peDNA:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        "============================================="
    output:
        "intermediate/donefiles/blast_peDNA.done"
    shell:
        """
        echo lalal
        """




rule subsample_reads:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        r1="input/reads/{sample}_true.mvome.reads1.fq.gz",
        r2="input/reads/{sample}_true.mvome.reads2.fq.gz"
    output:
        r1="intermediate/reads/{sample}_true.mvome.reads1.subsample.fq",
        r2="intermediate/reads/{sample}_true.mvome.reads2.subsample.fq"
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
        r1="intermediate/reads/{sample}_true.mvome.reads1.subsample.fq",
        r2="intermediate/reads/{sample}_true.mvome.reads2.subsample.fq"
    output:
        r1="intermediate/reads/{sample}_true.mvome.reads1.subsample.fq.gz",
        r2="intermediate/reads/{sample}_true.mvome.reads2.subsample.fq.gz"
    threads: 8
    shell:
        """
        gzip {input.r1}
        gzip {input.r2}
        """


rule concatenate_MAGs_by_type:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        top20_df=config["top20_df"]
    output:
        "intermediate/donefiles/MAG_sorting_done.log"
    threads: 2
    log: "log/concatenate_MAGs_by_type.log"
    shell:
        "Rscript scripts/concatenate_MAGs_by_type.R {input} && touch {output}"
    

rule map_subsampled_reads_to_sorted_MAGs:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        r1="intermediate/reads/{sample}_true.mvome.reads1.subsample.fq.gz",
        r2="intermediate/reads/{sample}_true.mvome.reads2.subsample.fq.gz",
        MAG_sorting_done="intermediate/donefiles/sorting_done.log"
    output:
        outfile="intermediate/donefiles/{sample}_mapping_done.log"
    params: minid=config["MIN_MAPPING_ID"]
    resources:
        mem=config["sbatch_mem"],
        cpus=config["sbatch_cpus"],
        time=config["sbatch_time"]
    shell:
        """
        for label in {ev_producer,gta,viral};
        do
            bbmap.sh \
            ref=intermediate/${station}_combined_${label}_contigs.fa \
            in={input.r1} in2={input.r2} \
            outm=intermediate/mapped_reads/${station}_mapped_to_${label}.reads1.fasta.gz \
            outm2=intermediate/mapped_reads/${station}_mapped_to_${label}.reads2.fasta.gz  \
            minidentity={params.minid} \
            nodisk
        done
        touch {output.outfile}
        """


