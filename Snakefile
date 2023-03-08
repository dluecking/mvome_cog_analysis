"""
command to run:
snakemake -pr -j 16 --use-conda --rerun-incomplete --conda-prefix ~/envs/ \
--cluster  "sbatch --time={resources.time} --mem={resources.mem} --cpus={resources.cpus}"
"""
 
configfile: "config.yml"
localrules: subsample_reads, zip_subsamples, concatenate_MAGs_by_type



# custom expand
def customExpand():
    import pandas as pd 
    data = pd.read_csv("input/top20_df.tsv")
    data = data[(data['label_after_cov'] != "unclear")]
    data = data[(data['label_after_cov'] != "remove")]
    data = pd.DataFrame(data.groupby(['Sample', 'label_after_cov']))[0]

    list_to_return = []

    for i in range(0, len(d)):
        sample = f[i][0]
        label = f[i][1]
        list_to_return.append(f'file/{sample}_vs_{label}.fa')
    return(list_to_return)


################################################################################
# all
rule all:
    input:
        "plots/cog_barplot.png",
        "plots/cog_barplot.svg"

################################################################################
# rules for all other rules
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
        expand("intermediate/blast_out/microbial/{sample}_r{read}.out", sample = config["samples"], read = [1,2]),
        customExpand()
    output:
        "intermediate/collected_data.tsv" 
    shell:
        """
        Rscript scripts/collect_data.R
        """


# microbial path ###############################################################

rule fraggenescane_microbial:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        r1="input/microbial_reads/{sample}_microbial.reads{read}.1M.fq"
        r2="input/microbial_reads/{sample}_microbial.reads{read}.1M.fq"
    output:
        fg1="intermediate/protein_fragments/microbial/{sample}_microbial_fraggene_r{read}.faa",
        fg2="intermediate/protein_fragments/microbial/{sample}_microbial_fraggene_r{read}.faa"
    resources:
        mem=config["sbatch_mem"],
        cpus=config["sbatch_cpus"],
        time=config["sbatch_time"]
    shell:
        """
        run_FragGeneScan.pl -genome={input.r1} -complete=0 -train=illumina_5 \
        -out {output.fg1} -thread {resources.cpus}

        run_FragGeneScan.pl -genome={input.r2} -complete=0 -train=illumina_5 \
        -out {output.fg2} -thread {resources.cpus}
        """

rule blast_microbial:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        "intermediate/protein_fragments/microbial/{sample}_microbial_fraggene_r{read}.faa"
    output:
        "intermediate/blast_out/microbial/{sample}_microbial_fraggene_r{read}_blast.out"
    params: COG_DB=config["COG_DB_PATH"]
    shell:
        """
        diamond blastp --query {input} \
        --db {params.COG_DB} \
        --out {output} \
        -f 6 \
        --max-target-seqs 1 \
        --query-cover 80 \
        --subject-cover 10
        """


# peDNA path ################################################################### 

rule blast_peDNA:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        "intermediate/protein_fragments/peDNA/{sample}_mapped_to_{label}_r{read_file}.faa"
    output:
        "intermediate/blast_out/peDNA/{sample}_mapped_to_{label}_r{read_file}.faa"
    params: COG_DB=config["COG_DB_PATH"]
    shell:
        """
        diamond blastp --query {input} \
        --db {params.COG_DB} \
        --out {output} \
        -f 6 \
        --max-target-seqs 1 \
        --query-cover 80 \
        --subject-cover 10
        """

rule fraggenescane_peDNA:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        "intermediate/reads/{sample}_mapped_to_{label}_r{read_file}.fq"
    output:
        "intermediate/protein_fragments/peDNA/{sample}_mapped_to_{label}_r{read_file}.faa"
    resources:
        mem=config["sbatch_mem"],
        cpus=config["sbatch_cpus"],
        time=config["sbatch_time"]
    shell:
        """
        run_FragGeneScan.pl -genome={input} -complete=0 -train=illumina_5 \
        -out {output} -thread {resources.cpus}
        """

rule unzip_mapped_reads:
    input:
        r1="intermediate/reads/{sample}_mapped_to_{label}_r1.fq.gz",
        r2="intermediate/reads/{sample}_mapped_to_{label}_r2.fq.gz"
    output:
        r1="intermediate/reads/{sample}_mapped_to_{label}_r1.fq",
        r2="intermediate/reads/{sample}_mapped_to_{label}_r2.fq"
    threads: 8
    shell:
        """
        gzip {input.r1}
        gzip {input.r2}
        """

rule map_subsampled_reads_to_sorted_MAGs:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        r1="intermediate/reads/{sample}_true.mvome.reads1.subsample.fq.gz",
        r2="intermediate/reads/{sample}_true.mvome.reads2.subsample.fq.gz",
        MAG_sorting_done="intermediate/donefiles/sorting_done.log"
    output:
        r1="intermediate/reads/{sample}_mapped_to_{label}_r1.fq.gz",
        r2="intermediate/reads/{sample}_mapped_to_{label}_r2.fq.gz"
    params: minid=config["MIN_MAPPING_ID"]
    resources:
        mem=config["sbatch_mem"],
        cpus=config["sbatch_cpus"],
        time=config["sbatch_time"]
    shell:
        """
        bbmap.sh \
        ref=intermediate/${station}_combined_${label}_contigs.fa \
        in={input.r1} in2={input.r2} \
        outm={output.r1} \
        outm2={output.r2}  \
        minidentity={params.minid} \
        nodisk
        """


# prep work ####################################################################
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

rule zip_peDNA_subsamples:
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
    



