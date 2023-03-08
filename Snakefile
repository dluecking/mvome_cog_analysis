"""
command to run:
snakemake -pr -j 16 --use-conda --conda-prefix ~/envs/ --rerun-incomplete \
--cluster  "sbatch --time={resources.time} --mem={resources.mem} \
--cpus-per-task={resources.cpus} --output=log/slurm/%x.%j.out --partition=CLUSTER" \
--latency-wait 60
"""
 
configfile: "config.yml"
localrules: all, plotting, collect_blast_out_data, unzip_mapped_reads, subsample_reads, zip_peDNA_subsamples, convert_fastq_to_fasta

# custom expand
def customExpand():
    import pandas as pd 
    data = pd.read_table("input/top20_df.tsv")
    data = data[(data['label_after_cov'] != "unclear")]
    data = data[(data['label_after_cov'] != "remove")]
    data = pd.DataFrame(data.groupby(['Sample', 'label_after_cov']))[0]

    list_to_return = []

    for i in range(0, len(data)):
        sample = data[i][0]
        label = data[i][1]
        r1=f'intermediate/blast_out/peDNA/{sample}_mapped_to_{label}_r1.out'
        r2=f'intermediate/blast_out/peDNA/{sample}_mapped_to_{label}_r2.out'
        list_to_return.append(r1)
        list_to_return.append(r2)
        
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
        expand("intermediate/blast_out/microbial/{sample}_r{read}_blast.out", sample = config["samples"], read = [1,2]),
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
        "input/microbial_reads/{sample}_microbial.reads{read}.1M.fa"
    output:
        "intermediate/protein_fragments/microbial/{sample}_microbial_fraggene_r{read}.faa.faa",
    resources:
        mem=config["sbatch_mem"],
        cpus=config["sbatch_cpus"],
        time=config["sbatch_time"]
    shell:
        """
        run_FragGeneScan.pl -genome={input} -complete=0 -train=illumina_5 \
        -out {output} -thread={resources.cpus}
        """

rule blast_microbial:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        "intermediate/protein_fragments/microbial/{sample}_microbial_fraggene_r{read}.faa.faa"
    output:
        "intermediate/blast_out/microbial/{sample}_r{read}_blast.out"
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
        "intermediate/protein_fragments/peDNA/{sample}_mapped_to_{label}_r{read_file}.faa.faa"
    output:
        "intermediate/blast_out/peDNA/{sample}_mapped_to_{label}_r{read_file}.out"
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
        "intermediate/mapped_reads/{sample}_mapped_to_{label}_r{reads}.fa"
    output:
        "intermediate/protein_fragments/peDNA/{sample}_mapped_to_{label}_r{reads}.faa.faa"
    resources:
        mem=config["sbatch_mem"],
        cpus=config["sbatch_cpus"],
        time=config["sbatch_time"]
    shell:
        """
        run_FragGeneScan.pl -genome={input} -complete=0 -train=illumina_5 \
        -out {output} -thread {resources.cpus}
        """

rule convert_fastq_to_fasta:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        "intermediate/mapped_reads/{sample}_mapped_to_{label}_r{reads}.fq",
    output:
        "intermediate/mapped_reads/{sample}_mapped_to_{label}_r{reads}.fa",
    threads: 8
    shell:
        """
        seqtk seq -a {input} > {output}
        """


rule unzip_mapped_reads:
    input:
        "intermediate/mapped_reads/{sample}_mapped_to_{label}_r{reads}.fq.gz",
    output:
        "intermediate/mapped_reads/{sample}_mapped_to_{label}_r{reads}.fq",
    threads: 8
    shell:
        """
        gunzip {input}
        """

rule map_subsampled_reads_to_sorted_MAGs:
    conda:
        "env/mvome_cog_analysis.yml"
    input:
        r1="intermediate/reads/{sample}_true.mvome.reads1.subsample.fq.gz",
        r2="intermediate/reads/{sample}_true.mvome.reads2.subsample.fq.gz",
        reference="input/MAGs_sorted/{sample}_combined_{label}_contigs.fa"
    output:
        r1="intermediate/mapped_reads/{sample}_mapped_to_{label}_r1.fq.gz",
        r2="intermediate/mapped_reads/{sample}_mapped_to_{label}_r2.fq.gz",
    params: minid=config["MIN_MAPPING_ID"]
    resources:
        mem=config["sbatch_mem"],
        cpus=config["sbatch_cpus"],
        time=config["sbatch_time"]
    shell:
        """
        bbmap.sh \
        ref={input.reference} \
        in={input.r1} in2={input.r2} \
        outm={output.r1} outm2={output.r2}  \
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