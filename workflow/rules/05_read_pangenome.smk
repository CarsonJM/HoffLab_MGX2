# -------------------------------------
# Read-Based Pangenome Module
# -------------------------------------
import pandas as pd


# Load sample information and validate
configfile: "config/config.yaml"


samples_df = pd.read_csv(config["samples_df"], sep="\t")
samples = samples_df["sample"]


# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# load report
report: "../report/workflow.rst"


# -------------------------------------
# Read-Based Pangenome Rules
# -------------------------------------
localrules:
    panphlan_download,
    panphlan_profile,


# -------------------------------------
# 01 Panphlan
# -------------------------------------
# download pangenome for specified species
rule panphlan_download:
    message:
        "Downloading pangenome for {wildcards.species}"
    output:
        resources + "panphlan/{species}/{species}_pangenome.tsv",
    params:
        out_dir=resources + "panphlan/",
    conda:
        "../envs/panphlan:3.1--py_0.yml"
    # container:
    #     "docker://quay.io/biocontainers/panphlan:3.1--py_0"
    benchmark:
        "benchmark/05_READ_BASED_PANGENOME/download_pangenome_{species}.tsv"
    resources:
        runtime="10m",
        mem_mb="10GB",
    shell:
        """
        # download panphlan pangenome
        panphlan_download_pangenome.py \
        --input_name {wildcards.species} \
        --output {params.out_dir} -v
        """


# map samples to pangenomes
rule panphlan_map:
    message:
        "Mapping {wildcards.sample} to {wildcards.species} pangenome"
    input:
        pangenome=resources + "panphlan/{species}/{species}_pangenome.tsv",
        merged=results + "03_READ_BASED_FUNCTION/01_merge_pairs/{sample}.fastq.gz",
    output:
        results + "05_READ_BASED_PANGENOME/01_panphlan/{species}/{sample}.tsv",
    params:
        bowtie2_db=resources + "panphlan/{species}/{species}",
        extra_args=config["read_pangenome"]["panphlan_map_arguments"],
    conda:
        "../envs/panphlan:3.1--py_0.yml"
    # container:
    #     "docker://quay.io/biocontainers/panphlan:3.1--py_0"
    benchmark:
        "benchmark/05_READ_BASED_PANGENOME/panphlan_map_{species}_{sample}.tsv"
    resources:
        runtime=config["read_pangenome"]["panphlan_runtime"],
        mem_mb=config["read_pangenome"]["panphlan_memory"],
    threads: config["read_pangenome"]["panphlan_threads"]
    shell:
        """
        # map to pangenomes
        panphlan_map.py \
        --input {input.merged} \
        --indexes {params.bowtie2_db} \
        --pangenome {input.pangenome} \
        --output {output} \
        {params.extra_args}
        """


# profile mapping results
rule panphlan_profile:
    message:
        "Profiling {wildcards.species} sample mapping"
    input:
        pangenome=resources + "panphlan/{species}/{species}_pangenome.tsv",
        maps=expand(
            results + "05_READ_BASED_PANGENOME/01_panphlan/{{species}}/{sample}.tsv",
            sample=samples,
        ),
    output:
        results + "05_READ_BASED_PANGENOME/01_panphlan/{species}_profile.tsv",
    params:
        in_dir=results + "05_READ_BASED_PANGENOME/01_panphlan/{species}",
        extra_args=config["read_pangenome"]["panphlan_profile_arguments"],
    conda:
        "../envs/panphlan:3.1--py_0.yml"
    # container:
    #     "docker://quay.io/biocontainers/panphlan:3.1--py_0"
    benchmark:
        "benchmark/05_READ_BASED_PANGENOME/panphlan_map_{species}.tsv"
    resources:
        runtime="10m",
        mem_mb="10GB",
    shell:
        """
        # profile based on mapping results
        panphlan_profiling.py \
        --i_dna {params.in_dir} \
        --o_matrix {output} \
        --pangenome {input.pangenome} \
        {params.extra_args}
        """
