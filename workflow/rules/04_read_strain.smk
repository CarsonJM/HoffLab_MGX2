# -------------------------------------
# Read-Based Phylogeny Module
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
# Read-Based Phylogeny Rules
# -------------------------------------
# -----------------------------------------------------
# 01 StrainPhlan
# -----------------------------------------------------
# reconstruct all species strains in each sample
rule sample2markers:
    input:
        results + "02_READ_BASED_TAXONOMY/01_metaphlan/{sample}.sam.bz2",
        spa_db=resources + "metaphlan/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl",
    output:
        results + "04_READ_BASED_STRAIN/01_consensus_markers/{sample}.pkl",
    params:
        out_dir=results + "04_READ_BASED_STRAIN/01_consensus_markers/",
    # conda:
    #     "../envs/humann:3.6--pyh7cba7a3_1.yml"
    container:
        "docker://quay.io/biocontainers/humann:3.6--pyh7cba7a3_1"
    benchmark:
        "benchmark/04_READ_BASED_STRAIN/sample2markers_{sample}.tsv"
    resources:
        runtime=config["read_strain"]["strainphlan_runtime"],
        mem_mb=config["read_strain"]["strainphlan_memory"],
    threads: config["read_strain"]["strainphlan_threads"]
    shell:
        """
        # combine paired end reads
        sample2markers.py -i {input} \
        -o {params.out_dir} \
        --nprocs {threads} \
        --input_format bz2 \
        --database {input.spa_db}
        """


# identify all clades for strainphlan
rule strainphlan_clades:
    input:
        expand(
            results + "04_READ_BASED_STRAIN/01_consensus_markers/{sample}.pkl",
            sample=samples,
        ),
        spa_db=resources + "metaphlan/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl",
    output:
        results + "04_READ_BASED_STRAIN/02_clades/strainphlan_clades",
    params:
        out_dir=results + "04_READ_BASED_STRAIN/02_clades/",
    # conda:
    #     "../envs/humann:3.6--pyh7cba7a3_1.yml"
    container:
        "docker://quay.io/biocontainers/humann:3.6--pyh7cba7a3_1"
    benchmark:
        "benchmark/04_READ_BASED_STRAIN/strainphlan_clades.tsv"
    resources:
        runtime=config["read_strain"]["strainphlan_runtime"],
        mem_mb=config["read_strain"]["strainphlan_memory"],
    threads: config["read_strain"]["strainphlan_threads"]
    shell:
        """
        # combine paired end reads
        strainphlan \
        --samples {input} \
        --database {input.spa_db} \
        --nprocs {threads} \
        --output_dir {params.out_dir} \
        --print_clades_only > {output}
        """


# extract markers
rule extract_markers:
    message:
        "Extracting strainphlan markers for {wildcards.species}"
    input:
        resources + "metaphlan/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl",
    output:
        resources + "strainphlan/markers/{species}.fna",
    params:
        out_dir=resources + "strainphlan/markers/",
    # conda:
    #     "../envs/humann:3.6--pyh7cba7a3_1.yml"
    container:
        "docker://quay.io/biocontainers/humann:3.6--pyh7cba7a3_1"
    benchmark:
        "benchmark/04_READ_BASED_STRAIN/extract_markers_{species}.tsv"
    resources:
        runtime="1h",
        mem_mb="10GB",
    threads: config["read_strain"]["strainphlan_threads"]
    shell:
        """
        # combine paired end reads
        extract_markers.py \
        --database {input} \
        --clades {wildcards.species} \
        --output_dir {params.out_dir}
        """


# run strainphlan
rule strainphlan:
    input:
        pkl=expand(
            results + "04_READ_BASED_STRAIN/01_consensus_markers/{sample}.pkl",
            sample=samples,
        ),
        db_markers=resources + "strainphlan/markers/{species}.fna",
        spa_db=resources + "metaphlan/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl",
    output:
        results + "04_READ_BASED_STRAIN/03_strainphlan/{species}.info",
    params:
        out_dir=results + "04_READ_BASED_STRAIN/03_strainphlan/",
        species=config["read_strain"]["strainphlan_species"],
    # conda:
    #     "../envs/humann:3.6--pyh7cba7a3_1.yml"
    container:
        "docker://quay.io/biocontainers/humann:3.6--pyh7cba7a3_1"
    benchmark:
        "benchmark/04_READ_BASED_STRAIN/strainphlan_{species}.tsv"
    resources:
        runtime=config["read_strain"]["strainphlan_runtime"],
        mem_mb=config["read_strain"]["strainphlan_memory"],
    threads: config["read_strain"]["strainphlan_threads"]
    shell:
        """
        # combine paired end reads
        strainphlan \
        --database {input.spa_db} \
        -s {input.pkl} \
        -m {input.db_markers} \
        -o {params.out_dir} \
        -n {threads} \
        -c {params.species} \
        --mutation_rates
        """
