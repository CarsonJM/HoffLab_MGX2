# -------------------------------------
# Read-Based Taxonomy Module
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
# Read-Based Taxonomy Rules
# -------------------------------------
localrules: download_metaphlan_db, merge_metaphlan_profiles, metaphlan_count_features, metaphlan_taxa_reduction


# -----------------------------------------------------
# 01 MetaPhlan
# -----------------------------------------------------
# download metaphlan database
rule download_metaphlan_db:
    output:
        resources + "metaphlan/mpa_vJan21_CHOCOPhlAnSGB_202103.1.bt2"
    log:
        results + "00_LOGS/02_download_metaphlan_db.log",
    params:
        mpa_dir=resources + "metaphlan/",
    conda:
        "../envs/metaphlan:4.0.2--pyhca03a8a_0.yml"
    benchmark:
        "benchmark/02_READ_BASED_TAXONOMY/download_metaphlan_db.tsv"
    resources:
        runtime="00:30:00",
        mem_mb="5000",
    threads:
        config["read_taxonomy"]["metaphlan_threads"]
    shell:
        """
        # download metaphlan database
        metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 \
        --nproc {threads} \
        --bowtie2db {params.mpa_dir} > {log} 2>&1
        """


# run metaphlan
rule metaphlan:
    input:
        mpa_db=resources + "metaphlan/mpa_vJan21_CHOCOPhlAnSGB_202103.1.bt2",
        R1=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_1.fastq",
        R2=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_2.fastq",
    output:
        mpa=results + "02_READ_BASED_TAXONOMY/01_metaphlan/{sample}_profile.tsv",
        sam=results + "02_READ_BASED_TAXONOMY/01_metaphlan/{sample}.sam.bz2",
    log:
        results + "00_LOGS/02_metaphlan.{sample}.log",
    params:
        extra_args=config["read_taxonomy"]["metaphlan_arguments"],
        mpa_dir=resources + "metaphlan/",
        index="mpa_vJan21_CHOCOPhlAnSGB_202103"
    conda:
        "../envs/metaphlan:4.0.2--pyhca03a8a_0.yml"
    benchmark:
        "benchmark/02_READ_BASED_TAXONOMY/metaphlan_{sample}.tsv"
    resources:
        runtime="10:00:00",
        mem_mb="25000",
    threads:
        config["read_taxonomy"]["metaphlan_threads"]
    shell:
        """
        # run metaphlan
        metaphlan {input.R1},{input.R2} \
        --input_type fastq \
        --nproc {threads} \
        --index {params.index} \
        --bowtie2db {params.mpa_dir} \
        --no-map \
        --samout {output.sam} \
        --output_file {output.mpa} \
        {params.extra_args} > {log} 2>&1
        """


# combine metaphlan profils across samples
rule merge_metaphlan_profiles:
    input:
        expand(results + "02_READ_BASED_TAXONOMY/01_metaphlan/{sample}_profile.tsv", sample=samples)
    output:
        results + "02_READ_BASED_TAXONOMY/metaphlan_merged_profiles.tsv"
    params:
        in_dir=results + "02_READ_BASED_TAXONOMY/01_metaphlan/"
    conda:
        "../envs/metaphlan:4.0.2--pyhca03a8a_0.yml"
    benchmark:
        "benchmark/02_READ_BASED_TAXONOMY/merge_metaphlan_profiles.tsv"
    resources:
        runtime="00:01:00",
        mem_mb="1000",
    shell:
        """
        # merge metaphlan profiles
        merge_metaphlan_tables.py {input} > {output}
        """


# combine metaphlan profils across samples
rule metaphlan_count_features:
    input:
        results + "02_READ_BASED_TAXONOMY/metaphlan_merged_profiles.tsv"
    output:
        results + "02_READ_BASED_TAXONOMY/metaphlan_species_counts.tsv"
    conda:
        "../envs/metaphlan:4.0.2--pyhca03a8a_0.yml"
    benchmark:
        "benchmark/02_READ_BASED_TAXONOMY/count_metaphlan_features.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # count metaphlan species
        count_features.py --input {input} \
        --output {output} \
        --include s__ \
        --filter t__ \
        --reduce-sample-name
        """


# reduce metaphlan output to each taxonomic level
rule metaphlan_taxa_reduction:
    input:
        results+"02_READ_BASED_TAXONOMY/metaphlan_merged_profiles.tsv"
    output:
        species=results+"02_READ_BASED_TAXONOMY/metaphlan_combined_profiles_species.tsv",
        genus=results+"02_READ_BASED_TAXONOMY/metaphlan_combined_profiles_genus.tsv",
        family=results+"02_READ_BASED_TAXONOMY/metaphlan_combined_profiles_family.tsv",
        order=results+"02_READ_BASED_TAXONOMY/metaphlan_combined_profiles_order.tsv",
        class_df=results+"02_READ_BASED_TAXONOMY/metaphlan_combined_profiles_class.tsv",
        phylum=results+"02_READ_BASED_TAXONOMY/metaphlan_combined_profiles_phylum.tsv",
        kingdom=results+"02_READ_BASED_TAXONOMY/metaphlan_combined_profiles_kingdom.tsv",
    benchmark:
        "benchmark/02_READ_BASED_TAXONOMY/metaphlan_taxa_reduction.tsv"
    resources:
        runtime="00:01:00",
        mem_mb="1000",
    script:
        "../scripts/02_metaphlan_taxa_reduction.py"
