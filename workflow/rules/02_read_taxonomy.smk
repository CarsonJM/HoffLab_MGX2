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
localrules:
    sgb_to_gtdb_taxonomy,
    merge_metaphlan_profiles,
    metaphlan_count_features,
    metaphlan_alpha_diversity,
    metaphlan_beta_diversity,
    metaphlan_taxa_reduction,


# -----------------------------------------------------
# 01 MetaPhlan
# -----------------------------------------------------
# download metaphlan database
rule download_metaphlan_db:
    message:
        "Downloading MetaPhlan4 database mpa_vJan21_CHOCOPhlAnSGB_202103"
    output:
        resources + "metaphlan/mpa_vJan21_CHOCOPhlAnSGB_202103.1.bt2l",
    params:
        mpa_dir=resources + "metaphlan/",
    # conda:
    #     "../envs/humann:3.6--pyh7cba7a3_1.yml"
    container:
        "docker://quay.io/biocontainers/humann:3.6--pyh7cba7a3_1"
    benchmark:
        "benchmark/02_READ_BASED_TAXONOMY/download_metaphlan_db.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="100000",
    threads: config["read_taxonomy"]["metaphlan_threads"]
    shell:
        """
        rm -rf {params.mpa_dir}

        # download metaphlan database
        metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 \
        --nproc {threads} \
        --bowtie2db {params.mpa_dir}
        """


# run metaphlan
rule metaphlan:
    message:
        "Running MetaPhlan4 on {wildcards.sample}"
    input:
        mpa_db=resources + "metaphlan/mpa_vJan21_CHOCOPhlAnSGB_202103.1.bt2l",
        R1=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_1.fastq.gz",
        R2=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_2.fastq.gz",
    output:
        mpa=results + "02_READ_BASED_TAXONOMY/01_metaphlan/{sample}_profile.tsv",
        sam=results + "02_READ_BASED_TAXONOMY/01_metaphlan/{sample}.sam.bz2",
    params:
        extra_args=config["read_taxonomy"]["metaphlan_arguments"],
        mpa_dir=resources + "metaphlan/",
        index="mpa_vJan21_CHOCOPhlAnSGB_202103",
    # conda:
    #     "../envs/humann:3.6--pyh7cba7a3_1.yml"
    container:
        "docker://quay.io/biocontainers/humann:3.6--pyh7cba7a3_1"
    benchmark:
        "benchmark/02_READ_BASED_TAXONOMY/metaphlan_{sample}.tsv"
    resources:
        runtime="10:00:00",
        mem_mb="25000",
    threads: config["read_taxonomy"]["metaphlan_threads"]
    shell:
        """
        # run metaphlan
        metaphlan {input.R1},{input.R2} \
        --input_type fastq \
        --nproc {threads} \
        --index {params.index} \
        --bowtie2db {params.mpa_dir} \
        --no_map \
        --samout {output.sam} \
        --output_file {output.mpa} \
        {params.extra_args}
        """


# run convert SGB to GTDB taxonomy
rule sgb_to_gtdb_taxonomy:
    message:
        "Converting MetaPhlan4 SGB taxonomy to standardized GTDB taxonomy for {wildcards.sample}"
    input:
        results + "02_READ_BASED_TAXONOMY/01_metaphlan/{sample}_profile.tsv",
    output:
        results + "02_READ_BASED_TAXONOMY/01_metaphlan/{sample}_profile_gtdb.tsv",
    # conda:
    #     "../envs/humann:3.6--pyh7cba7a3_1.yml"
    container:
        "docker://quay.io/biocontainers/humann:3.6--pyh7cba7a3_1"
    benchmark:
        "benchmark/02_READ_BASED_TAXONOMY/sgb_to_gtdb_taxonomy_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="5000",
    shell:
        """
        # convert sgb to gtdb taxonomy
        sgb_to_gtdb_profile.py \
        --input {input} \
        --output {output}
        """


# combine metaphlan profils across samples
rule merge_metaphlan_profiles:
    message:
        "Combining MetaPhlan4 results across all samples"
    input:
        sgb=expand(
            results + "02_READ_BASED_TAXONOMY/01_metaphlan/{sample}_profile.tsv",
            sample=samples,
        ),
        gtdb=expand(
            results + "02_READ_BASED_TAXONOMY/01_metaphlan/{sample}_profile_gtdb.tsv",
            sample=samples,
        ),
    output:
        sgb=results + "02_READ_BASED_TAXONOMY/metaphlan_merged_profiles.tsv",
        gtdb=results + "02_READ_BASED_TAXONOMY/metaphlan_merged_profiles_gtdb.tsv",
    params:
        in_dir=results + "02_READ_BASED_TAXONOMY/01_metaphlan/",
    # conda:
    #     "../envs/humann:3.6--pyh7cba7a3_1.yml"
    container:
        "docker://quay.io/biocontainers/humann:3.6--pyh7cba7a3_1"
    benchmark:
        "benchmark/02_READ_BASED_TAXONOMY/merge_metaphlan_profiles.tsv"
    resources:
        runtime="00:01:00",
        mem_mb="1000",
    shell:
        """
        # merge metaphlan profiles
        merge_metaphlan_tables.py {input.sgb} > {output.sgb}

        # merge metaphlan gtdb profiles
        merge_metaphlan_tables.py --gtdb_profiles {input.gtdb} > {output.gtdb}
        """


# -----------------------------------------------------
# 02 Analysis
# -----------------------------------------------------
# combine metaphlan profils across samples
rule metaphlan_count_features:
    message:
        "Counting the number of species MetaPhlan4 profiles"
    input:
        results + "02_READ_BASED_TAXONOMY/metaphlan_merged_profiles.tsv",
    output:
        results + "02_READ_BASED_TAXONOMY/metaphlan_species_counts.tsv",
    params:
        count_features_script="workflow/scripts/02_metaphlan_count_features.py",
    # conda:
    #     "../envs/humann:3.6--pyh7cba7a3_1.yml"
    container:
        "docker://quay.io/biocontainers/humann:3.6--pyh7cba7a3_1"
    benchmark:
        "benchmark/02_READ_BASED_TAXONOMY/count_metaphlan_features.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # count metaphlan species
        python {params.count_features_script} --input {input} \
        --output {output} \
        --include s__ \
        --filter t__ \
        --reduce-sample-name
        """


# calculate alpha diversity using metaphlan scripts
rule metaphlan_alpha_diversity:
    message:
        (
            "Calculating alpha diversity using "
            + config["read_taxonomy"]["alpha_diversity"]
        )
    input:
        results + "02_READ_BASED_TAXONOMY/metaphlan_merged_profiles.tsv",
    output:
        results
        + "02_READ_BASED_TAXONOMY/02_diversity_metrics/metaphlan_merged_profiles_"
        + config["read_taxonomy"]["alpha_diversity"]
        + ".tsv",
    params:
        diversity_script="workflow/scripts/02_metaphlan_diversity.R",
        alpha_metric=config["read_taxonomy"]["alpha_diversity"],
        out_dir=results + "02_READ_BASED_TAXONOMY/02_diversity_metrics",
    container:
        "/gscratch/pedslabs/hofflab/carsonjm/apptainer/metaphlan_diversity.sif"
    benchmark:
        "benchmark/02_READ_BASED_TAXONOMY/metaphlan_alpha_diversity.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # calculate alpha diversity using MetaPhlan4 script
        Rscript {params.diversity_script} \
        --file {input} \
        --out_directory {params.out_dir} \
        --diversity alpha \
        --metric {params.alpha_metric} \
        --taxon_separator t__
        """


# calculate beta diversity using metaphlan scripts
rule metaphlan_beta_diversity:
    message:
        (
            "Calculating beta diversity using "
            + config["read_taxonomy"]["beta_diversity"]
        )
    input:
        results + "02_READ_BASED_TAXONOMY/metaphlan_merged_profiles.tsv",
    output:
        results
        + "02_READ_BASED_TAXONOMY/02_diversity_metrics/metaphlan_merged_profiles_"
        + config["read_taxonomy"]["beta_diversity"]
        + ".tsv",
    params:
        diversity_script="workflow/scripts/02_metaphlan_diversity.R",
        beta_metric=config["read_taxonomy"]["beta_diversity"],
        out_dir=results + "02_READ_BASED_TAXONOMY/02_diversity_metrics",
    container:
        "/gscratch/pedslabs/hofflab/carsonjm/apptainer/metaphlan_diversity.sif"
    benchmark:
        "benchmark/02_READ_BASED_TAXONOMY/metaphlan_beta_diversity.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # calculate beta diversity using MetaPhlan4 script
        Rscript {params.diversity_script} \
        --file {input} \
        --out_directory {params.out_dir} \
        --diversity beta \
        --metric {params.beta_metric} \
        --taxon_separator t__
        """


# reduce metaphlan output to each taxonomic level
rule metaphlan_taxa_reduction:
    message:
        "Organizing MetaPhlan4 outputs by taxonomic level"
    input:
        results + "02_READ_BASED_TAXONOMY/metaphlan_merged_profiles.tsv",
    output:
        species=results + "02_READ_BASED_TAXONOMY/metaphlan_merged_profiles_species.tsv",
        genus=results + "02_READ_BASED_TAXONOMY/metaphlan_merged_profiles_genus.tsv",
        phylum=results + "02_READ_BASED_TAXONOMY/metaphlan_merged_profiles_phylum.tsv",
    # conda:
    #     "../envs/scripts.yml"
    container:
        "/gscratch/pedslabs/hofflab/carsonjm/apptainer/scripts.sif"
    benchmark:
        "benchmark/02_READ_BASED_TAXONOMY/metaphlan_taxa_reduction.tsv"
    resources:
        runtime="00:01:00",
        mem_mb="1000",
    script:
        "../scripts/02_metaphlan_taxa_reduction.py"
