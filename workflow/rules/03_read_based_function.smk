# -------------------------------------
# Read-Based Function Module
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
# Read-Based Function Rules
# -------------------------------------
# -----------------------------------------------------
# 01 HUManN
# -----------------------------------------------------
# download humann database
rule download_humann_db:
    output:
        resources + "uniref/uniref90_201901.dmnd"
    log:
        results + "00_LOGS/03_download_humann_db.log",
    params:
        humann_dir=resources + "humann/",
    conda:
        "../envs/biobakery_workflows.yml"
    container: 
        "docker://biobakery/workflows:3.0.0.a.7"
    benchmark:
        "benchmark/03_READ_BASED_FUNCTION/download_humann_db.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="10000",
    shell:
        """
        # download humann chocophlan database
        humann_databases --download chocophlan full {params.humann_dir}

        # download humann uniref90 database
        humann_databases --download uniref uniref90_diamond {params.humann_dir}

        # download humann utility mapping database
        humann_databases --download utility_mapping full {params.humann_dir}
        """


# merge read pairs
rule merge_read_pairs:
    input:
        R1=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_1.fastq",
        R2=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_2.fastq",
    output:
        results + "03_READ_BASED_FUNCTION/01_merge_pairs/{sample}.fastq"
    log:
        results + "00_LOGS/03_merge_pairs.{sample}.log"
    benchmark:
        "benchmark/03_READ_BASED_FUNCTION/merge_pairs_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="5000",
    shell:
        """
        # merge read pairs for kneaddata
        cat {input.R1} {input.R2} > {output}
        """


# run humann on all samples
rule humann:
    input:
        seq=results + "03_READ_BASED_FUNCTION/01_merge_pairs/{sample}.fastq",
        mpa=results + "02_READ_BASED_TAXONOMY/01_metaphlan/{sample}_profile.tsv",
        humann_db=resources + "uniref/uniref90_201901.dmnd"
    output:
        gf=results + "03_READ_BASED_FUNCTION/02_humann/{sample}_genefamilies.tsv",
        pa=results + "03_READ_BASED_FUNCTION/02_humann/{sample}_pathabundance.tsv",
        pc=results + "03_READ_BASED_FUNCTION/02_humann/{sample}_pathcoverage.tsv",
        log=results + "03_READ_BASED_FUNCTION/03_humann_logs/{sample}.log",
    log:
        results + "00_LOGS/03_humann.{sample}.log",
    params:
        out_dir=results + "03_READ_BASED_FUNCTION/02_humann/",
        cpa_dir=resources + "humann/chocophlan/",
        uniref_dir=resources + "humann/uniref/",
        map_dir=resources + "humann/utility_mapping/"
    conda:
        "../envs/biobakery_workflows.yml"
    container: 
        "docker://biobakery/workflows:3.0.0.a.7"
    benchmark:
        "benchmark/03_READ_BASED_FUNCTION/humann_{sample}.tsv"
    resources:
        runtime="48:00:00",
        mem_mb="50000",
    threads:
        config["read_function"]["humann_threads"]
    shell:
        """
        # run humann on samples
        humann --input {input.seq} \
        --output {params.out_dir} \
        --threads {threads} \
        --taxonomic-profile {input.mpa} \
        --input_format fastq \
        --o-log {output.log} \
        --nucleotide-database {params.cpa_dir} \
        --protein-database {params.uniref_dir} \
        --remove-temp-output \
        {params.extra_args}

        # copy log to logs dir
        cp {output.log} {log}
        """


# get counts from humann
rule get_counts_from_humann_logs:
    input:
        expand(results + "03_READ_BASED_FUNCTION/03_humann_logs/{sample}.log", sample=samples)
    output:
        results + "03_READ_BASED_FUNCTION/humann_read_alignments.tsv"
    params:
        in_dir=results + "03_READ_BASED_FUNCTION/03_humann_logs/"
    conda:
        "../envs/biobakery_workflows.yml"
    container: 
        "docker://biobakery/workflows:3.0.0.a.7"
    benchmark:
        "benchmark/03_READ_BASED_FUNCTION/get_counts_from_humann_logs.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # count reads aligned to species
        get_counts_from_humann_logs.py --input {params.in_dir} \
        --output {output}
        """


# group gene families at KO level
rule humann_regroup_gene_families:
    input:
        results + "03_READ_BASED_FUNCTION/02_humann/{sample}_genefamilies.tsv",
    output:
        results + "03_READ_FUNCTION/02_humann/{sample}_ecs.tsv"
    log:
        results + "00_LOGS/03_humann_regroup_table.{sample}.log",
    conda:
        "../envs/biobakery_workflows.yml"
    container: 
        "docker://biobakery/workflows:3.0.0.a.7"
    benchmark:
        "benchmark/03_READ_FUNCTION/humann_regoup_table_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="5000",
    shell:
        """
        # regroup gene families
        humann_regroup_table --input {input} \
        --output {output} \
        --groups uniref90_level4ec
        """


# merge unnormalized tables
rule humann_join_unnormalized_tables:
    input:
        gf=expand(results + "03_READ_BASED_FUNCTION/02_humann/{sample}_genefamilies.tsv", sample=samples),
        ecs=expand(results + "03_READ_FUNCTION/02_humann/{sample}_ecs.tsv", sample=samples),
        pa=expand(results + "03_READ_BASED_FUNCTION/02_humann/{sample}_pathabundance.tsv", sample=samples),
    output:
        gf=results + "03_READ_BASED_FUNCTION/unnormalized_genefamilies.tsv",
        ecs=results + "03_READ_FUNCTION/unnormalized_ecs.tsv",
        pa=results + "03_READ_BASED_FUNCTION/unnormalized_pathabundance.tsv",
    params:
        in_dir=results + "03_READ_BASED_FUNCTION/02_humann/"
    log:
        results + "00_LOGS/03_humann_join_unnormalized.log",
    conda:
        "../envs/biobakery_workflows.yml"
    container: 
        "docker://biobakery/workflows:3.0.0.a.7"
    benchmark:
        "benchmark/03_READ_BASED_FUNCTION/humann_join_unnormalized.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # run join tables
        humann_join_tables --input {params.in_dir} \
        --output {output.gf} \
        --file-name genefamilies

        # run join tables
        humann_join_tables --input {params.in_dir} \
        --output {output.ecs} \
        --file-name ecs

        # run join tables
        humann_join_tables --input {params.in_dir} \
        --output {output.pa} \
        --file-name pathabundance
        """


# humann normalize tables
rule humann_renorm_tables:
    input:
        gf=results + "03_READ_BASED_FUNCTION/02_humann/{sample}_genefamilies.tsv",
        ecs=results + "03_READ_FUNCTION/02_humann/{sample}_ecs.tsv",
        pa=results + "03_READ_BASED_FUNCTION/02_humann/{sample}_pathabundance.tsv",
    output:
        gf=results + "03_READ_BASED_FUNCTION/04_humann_renorm/{sample}_genefamilies_renorm.tsv",
        ecs=results + "03_READ_FUNCTION/04_humann_renorm/{sample}_ecs_renorm.tsv",
        pa=results + "03_READ_BASED_FUNCTION/04_humann_renorm/{sample}_pathabundance_renorm.tsv",
    log:
        results + "00_LOGS/03_humann_renorm.log",
    conda:
        "../envs/biobakery_workflows.yml"
    container: 
        "docker://biobakery/workflows:3.0.0.a.7"
    benchmark:
        "benchmark/03_READ_BASED_FUNCTION/humann_renorm_{sample}.tsv"
    resources:
        runtime="01:00:00",
        mem_mb="5000",
    shell:
        """
        # run renorm tables
        humann_renorm_table --input {input.gf} \
        --output {output.gf} \
        --units relab \
        --special n

        # run renorm tables
        humann_renorm_table --input {input.ecs} \
        --output {output.ecs} \
        --units relab \
        --special n

        # run renorm tables
        humann_renorm_table --input {input.pa} \
        --output {output.pa} \
        --units relab \
        --special n
        """


# merge unnormalized tables
rule humann_join_normalized_tables:
    input:
        gf=expand(results + "03_READ_BASED_FUNCTION/04_humann_renorm/{sample}_genefamilies_renorm.tsv", sample=samples),
        ecs=expand(results + "03_READ_FUNCTION/04_humann_renorm/{sample}_ecs_renorm.tsv", sample=samples),
        pa=expand(results + "03_READ_BASED_FUNCTION/04_humann_renorm/{sample}_pathabundance_renorm.tsv", sample=samples),
    output:
        gf=results + "03_READ_BASED_FUNCTION/genefamilies_relab.tsv",
        ecs=results + "03_READ_FUNCTION/ecs_relab.tsv",
        pa=results + "03_READ_BASED_FUNCTION/pathabundance_relab.tsv",
    params:
        in_dir=results + "03_READ_BASED_FUNCTION/04_humann_renorm/"
    log:
        results + "00_LOGS/03_humann_join_normalized.log",
    conda:
        "../envs/biobakery_workflows.yml"
    container: 
        "docker://biobakery/workflows:3.0.0.a.7"
    benchmark:
        "benchmark/03_READ_BASED_FUNCTION/humann_join_normalized.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # run join tables
        humann_join_tables --input {params.in_dir} \
        --output {output.gf} \
        --file-name genefamilies

        # run join tables
        humann_join_tables --input {params.in_dir} \
        --output {output.ecs} \
        --file-name ecs

        # run join tables
        humann_join_tables --input {params.in_dir} \
        --output {output.pa} \
        --file-name pathabundance
        """


# determine feature counts
rule humann_count_features:
    input:
        gf=results + "03_READ_BASED_FUNCTION/genefamilies_relab.tsv",
        ecs=results + "03_READ_FUNCTION/ecs_relab.tsv",
        pa=results + "03_READ_BASED_FUNCTION/pathabundance_relab.tsv",
    output:
        gf=results + "03_READ_BASED_FUNCTION/05_humann_feature_counts/genefamilies_relab_counts.tsv",
        ecs=results + "03_READ_FUNCTION/05_humann_feature_counts/ecs_relab_counts.tsv",
        pa=results + "03_READ_BASED_FUNCTION/05_humann_feature_counts/pathabundance_relab_counts.tsv",
    log:
        results + "00_LOGS/03_humann_count_features.log",
    conda:
        "../envs/biobakery_workflows.yml"
    container: 
        "docker://biobakery/workflows:3.0.0.a.7"
    benchmark:
        "benchmark/03_READ_BASED_FUNCTION/humann_count_features.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # count features
        count_features.py --input {params.input.gf} \
        --output {output.gf} \
        --reduce-sample-name \
        --ignore-un-features \
        --ignore-stratification

        # count features
        count_features.py --input {params.input.ecs} \
        --output {output.ecs} \
        --reduce-sample-name \
        --ignore-un-features \
        --ignore-stratification

        # count features
        count_features.py --input {params.input.pa} \
        --output {output.pa} \
        --reduce-sample-name \
        --ignore-un-features \
        --ignore-stratification
        """


# merge feature tables
rule humann_join_feature_counts_tables:
    input:
        gf=results + "03_READ_BASED_FUNCTION/05_humann_feature_counts/genefamilies_relab_counts.tsv",
        ecs=results + "03_READ_FUNCTION/05_humann_feature_counts/ecs_relab_counts.tsv",
        pa=results + "03_READ_BASED_FUNCTION/05_humann_feature_counts/pathabundance_relab_counts.tsv",
    output:
        results + "03_READ_BASED_FUNCTION/merged_relab_counts.tsv",
    params:
        in_dir=results + "03_READ_BASED_FUNCTION/05_humann_feature_counts/"
    log:
        results + "00_LOGS/03_humann_join_feature_counts_tables.log",
    conda:
        "../envs/biobakery_workflows.yml"
    container: 
        "docker://biobakery/workflows:3.0.0.a.7"
    benchmark:
        "benchmark/03_READ_BASED_FUNCTION/humann_join_feature_counts_tables.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # run join tables
        humann_join_tables --input {params.in_dir} \
        --output {output} \
        --file-name relab_counts
        """