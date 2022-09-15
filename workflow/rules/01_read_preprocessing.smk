# -------------------------------------
# Read Preprocessing Module
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
# Preprocessing Rules
# -------------------------------------
localrules: symlink_reads, merge_replicates, download_kneaddata_db, kneaddata_read_counts, combine_read_counts


# -----------------------------------------------------
# 00 Symlink Reads
# -----------------------------------------------------
# symlink input paths to new paths
rule symlink_reads:
    input:
        R1=lambda wildcards: samples_df[
            (+samples_df["sample"] + "_" + samples_df["replicate"])
            == wildcards.sample_replicate
        ]["R1"].iloc[0],
        R2=lambda wildcards: samples_df[
            (+samples_df["sample"] + "_" + samples_df["replicate"])
            == wildcards.sample_replicate
        ]["R2"].iloc[0],
    output:
        R1=results + "00_INPUT/{sample_replicate}_1.fastq.gz",
        R2=results + "00_INPUT/{sample_replicate}_2.fastq.gz",
    benchmark:
        "benchmark/01_READ_PREPROCESSING/symlink_reads_{sample_replicate}.tsv"
    resources:
        runtime="00:00:10",
        mem_mb="100",
    shell:
        """
        # symlink input paths to renamed files
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """


# -----------------------------------------------------
# 01 Merge Replicates
# -----------------------------------------------------
# identify replicates
sample_replicate = samples_df[["sample", "replicate"]]
sample_replicate_dictionary = sample_replicate.set_index("sample").to_dict()[
    "replicate"
]


# merge replicate files into single file
rule merge_replicates:
    input:
        R1=lambda wildcards: expand(
            results + "00_INPUT/{{sample}}_{replicate}_1.fastq.gz",
            replicate=sample_replicate_dictionary[wildcards.sample],
        ),
        R2=lambda wildcards: expand(
            results + "00_INPUT/{{sample}}_{replicate}_2.fastq.gz",
            replicate=sample_replicate_dictionary[wildcards.sample],
        ),
    output:
        R1=temp(results + "01_READ_PREPROCESSING/01_merge_replicates/{sample}.R1.fastq.gz"),
        R2=temp(results + "01_READ_PREPROCESSING/01_merge_replicates/{sample}.R2.fastq.gz"),
    benchmark:
        "benchmark/01_READ_PREPROCESSING/merge_replicates_{sample}.tsv"
    resources:
        runtime="00:00:10",
        mem_mb="100",
    shell:
        """
        # symlink input paths to renamed files
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """


# -----------------------------------------------------
# 02 Clumpify
# -----------------------------------------------------
# deduplicate reads with clumpify
rule clumpify:
    input:
        R1=results + "01_READ_PREPROCESSING/01_merge_replicates/{sample}.R1.fastq.gz",
        R2=results + "01_READ_PREPROCESSING/01_merge_replicates/{sample}.R2.fastq.gz",
    output:
        R1=results + "01_READ_PREPROCESSING/02_clumpify/{sample}.R1.fastq",
        R2=results + "01_READ_PREPROCESSING/02_clumpify/{sample}.R2.fastq",
        log=results + "01_READ_PREPROCESSING/02_clumpify/{sample}.log",
    log:
        results + "00_LOGS/01_clumpify.{sample}.log"
    params:
        extra_args=config['read_preprocessing']['clumpify_arguments'],
        R1=results + "01_READ_PREPROCESSING/01_merge_replicates/{sample}.R1.fastq",
        R2=results + "01_READ_PREPROCESSING/01_merge_replicates/{sample}.R2.fastq",
        log_dir=results +"00_LOGS/",
    conda:
        "../envs/bbmap.yml"
    container:
        "docker://quay.io/biocontainers/bbmap:38.95--he522d1c_0"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/clumpify_{sample}.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="5000",
    shell:
        """
        # unzip input files
        gunzip {input.R1}
        gunzip {input.R2}

        # run clumpify
        clumpify.sh \
        in={params.R1} \
        in2={params.R2} \
        out={output.R1} \
        out2={output.R2} \
        {params.extra_args} > {log} 2>&1

        # copy log to log dir
        cp {log} {output.log}
        """

# -----------------------------------------------------
# 03 KneadData
# -----------------------------------------------------
# download biobakery workflow databases
rule download_kneaddata_db:
    output:
        resources + "kneaddata/hg37dec_v0.1.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.2.bt2",
        resources + "kneaddata/hg37dec_v0.1.3.bt2",
        resources + "kneaddata/hg37dec_v0.1.4.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.2.bt2",
    log:
        results + "00_LOGS/01_download_kneaddata_database.log",
    params:
        kneaddata_db=resources + "kneaddata/",
    conda:
        "../envs/biobakery_workflows.yml"
    container: 
        "/gscratch/stf/carsonjm/apptainer/workflows_3.0.0.a.7.sif"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/download_kneaddata.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # download human genome reference to desired directory
        kneaddata_database --download human_genome bowtie2 {params.kneaddata_db} > {log} 2>&1
        """


# Quality filter and remove human reads with kneaddata
rule kneaddata:
    input:
        resources + "kneaddata/hg37dec_v0.1.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.2.bt2",
        resources + "kneaddata/hg37dec_v0.1.3.bt2",
        resources + "kneaddata/hg37dec_v0.1.4.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.2.bt2",
        R1=results + "01_READ_PREPROCESSING/02_clumpify/{sample}.R1.fastq",
        R2=results + "01_READ_PREPROCESSING/02_clumpify/{sample}.R2.fastq",
    output:
        R1=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_1.fastq",
        R2=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_2.fastq",
        log=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}.log",
    log:
        results + "00_LOGS/01_kneaddata.{sample}.log"
    params:
        out_dir=results + "01_READ_PREPROCESSING/03_kneaddata/",
        human_db=resources + "kneaddata/",
        extra_args=config["read_preprocessing"]["kneaddata_arguments"],
        prefix="{sample}",
        log_dir=results + "00_LOGS/",
    conda:
        "../envs/biobakery_workflows.yml"
    container: 
        "/gscratch/stf/carsonjm/apptainer/workflows_3.0.0.a.7.sif"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/kneaddata_{sample}.tsv"
    resources:
        runtime="4:00:00",
        mem_mb="5000",
    threads: config["read_preprocessing"]["kneaddata_threads"]
    shell:
        """
        # run kneaddata to quality filter and remove host reads
        kneaddata --input {input.R1} --input {input.R2} \
        --output {params.out_dir} \
        --output-prefix {params.prefix} \
        --reference-db {params.human_db} \
        --threads {threads} \
        --run-trf \
        --serial \
        --log {log} \
        {params.extra_args}

        # copy log to log dir
        cp {log} {output.log}
        """


# -----------------------------------------------------
# 04 Read Counts
# -----------------------------------------------------
# determine clumpify read counts
rule clumpify_read_counts:
    input:
        expand(results + "01_READ_PREPROCESSING/02_clumpify/{sample}.log", sample=samples),
    output:
        results + "01_READ_PREPROCESSING/02_clumpify/read_counts.tsv",
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    notebook:
        "../notebooks/01_clumpify_read_counts.py.ipynb"


# determine read counts using kneaddata utils
rule kneaddata_read_counts:
    input:
        expand(results + "01_READ_PREPROCESSING/03_kneaddata/{sample}.log", sample=samples),
    output:
        results + "01_READ_PREPROCESSING/03_kneaddata/read_counts.tsv",
    params:
        log_dir=results + "01_READ_PREPROCESSING/03_kneaddata/",
    conda:
        "../envs/biobakery_workflows.yml"
    container: 
        "/gscratch/stf/carsonjm/apptainer/workflows_3.0.0.a.7.sif"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/kneaddata_read_counts.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # generate read counts from kneaddata log files
        kneaddata_read_count_table \
        --input {params.log_dir} \
        --output {output}
        """


# combine clumpify and kneaddata read counts
rule combine_read_counts:
    input:
        clumpify=results + "01_READ_PREPROCESSING/02_clumpify/read_counts.tsv",
        kneaddata=results + "01_READ_PREPROCESSING/03_kneaddata/read_counts.tsv",
    output:
        results + "01_READ_PREPROCESSING/read_preprocessing_report.tsv",
    benchmark:
        "benchmark/01_READ_PREPROCESSING/combine_read_counts.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    notebook:
        "../notebooks/01_combine_read_counts.py.ipynb"
