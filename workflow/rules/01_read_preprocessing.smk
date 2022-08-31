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
# -----------------------------------------------------
# 00 Symlink Reads
# -----------------------------------------------------
from snakemake.remote import AUTO

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
        partition="ckpt",
        mem_mb="1000",
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
        R1=results + "01_READ_PREPROCESSING/01_merge_replicates/{sample}_1.fastq.gz",
        R2=results + "01_READ_PREPROCESSING/01_merge_replicates/{sample}_2.fastq.gz",
    benchmark:
        "benchmark/01_READ_PREPROCESSING/merge_replicates_{sample}.tsv"
    resources:
        runtime="00:00:10",
        partition="ckpt",
        mem_mb="1000",
    shell:
        """
        # symlink input paths to renamed files
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """


# -----------------------------------------------------
# 02 Raw FASTQC & MULTIQC
# -----------------------------------------------------
# Run FASTQC on raw R1 reads
rule raw_fastqc_R1:
    input:
        results + "01_READ_PREPROCESSING/01_merge_replicates/{sample}_1.fastq.gz",
    output:
        results + "01_READ_PREPROCESSING/02_raw_fastqc/{sample}_1_fastqc.html",
    log:
        results + "00_LOGS/01_raw_fastqc_R1.{sample}.log",
    params:
        out_dir=results + "01_READ_PREPROCESSING/02_raw_fastqc/",
    conda:
        "../envs/fastqc.yml"
    container: 
        "docker://biocontainers/fastqc"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/raw_fastqc_R1_{sample}.tsv"
    resources:
        runtime="00:05:00",
        partition="ckpt",
        mem_mb="1000",
    shell:
        """
        # generate fastqc report for forward reads
        fastqc {input} --outdir {params.out_dir} > {log} 2>&1
        """


# Run FASTQC on raw R2 reads
rule raw_fastqc_R2:
    input:
        results + "01_READ_PREPROCESSING/01_merge_replicates/{sample}_2.fastq.gz",
    output:
        results + "01_READ_PREPROCESSING/02_raw_fastqc/{sample}_2_fastqc.html",
    log:
        results + "00_LOGS/01_raw_fastqc_R2.{sample}.log",
    params:
        out_dir=results + "01_READ_PREPROCESSING/02_raw_fastqc/",
    conda:
        "../envs/fastqc.yml"
    container: 
        "docker://biocontainers/fastqc"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/raw_fastqc_R2_{sample}.tsv"
    resources:
        runtime="00:05:00",
        partition="ckpt",
        mem_mb="1000",
    shell:
        """
        # generate fastqc report for reverse reads
        fastqc {input} --outdir {params.out_dir} > {log} 2>&1
        """


# Combine FASTQC results using multiQC
rule raw_multiqc:
    input:
        expand(
            results
            + "01_READ_PREPROCESSING/02_raw_fastqc/{sample}_{read}_fastqc.html",
            sample=samples,
            read=["1", "2"],
        ),
    output:
        report(
            results + "01_READ_PREPROCESSING/raw_multiqc_report.html",
            caption="../report/01_read_preprocessing_analysis.rst",
            category="Step 01: Read preprocessing",
        ),
    log:
        results + "00_LOGS/01_raw_multiqc.log",
    params:
        in_dir=results + "01_READ_PREPROCESSING/02_raw_fastqc/",
        out_dir=results + "01_READ_PREPROCESSING/single/",
        out_file=results + "01_READ_PREPROCESSING/single/multiqc_report.html",
    conda:
        "../envs/multiqc.yml"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/raw_multiqc.tsv"
    resources:
        runtime="00:10:00",
        partition="ckpt",
        mem_mb="1000",
    shell:
        """
        # remove multiqc output
        rm -f {params.out_file}

        # generate mutliqc report
        multiqc {params.in_dir} -o {params.out_dir} > {log} 2>&1
        mv {params.out_file} {output}

        rm -r {params.out_dir}
        """


# -----------------------------------------------------
# 03 KneadData
# -----------------------------------------------------
# build kneaddata bowtie2 database
rule download_kneaddata_database:
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
        "../envs/kneaddata.yml"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/download_kneaddata.tsv"
    resources:
        runtime="00:10:00",
        partition="ckpt",
        mem_mb="10000",
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
        R1=results + "01_READ_PREPROCESSING/01_merge_replicates/{sample}_1.fastq.gz",
        R2=results + "01_READ_PREPROCESSING/01_merge_replicates/{sample}_2.fastq.gz",
    output:
        log=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}.log",
        R1=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_1.fastq",
        R2=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_2.fastq",
    log:
        results + "00_LOGS/01_kneaddata.{sample}.log",
    params:
        out_dir=results + "01_READ_PREPROCESSING/03_kneaddata/",
        human_db=resources + "kneaddata/",
        extra_args=config["read_preprocessing"]["kneaddata_arguments"],
        prefix="{sample}",
    conda:
        "../envs/kneaddata.yml"
    threads: config["read_preprocessing"]["kneaddata_threads"]
    benchmark:
        "benchmark/01_READ_PREPROCESSING/kneaddata_{sample}.tsv"
    resources:
        runtime="01:00:00",
        partition="ckpt",
        mem_mb="10000",
    shell:
        """
        # run kneaddata to quality filter and remove host reads
        kneaddata --input {input.R1} --input {input.R2} \
        --output {params.out_dir} \
        --output-prefix {params.prefix} \
        --reference-db {params.human_db} \
        --threads {threads} \
        {params.extra_args}

        # copy log output to log directory
        cp {output.log} {log}
        """


# -----------------------------------------------------
# 04 Preprocessed FASTQC & MULTIQC
# -----------------------------------------------------
# Run FASTQC on preprocessed R1 reads
rule preprocessed_fastqc_R1:
    input:
        results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_1.fastq",
    output:
        results
        + "01_READ_PREPROCESSING/04_preprocessed_fastqc/{sample}_paired_1_fastqc.html",
    log:
        results + "00_LOGS/01_preprocessed_fastqc_R1.{sample}.log",
    params:
        out_dir=results + "01_READ_PREPROCESSING/04_preprocessed_fastqc/",
    conda:
        "../envs/fastqc.yml"
    container: 
        "docker://biocontainers/fastqc"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/preprocessed_fastqc_R1_{sample}.tsv"
    resources:
        runtime="00:05:00",
        partition="ckpt",
        mem_mb="1000",
    shell:
        """
        # generate fastqc report for forward reads
        fastqc {input} --outdir {params.out_dir} > {log} 2>&1
        """


# Run FASTQC on preprocessed R2 reads
rule preprocessed_fastqc_R2:
    input:
        results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_2.fastq",
    output:
        results
        + "01_READ_PREPROCESSING/04_preprocessed_fastqc/{sample}_paired_2_fastqc.html",
    log:
        results + "00_LOGS/01_preprocessed_fastqc_R2.{sample}.log",
    params:
        out_dir=results + "01_READ_PREPROCESSING/04_preprocessed_fastqc/",
    conda:
        "../envs/fastqc.yml"
    container: 
        "docker://biocontainers/fastqc"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/preprocessed_fastqc_R2_{sample}.tsv"
    resources:
        runtime="00:05:00",
        partition="ckpt",
        mem_mb="1000",
    shell:
        """
        # generate fastqc report for reverse reads
        fastqc {input} --outdir {params.out_dir} > {log} 2>&1
        """


# Combine FASTQC results using multiQC
rule preprocessed_multiqc:
    input:
        expand(
            results
            + "01_READ_PREPROCESSING/04_preprocessed_fastqc/{sample}_paired_{read}_fastqc.html",
            sample=samples,
            read=["1", "2"],
        ),
    output:
        report(
            results + "01_READ_PREPROCESSING/preprocessed_multiqc_report.html",
            caption="../report/01_read_preprocessing_analysis.rst",
            category="Step 01: Read preprocessing",
        ),
    log:
        results + "00_LOGS/01_preprocessed_multiqc.log",
    params:
        in_dir=results + "01_READ_PREPROCESSING/04_preprocessed_fastqc/",
        out_file=results
        + "01_READ_PREPROCESSING/04_preprocessed_fastqc/multiqc_report.html",
    conda:
        "../envs/multiqc.yml"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/preprocessed_multiqc.tsv"
    resources:
        runtime="00:10:00",
        partition="ckpt",
        mem_mb="1000",
    shell:
        """
        # remove multiqc output
        rm -f {params.out_file}

        # generate mutliqc report
        multiqc {params.in_dir} -o {params.in_dir} > {log} 2>&1
        mv {params.out_file} {output}
        """
