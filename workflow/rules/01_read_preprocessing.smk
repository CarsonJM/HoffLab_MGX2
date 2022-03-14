# -------------------------------------
# Read preprocessing (only runs if config['input_type'] == "reads")
# -------------------------------------
import pandas as pd
import os


# Load sample information and validate
configfile: "config/config.yaml"


samples_df = pd.read_csv(config["sample_info"]["sample_list_path"], sep="\t")

# get current working directory so absolute paths can be used for input/output files
results = os.getcwd()


# load report
report: "report/workflow.rst"


# load resources folder path
resources = config["resources_path"]

# load sample information to be used in workflow
samples = samples_df["sample"]
R1_files = samples_df["R1"]
R2_files = samples_df["R2"]
contig_files = samples_df["contigs"]

# -------------------------------------
# Read preprocessing rules
# -------------------------------------


# create a system link for input reads (makes downstream processing much easier)
rule link_reads:
    input:
        R1=expand("{R1}", R1=R1_files),
        R2=expand("{R2}", R2=R2_files),
    output:
        R1=results + "/00_INPUT_DATA/01_reads/{sample}_R1.fastq.gz",
        R2=results + "/00_INPUT_DATA/01_reads/{sample}_R2.fastq.gz",
    params:
        R1=lambda wildcards: samples_df[samples_df["sample"] == wildcards.sample].iloc[
            0
        ]["R1"],
        R2=lambda wildcards: samples_df[samples_df["sample"] == wildcards.sample].iloc[
            0
        ]["R2"],
    shell:
        """
        ln -s {params.R1} {output.R1}
        ln -s {params.R2} {output.R2}
        """


# build kneaddata bowtie2 database
rule build_kneaddata_database:
    output:
        resources + "kneaddata/hg37dec_v0.1.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.2.bt2",
        resources + "kneaddata/hg37dec_v0.1.3.bt2",
        resources + "kneaddata/hg37dec_v0.1.4.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.2.bt2",
    params:
        human_db_dir=resources + "kneaddata/",
    threads: 1
    conda:
        "../envs/kneaddata.yml"
    shell:
        """
        # download human genome reference to desired directory
        kneaddata_database --download human_genome bowtie2 {params.human_db_dir}
        """


# quality filter reads with kneaddata
rule kneaddata:
    input:
        resources + "kneaddata/hg37dec_v0.1.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.2.bt2",
        resources + "kneaddata/hg37dec_v0.1.3.bt2",
        resources + "kneaddata/hg37dec_v0.1.4.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.2.bt2",
        R1=results + "/00_INPUT_DATA/01_reads/{sample}_R1.fastq.gz",
        R2=results + "/00_INPUT_DATA/01_reads/{sample}_R2.fastq.gz",
    output:
        kd_R1=results
        + "/01_READ_PREPROCESSING/01_kneaddata/{sample}/kd_paired_1.fastq.gz",
        kd_R2=results
        + "/01_READ_PREPROCESSING/01_kneaddata/{sample}/kd_paired_2.fastq.gz",
        kd_R1S=results
        + "/01_READ_PREPROCESSING/01_kneaddata/{sample}/kd_unmatched_1.fastq.gz",
        kd_R2S=results
        + "/01_READ_PREPROCESSING/01_kneaddata/{sample}/kd_unmatched_2.fastq.gz",
    params:
        kd_dir=directory(results + "/01_READ_PREPROCESSING/01_kneaddata/{sample}/"),
        human_db=resources + "kneaddata/",
        extra_args=config["kneaddata"]["extra_arguments"],
    log:
        results + "/00_LOGS/01_read_preprocessing_{sample}.kneaddata.log",
    threads: 10
    priority: 2
    conda:
        "../envs/kneaddata.yml"
    shell:
        """
        # run kneaddata to quality filter and remove host reads
        kneaddata --input {input.R1} --input {input.R2} \
        --output {params.kd_dir} \
        --output-prefix kd\
        --reference-db {params.human_db} \
        --threads {threads} --log {log} \
        {params.extra_args}

        # gzip all output fastq files
        gzip {params.kd_dir}/*.fastq
        """


# count reads at each step
rule read_counts:
    input:
        pre_R1=results + "/00_INPUT_DATA/01_reads/{sample}_R1.fastq.gz",
        kd_R1=results
        + "/01_READ_PREPROCESSING/01_kneaddata/{sample}/kd_paired_1.fastq.gz",
        kd_R1S=results
        + "/01_READ_PREPROCESSING/01_kneaddata/{sample}/kd_unmatched_1.fastq.gz",
        kd_R2S=results
        + "/01_READ_PREPROCESSING/01_kneaddata/{sample}/kd_unmatched_2.fastq.gz",
    output:
        read_counts=results
        + "/01_READ_PREPROCESSING/02_read_counts/{sample}_read_counts.csv",
    params:
        unzip_pre_R1=results + "/00_INPUT_DATA/01_reads/{sample}_R1.fastq",
        unzip_kd_R1=results
        + "/01_READ_PREPROCESSING/01_kneaddata/{sample}/kd_paired_1.fastq",
        unzip_kd_R1S=results
        + "/01_READ_PREPROCESSING/01_kneaddata/{sample}/kd_unmatched_1.fastq",
        unzip_kd_R2S=results
        + "/01_READ_PREPROCESSING/01_kneaddata/{sample}/kd_unmatched_2.fastq",
    threads: 1
    shell:
        """
        # unzip read files
        gunzip -f {input.pre_R1}
        gunzip {input.kd_R1}
        gunzip {input.kd_R1S}
        gunzip {input.kd_R2S}

        # count reads in each file
        touch {output.read_counts}
        echo "sample,unprocessed_paired,kneaddata_paired,kneaddata_R1S,kneaddata_R2S" > {output.read_counts}
        echo "{wildcards.sample},\
        $(($(cat {params.unzip_pre_R1}|wc -l)/4)),\
        $(($(cat {params.unzip_kd_R1}|wc -l)/4)),\
        $(($(cat {params.unzip_kd_R1S}|wc -l)/4)),\
        $(($(cat {params.unzip_kd_R2S}|wc -l)/4))" >> {output.read_counts}

        # gzip read files
        gzip {params.unzip_pre_R1}
        gzip {params.unzip_kd_R1}
        gzip {params.unzip_kd_R1S}
        gzip {params.unzip_kd_R2S}
        """


# combine read count files
rule combine_read_counts:
    input:
        expand(
            results + "/01_READ_PREPROCESSING/02_read_counts/{sample}_read_counts.csv",
            sample=samples,
        ),
    output:
        results + "/01_READ_PREPROCESSING/read_preprocessing_report.csv",
    threads: 1
    shell:
        """
        # combine read counts files, and only keep header row from one file
        awk 'FNR>1 || NR==1' {input} > {output}
        """


# visualize read counts
rule read_preprocessing_analysis:
    input:
        results + "/01_READ_PREPROCESSING/read_preprocessing_report.csv",
    output:
        report(
            results + "/01_READ_PREPROCESSING/read_preprocessing_figure.png",
            caption="../report/read_preprocessing_analysis.rst",
            category="Step 01: Read preprocessing",
        ),
    threads: 1
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/read_preprocessing_analysis.py.ipynb"
