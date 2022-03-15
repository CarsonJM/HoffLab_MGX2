# -------------------------------------
# Read Preprocessing Module
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


# -------------------------------------
# Read Preprocessing Rules
# -------------------------------------
### Set up workflow ###
# create a symlink for each reads file (helps with downstream processing)
rule symlink_reads:
    input:
        R1=lambda wildcards: samples_df[(samples_df["sample"].astype("string") + "_" + samples_df["assembly"].astype("string") + "_" + samples_df["replicate"].astype("string")) == wildcards.sample_assembly_replicate].iloc[
            0
        ]["R1"],
        R2=lambda wildcards: samples_df[(samples_df["sample"].astype("string") + "_" + samples_df["assembly"].astype("string") + "_" + samples_df["replicate"].astype("string")) == wildcards.sample_assembly_replicate].iloc[
            0
        ]["R2"],
    output:
        R1=results + "/00_INPUT_DATA/01_reads/{sample_assembly_replicate}_R1.fastq.gz",
        R2=results + "/00_INPUT_DATA/01_reads/{sample_assembly_replicate}_R2.fastq.gz",
    priority: 3
    shell:
        """
        # create a symlink for each file
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """

# merge replicates from each sample
rule merge_replicates:
    input:
        R1=lambda wildcards: expand(results + "/00_INPUT_DATA/01_reads/{{sample_assembly}}_{replicate}_R1.fastq.gz", replicate=samples_df[samples_df["sample"].astype("string") + "_" + samples_df["assembly"].astype("string") == wildcards.sample_assembly]["replicate"]),
        R2=lambda wildcards: expand(results + "/00_INPUT_DATA/01_reads/{{sample_assembly}}_{replicate}_R2.fastq.gz", replicate=samples_df[samples_df["sample"].astype("string") + "_" + samples_df["assembly"].astype("string") == wildcards.sample_assembly]["replicate"]),
    output:
        R1=results + "/01_READ_PREPROCESSING/01_merge_replicates/{sample_assembly}_R1.fastq.gz",
        R2=results + "/01_READ_PREPROCESSING/01_merge_replicates/{sample_assembly}_R2.fastq.gz",
    threads: 1
    priority: 3
    conda:
        "../envs/clumpify.yml"
    shell:
        """
        # merge replicates for each sample and unzip
        cat {input.R1} > {output.R1}
        cat {input.R2} > {output.R2}
        """

### Deduplicate reads using clumpify ###
rule clumpify:
    input:
        R1=results + "/01_READ_PREPROCESSING/01_merge_replicates/{sample_assembly}_R1.fastq.gz",
        R2=results + "/01_READ_PREPROCESSING/01_merge_replicates/{sample_assembly}_R2.fastq.gz",
    output:
        R1=results + "/01_READ_PREPROCESSING/02_clumpify/{sample_assembly}_R1.fastq.gz",
        R2=results + "/01_READ_PREPROCESSING/02_clumpify/{sample_assembly}_R2.fastq.gz",
    params:
        extra_args=config['clumpify']['extra_args'],
        R1=results + "/01_READ_PREPROCESSING/01_merge_replicates/{sample_assembly}_R1.fastq",
        R2=results + "/01_READ_PREPROCESSING/01_merge_replicates/{sample_assembly}_R2.fastq",
    threads: 1
    priority: 3
    conda:
        "../envs/clumpify.yml"
    shell:
        """
        # merge replicates for each sample and unzip
        gunzip {input.R1}
        gunzip {input.R2}

        # run clumpify
        clumpify.sh \
        in={params.R1} \
        in2={params.R2} \
        out={output.R1} \
        out2={output.R2} \
        dedupe=t \
        optical=t spany=t adjacent=t \
        {params.extra_args}

        # remove intermediate files
        gzip {params.R1}
        gzip {params.R2}
        """

### Quality filter reads using KneadData ###
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
    priority: 1
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
        R1=results + "/01_READ_PREPROCESSING/02_clumpify/{sample_assembly}_R1.fastq.gz",
        R2=results + "/01_READ_PREPROCESSING/02_clumpify/{sample_assembly}_R2.fastq.gz",
    output:
        R1=results
        + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_paired_1.fastq.gz",
        R2=results
        + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_paired_2.fastq.gz",
        R1S=results
        + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_unmatched_1.fastq.gz",
        R2S=results
        + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_unmatched_2.fastq.gz",
    params:
        kd_dir=directory(results + "/01_READ_PREPROCESSING/03_kneaddata/"),
        human_db=resources + "kneaddata/",
        extra_args=config["kneaddata"]["extra_args"],
        prefix="{sample_assembly}",
    threads: 10
    priority: 2
    conda:
        "../envs/kneaddata.yml"
    shell:
        """
        # run kneaddata to quality filter and remove host reads
        kneaddata --input {input.R1} --input {input.R2} \
        --output {params.kd_dir} \
        --output-prefix {params.prefix} \
        --reference-db {params.human_db} \
        --threads {threads} \
        {params.extra_args}

        # gzip all output fastq files
        gzip {params.kd_dir}/*.fastq
        """

### Count reads at each step ###
rule read_counts:
    input:
        pre_R1=results + "/01_READ_PREPROCESSING/01_merge_replicates/{sample_assembly}_R1.fastq.gz",
        clump_R1=results + "/01_READ_PREPROCESSING/02_clumpify/{sample_assembly}_R1.fastq.gz",
        kd_R1=results
        + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_paired_1.fastq.gz",
        kd_R1S=results
        + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_unmatched_1.fastq.gz",
        kd_R2S=results
        + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_unmatched_2.fastq.gz",
    output:
        read_counts=results
        + "/01_READ_PREPROCESSING/04_read_counts/{sample_assembly}_read_counts.csv",
    params:
        unzip_pre_R1=results + "/01_READ_PREPROCESSING/01_merge_replicates/{sample_assembly}_R1.fastq",
        unzip_clump_R1=results + "/01_READ_PREPROCESSING/02_clumpify/{sample_assembly}_R1.fastq",
        unzip_kd_R1=results
        + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_paired_1.fastq",
        unzip_kd_R1S=results
        + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_unmatched_1.fastq",
        unzip_kd_R2S=results
        + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_unmatched_2.fastq",
    threads: 1
    shell:
        """
        # unzip read files
        gunzip -f {input.pre_R1}
        gunzip -f {input.clump_R1}
        gunzip -f {input.kd_R1}
        gunzip -f {input.kd_R1S}
        gunzip -f {input.kd_R2S}

        # count reads in each file
        touch {output.read_counts}
        echo "sample,unprocessed,deduplicated,kneaddata_paired,kneaddata_R1S,kneaddata_R2S" > {output.read_counts}
        echo "{wildcards.sample_assembly},\
        $(($(cat {params.unzip_pre_R1}|wc -l)/4)),\
        $(($(cat {params.unzip_clump_R1}|wc -l)/4)),\
        $(($(cat {params.unzip_kd_R1}|wc -l)/4)),\
        $(($(cat {params.unzip_kd_R1S}|wc -l)/4)),\
        $(($(cat {params.unzip_kd_R2S}|wc -l)/4))" >> {output.read_counts}

        # gzip read files
        gzip -f {params.unzip_pre_R1}
        gzip -f {params.unzip_clump_R1}
        gzip -f {params.unzip_kd_R1}
        gzip -f {params.unzip_kd_R1S}
        gzip -f {params.unzip_kd_R2S}
        """


# # combine read count files
# rule combine_read_counts:
#     input:
#         expand(
#             results + "/01_READ_PREPROCESSING/02_read_counts/{sample}_read_counts.csv",
#             sample=samples,
#         ),
#     output:
#         results + "/01_READ_PREPROCESSING/read_preprocessing_report.csv",
#     threads: 1
#     shell:
#         """
#         # combine read counts files, and only keep header row from one file
#         awk 'FNR>1 || NR==1' {input} > {output}
#         """


# # visualize read counts
# rule read_preprocessing_analysis:
#     input:
#         results + "/01_READ_PREPROCESSING/read_preprocessing_report.csv",
#     output:
#         report(
#             results + "/01_READ_PREPROCESSING/read_preprocessing_figure.png",
#             caption="../report/read_preprocessing_analysis.rst",
#             category="Step 01: Read preprocessing",
#         ),
#     threads: 1
#     conda:
#         "../envs/jupyter.yml"
#     notebook:
#         "../notebooks/read_preprocessing_analysis.py.ipynb"
