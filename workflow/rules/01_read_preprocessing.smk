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

# load resources folder path
resources = config["resources_path"]

# load report
report: "../report/workflow.rst"

# load sample information to be used in workflow
samples = samples_df["sample"]
samples_assemblies=list(set(samples_df["sample"].astype("string") + "_" + samples_df["assembly"].astype("string")))
R1_files=samples_df["R1"]
R2_files=samples_df["R2"]

# -------------------------------------
# Read Preprocessing Rules
# -------------------------------------
### Set up workflow ###
# create a symlink for each reads file (helps with downstream processing)
rule symlink_reads:
    input:
        R1_display=expand("{R1}", R1=R1_files),
        R2_display=expand("{R2}", R2=R2_files),
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
        R1_display=expand(results + "/00_INPUT_DATA/01_reads/{sample_assembly_replicate}_R1.fastq.gz", sample_assembly_replicate=samples_assemblies_replicates),
        R2_display=expand(results + "/00_INPUT_DATA/01_reads/{sample_assembly_replicate}_R2.fastq.gz", sample_assembly_replicate=samples_assemblies_replicates),
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
        log=results + "/01_READ_PREPROCESSING/02_clumpify/{sample_assembly}.log",
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
        {params.extra_args} > {output.log} 2>&1

        # remove intermediate files
        gzip {params.R1}
        gzip {params.R2}
        """

### Quality filter reads using KneadData ###
# build kneaddata bowtie2 database
rule build_kneaddata_database:
    output:
        resources + "/kneaddata/hg37dec_v0.1.1.bt2",
        resources + "/kneaddata/hg37dec_v0.1.2.bt2",
        resources + "/kneaddata/hg37dec_v0.1.3.bt2",
        resources + "/kneaddata/hg37dec_v0.1.4.bt2",
        resources + "/kneaddata/hg37dec_v0.1.rev.1.bt2",
        resources + "/kneaddata/hg37dec_v0.1.rev.2.bt2",
    params:
        human_db_dir=resources + "/kneaddata/",
    threads: 1
    priority: 1
    conda:
        "../envs/kneaddata.yml"
    shell:
        """
        # download human genome reference to desired directory
        kneaddata_database --download human_genome bowtie2 {params.human_db_dir}
        """

# Determine inputs for kneaddata (changes if clumpify will not be run)
if config['clumpify']['run_clumpify']:
    kneaddata_R1_input=results + "/01_READ_PREPROCESSING/02_clumpify/{sample_assembly}_R1.fastq.gz"
    kneaddata_R2_input=results + "/01_READ_PREPROCESSING/02_clumpify/{sample_assembly}_R2.fastq.gz"
else:
    kneaddata_R1_input=results + "/01_READ_PREPROCESSING/01_merge_replicates/{sample_assembly}_R1.fastq.gz"
    kneaddata_R2_input=results + "/01_READ_PREPROCESSING/01_merge_replicates/{sample_assembly}_R2.fastq.gz"


# quality filter reads with kneaddata
rule kneaddata:
    input:
        resources + "/kneaddata/hg37dec_v0.1.1.bt2",
        resources + "/kneaddata/hg37dec_v0.1.2.bt2",
        resources + "/kneaddata/hg37dec_v0.1.3.bt2",
        resources + "/kneaddata/hg37dec_v0.1.4.bt2",
        resources + "/kneaddata/hg37dec_v0.1.rev.1.bt2",
        resources + "/kneaddata/hg37dec_v0.1.rev.2.bt2",
        R1=kneaddata_R1_input,
        R2=kneaddata_R2_input,
    output:
        R1=results
        + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_paired_1.fastq.gz",
        R2=results
        + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_paired_2.fastq.gz",
        R1S=results
        + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_unmatched_1.fastq.gz",
        R2S=results
        + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_unmatched_2.fastq.gz",
        log=results
        + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}.log",
    params:
        output_dir=directory(results + "/01_READ_PREPROCESSING/03_kneaddata/"),
        human_db=resources + "/kneaddata/",
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
        --output {params.output_dir} \
        --output-prefix {params.prefix} \
        --reference-db {params.human_db} \
        --threads {threads} \
        {params.extra_args}

        # gzip all output fastq files
        gzip {params.output_dir}/*.fastq
        """

### Analyze preprocessing results ###
# determine clumpify read counts
rule clumpify_read_counts:
    input:
        expand(results
        + "/01_READ_PREPROCESSING/02_clumpify/{sample_assembly}.log", sample_assembly=samples_assemblies),
    output:
        results
        + "/01_READ_PREPROCESSING/02_clumpify/combined_read_counts.tsv",
    threads: 1
    priority: 3
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/01_read_preprocessing_clumpify_read_counts.py.ipynb"


# determine read counts using kneaddata utils
rule kneaddata_read_counts:
    input:
        expand(results
        + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}.log", sample_assembly=samples_assemblies),
    output:
        results
        + "/01_READ_PREPROCESSING/03_kneaddata/combined_read_counts.tsv",
    params:
        log_dir=results + "/01_READ_PREPROCESSING/03_kneaddata",
    threads: 1
    priority: 3
    conda:
        "../envs/kneaddata.yml"
    shell:
        """
        # generate read counts from kneaddata log files
        kneaddata_read_count_table \
        --input {params.log_dir} \
        --output {output}
        """

# Determine inputs for read counts analysis
if config['clumpify']['run_clumpify']:
    clumpify_read_count_input=results + "/01_READ_PREPROCESSING/02_clumpify/combined_read_counts.tsv",
else:
    clumpify_read_count_input=results + "/01_READ_PREPROCESSING/03_kneaddata/combined_read_counts.tsv",

### visualize read counts ###
rule read_count_analysis:
    input:
        clumpify=clumpify_read_count_input,
        kneaddata=results + "/01_READ_PREPROCESSING/03_kneaddata/combined_read_counts.tsv",
    output:
        figure=report(
            results + "/01_READ_PREPROCESSING/read_count_figure.png",
            caption="../report/01_read_preprocessing_read_count_analysis.rst",
            category="Step 01: Read preprocessing",
        ),
        report=results + "/01_READ_PREPROCESSING/read_count_report.tsv",
    params:
        run_clumpify=config['clumpify']['run_clumpify'],
    threads: 1
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/01_read_preprocessing_read_count_analysis.py.ipynb"
