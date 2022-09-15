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
# reconstruct all species strains in each sample
rule sample2markers:
    input:
        results + "02_READ_BASED_TAXONOMY/01_metaphlan/{sample}.sam.bz2"
    output:
        resources + "metaphlan/mpa_v30_CHOCOPhlAn_201901.1.bt2"
    log:
        results + "00_LOGS/04_sample2markers.{sample}.log"
    params:
        out_dir=results + "03_READ_BASED_PHYLOGENY/01_strainphlan/"
    conda:
        "../envs/humann.yml"
    container: 
        "docker://quay.io/biocontainers/humann:3.0.1--pyh5e36f6f_0"
    benchmark:
        "benchmark/02_READ_BASED_TAXONOMY/download_humann_db.tsv"
    resources:
        runtime="04:00:00",
        partition="ckpt",
        mem_mb="10000",
    threads:
        config["read_function"]["humann_threads"]
    shell:
        """
        # combine paired end reads
        cat {input.R1} {input.R2} > {params.combined}

        # run humann on samples
        humann --input {params.combined} \
        --output {params.out_dir} \
        --threads {threads} \
        --taxonomic-profile {input.mpa} \
        --input_format fastq \
        --o-log {log} \
        --nucleotide-database {params.cpa_dir} \
        --protein-database {params.uniref_dir} \
        --pathways {params.map_dir} \
        --remove-temp-output \
        {params.extra_args}
        """
