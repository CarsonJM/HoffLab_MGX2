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
    output:
        results + "04_READ_BASED_STRAIN/01_consensus_markers/{sample}.pkl",
    params:
        out_dir=results + "04_READ_BASED_STRAIN/01_consensus_markers/",
        extra_args=config["read_strain"]["samples2markers_arguments"],
        mpa=resources + "metaphlan/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl",
    # conda:
    #     "../envs/humann:3.6--pyh7cba7a3_1.yml"
    container:
        "docker://quay.io/biocontainers/humann:3.6--pyh7cba7a3_1"
    benchmark:
        "benchmark/04_READ_BASED_STRAIN/sample2markers_{sample}.tsv"
    resources:
        runtime="08:00:00",
        mem_mb="100000",
    threads: config["read_strain"]["samples2markers_threads"]
    shell:
        """
        # combine paired end reads
        sample2markers.py -i {input} \
        -o {params.out_dir} \
        --nprocs {threads} \
        --input_format bz2 \
        --database {params.mpa} \
        {params.extra_args}
        """


# identify all clades for strainphlan
rule strainphlan_clades:
    input:
        expand(
            results + "04_READ_BASED_STRAIN/01_consensus_markers/{sample}.pkl",
            sample=samples,
        ),
    output:
        results + "04_READ_BASED_STRAIN/02_clades/strainphlan_clades",
    params:
        out_dir=results + "04_READ_BASED_STRAIN/02_clades/",
        spa_db=resources + "metaphlan/mpa_v30_CHOCOPhlAn_201901.pkl",
    # conda:
    #     "../envs/humann:3.6--pyh7cba7a3_1.yml"
    container:
        "docker://quay.io/biocontainers/humann:3.6--pyh7cba7a3_1"
    benchmark:
        "benchmark/04_READ_BASED_STRAIN/strainphlan_clades.tsv"
    resources:
        runtime="01:00:00",
        mem_mb="1000",
    shell:
        """
        # combine paired end reads
        strainphlan \
        --samples {input} \
        --database {params.spa_db} \
        --output_dir {params.out_dir} \
        --print_clades_only > {output}
        """


# # identify all clades for strainphlan
# checkpoint order_clades_by_relab:
#     input:
#         results + "04_READ_BASED_PHYLOGENY/02_clades/strainphlan_clades"
#     output:
#         results + "04_READ_BASED_PHYLOGENY/02_clades/strainphlan_clades_ordered_by_relab.tsv"
#     log:
#         results + "00_LOGS/04_strainphlan_order_clades_by_relab.log"
#     params:
#         out_dir=results + "04_READ_BASED_PHYLOGENY/02_clades/",
#         spa_db=resources + "metaphlan/mpa_v30_CHOCOPhlAn_201901.pkl",
#     conda:
#         "../envs/biobakery_workflows.yml"
#     container:
#         "docker://biobakery/workflows:3.0.0.a.7"
#     benchmark:
#         "benchmark/04_READ_BASED_PHYLOGENY/strainphlan_clades.tsv"
#     resources:
#         runtime="01:00:00",
#         mem_mb="1000",
#     notebook:
# # extract markers for each clade
# rule extract_markers:
#     output:
#         results + "04_READ_BASED_PHYLOGENY/02_db_markers/s__Bacteroides_caccae.fna"
#     log:
#         results + "00_LOGS/04_extract_markers.s__Bacteroides_caccae.log"
#     params:
#         out_dir=results + "04_READ_BASED_PHYLOGENY/02_db_markers/",
#         spa_db=resources + "metaphlan/mpa_v30_CHOCOPhlAn_201901.pkl",
#     conda:
#         "../envs/humann.yml"
#     container:
#         "docker://quay.io/biocontainers/humann:3.0.1--pyh5e36f6f_0"
#     benchmark:
#         "benchmark/04_READ_BASED_PHYLOGENY/extract_markers_s__Bacteroides_caccae.tsv"
#     resources:
#         runtime="00:10:00",
#         partition="ckpt",
#         mem_mb="1000",
#     shell:
#         """
#         # extract markers
#         extract_markers.py --database {params.spa_db} -c s__Bacteroides_caccae -o {params.out_dir}
#         """
# # def determine_clades(wildcards):
# #     clades = []
# #     with checkpoints.strainphlan_clades.get(**wildcards).output[0].open() as file:
# #         for line in file:
# #             if "s__" in line:
# #                 clades.add(line.strip().split("\t")[1].split(": in ")[0])
# #     return expand(results + "04_READ_BASED_PHYLOGENY/04_strainphlan/{clade}.info", clade=clades)
# # run strainphlan
# rule strainphlan:
#     input:
#         pkl=expand(results + "04_READ_BASED_PHYLOGENY/01_consensus_markers/{sample}.pkl", sample=samples),
#         db_markers=results + "04_READ_BASED_PHYLOGENY/02_db_markers/s__Bacteroides_caccae.fna"
#     output:
#         results + "04_READ_BASED_PHYLOGENY/04_strainphlan/s__Bacteroides_caccae.info"
#     log:
#         results + "00_LOGS/04_strainphlan.s__Bacteroides_caccae.log"
#     params:
#         out_dir=results + "04_READ_BASED_PHYLOGENY/04_strainphlan/",
#         spa_db=resources + "metaphlan/mpa_v30_CHOCOPhlAn_201901.pkl",
#     conda:
#         "../envs/humann.yml"
#     container:
#         "docker://quay.io/biocontainers/humann:3.0.1--pyh5e36f6f_0"
#     benchmark:
#         "benchmark/04_READ_BASED_PHYLOGENY/strainphlan_s__Bacteroides_caccae.tsv"
#     resources:
#         runtime="04:00:00",
#         partition="ckpt",
#         mem_mb="10000",
#     threads:
#         config["read_phylogeny"]["strainphlan_threads"]
#     shell:
#         """
#         # combine paired end reads
#         strainphlan \
#         --database {params.spa_db} \
#         -s {input.pkl} \
#         -m {input.db_markers} \
#         -o {params.out_dir} \
#         -n {threads} \
#         -c s__Bacteroides_caccae \
#         --mutation_rates
#         """
