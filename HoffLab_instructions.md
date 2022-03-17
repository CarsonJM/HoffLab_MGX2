# HoffLab Metagenomics Pipeline Instructions

## 1. Start a screen session from the head node

`screen -S HoffLab_MGX`

- to detach from the screen use: `Ctrl+A+D`
- to reattach to the screen use: `screen -r HoffLab_MGX`
- to end the screen use when not on it: `screen -X -S HoffLab_MGX kill`

## 2. Log into a compute node

`qrsh -q new.q`

## 3. Make and activate snakemake conda environment using following command

`wget https://github.com/CarsonJM/HoffLab_MGX/blob/master/environment.yml`

`conda env create -f environment.yml`

This may not work, so another option is:

`conda create -y -n HoffLab_MGX -c conda-forge -c bioconda mamba=0.17.0 snakemake=7.2.1 snakefmt=0.4.4 snakedeploy=0.3.0 git=2.34.1`

After creating the environment, run:

`conda activate HoffLab_MGX`

## 4. Make a new directory where you want the analysis to take place

`mkdir <insert directory name here>`

## 5. Change to the location of the new directory

`cd <insert directory name here>`

## 6. Run the following command to deploy the HoffLab_MGX workflow in the specified directory

`snakedeploy deploy-workflow https://github.com/CarsonJM/HoffLab_MGX.git . --tag v0.1-beta`

*You should see 'config' and 'workflow' directories now*

## 7. Modify the config/config.yaml file so that:

- the location of the samples_list.tsv matches yours
- the resource directory is where you want to store large downloaded databases
- other workflow parameters are set as desired

## 7. Modify the sample.tsv file so that it matches the desired sample names and paths

- the "sample" column should reflect the patient & timepoint of the sample
- the "replicate" column is included so that multiple runs from the same sample can be merged
(i.e. "sample" 1 "replicate" 1 will be merged with "sample" 1 "replicate" 2)
- the "assembly" column should reflect samples that you want to be co-assembled together
(i.e. "sample" 1 "assembly" 1 will be co-assembled with "sample" 2 "assembly" 1)

## 8. Run the following command to generate a pdf preview of the workflow

`snakemake --dag | dot -Tpdf > dag.pdf`

## 9. Run a dry run of the workflow to verify everything is set up correctly

`snakemake --dry-run`

## 10. Set up a cruncher profile to use the qsub system

*if you want to run the workflow on cruncher with qsub submissions, copy and paste the following text into a file
at /home/<your_username>/.config/snakemake/cruncher/config.yaml*

`cluster: "qsub -V -cwd -q new.q -j y -o config/cluster_logs/{rulename}.{jobid}.log"`

`jobs: 20`

`latency-wait: 60`

`use-conda: True`
  
## 11. Then run the following command to run the workflow:

`snakemake --profile cruncher`

## 12. Then run the following to generate a report from the workflow:

`snakemake --profile cruncher --report <path to desired report location>`
