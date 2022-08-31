# Metagenomics pipeline instructions

## 1. Start a screen session from the head node

`screen -S phide_piper`

- to detach from the screen use: `Ctrl+A+D`
- to reattach to the screen use: `screen -r phide_piper`
- to end the screen use when not on it: `screen -X -S phide_piper kill`

## 2. Log into a compute node

`qrsh -q new.q`

## 3. Make and activate snakemake conda environment using following command

`conda create -n phide_piper`

`conda activate phide_piper`

`conda install -c conda-forge -c bioconda mamba snakemake snakefmt snakedeploy git -y`

## 4. Make a new directory where you want the analysis to take place

`mkdir <insert directory name here>`

## 5. Change to the location of the cloned directory

`cd <insert directory name here>`

## 6. Run the following command to deploy the phide_piper workflow in the specified directory

`snakedeploy deploy-workflow https://github.com/CarsonJM/phide_piper_dev.git . --branch master`

*You should see a 'config' and 'workflow' directories now*

## 7. Modify the config/config.yaml file so that:

- the resource directory is where you want to store large downloaded databases
- other workflow parameters are set as desired

## 7. Modify the sample.tsv file so that it matches the desired sample names and paths

## 8. Run the following command to generate a pdf preview of the workflow

`snakemake --configfile <path to config file> --dag | dot -Tpdf > dag.pdf`

## 9. Run a dry run of the workflow to verify everything is set up correctly

`snakemake --configfile <path to config file> --dry-run`

## 10. Set up a cruncher profile to use the qsub system

*if you want to run the workflow on cruncher with qsub submissions, copy and paste the following text into a file
at /home/<your_username>/.config/snakemake/cruncher/config.yaml*

`cluster: "qsub -V -cwd -q new.q -j y -o config/cluster_logs/{rulename}.{jobid}.log"
jobs: 20
latency-wait: 60
use-conda: True`
  
## 11. Then run the following command to run the workflow:

`snakemake --profile cruncher --configfile <path to config file>`

## 12. Then run the following to generate a report from the workflow:

`snakemake --profile cruncher --configfile <path to config file> --report <path to desired report location>`