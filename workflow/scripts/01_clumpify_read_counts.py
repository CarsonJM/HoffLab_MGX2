import pandas as pd

# create list to store information
sample = []
raw_reads = []

# iterate through each file in snakemake input
for log in snakemake.input:
    # read from log file
    read_log = open(str(log), 'r')
    lines = read_log.readlines()

    # save file name as sample
    relative_path = log.rpartition('/')[2]
    sample.append(relative_path.split('.log')[0])

    # save read count information
    for line in lines:
        if 'Reads In:' in line:
            line_strip = line.strip()
            reads_in_count = line_strip.split()[2]
            raw_reads.append(int(reads_in_count)/2)
            continue

# create pandas df of data and save as tsv
clumpify_read_counts_df = pd.DataFrame(data={'Sample':sample, 'raw reads': raw_reads})
clumpify_read_counts_df.to_csv(str(snakemake.output), sep='\t', index=False)