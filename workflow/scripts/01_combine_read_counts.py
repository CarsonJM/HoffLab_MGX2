import pandas as pd

# read kneaddata read count data
kneaddata_read_counts = pd.read_csv(str(snakemake.input.kneaddata), sep='\t')

# incorporate clumpify read counts
clumpify_read_counts = pd.read_csv(str(snakemake.input.clumpify), sep='\t')
read_counts = clumpify_read_counts.merge(kneaddata_read_counts, on='Sample', how='right')
read_counts.rename(columns={'raw pair1': 'deduplicated pair1', 'raw pair2': 'deduplicated pair2', 'raw reads':'raw pair1'}, inplace=True)

# save combined read counts file
read_counts.to_csv(str(snakemake.output), index=False, sep='\t')