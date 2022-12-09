import pandas as pd

# read metaphlan data
metaphlan_data = pd.read_csv(str(snakemake.input), header=1, sep='\t')

# format metaphlan data for plotting kingdom
def format_taxon_level(metaphlan_profile, taxonomic_level, output, number, delimiter):
    level_filtered = metaphlan_profile[metaphlan_profile['clade_name'].str.count('\|') == number]
    level_filtered[taxonomic_level] = level_filtered['clade_name'].str.split(delimiter, expand=True)[1]
    column = level_filtered[taxonomic_level]
    level_filtered.drop(taxonomic_level, inplace=True, axis=1)
    level_filtered.insert(loc=0, column=taxonomic_level, value=column)
    level_filtered.drop('clade_name', inplace=True, axis=1)
    level_filtered.to_csv(output, sep='\t', index=False)

metaphlan_phylum = format_taxon_level(metaphlan_data, 'phylum', str(snakemake.output.phylum), 1, 'p__')
metaphlan_genus = format_taxon_level(metaphlan_data, 'genus', str(snakemake.output.genus), 5, 'g__')
metaphlan_species = format_taxon_level(metaphlan_data, 'species', str(snakemake.output.species), 6, 's__')