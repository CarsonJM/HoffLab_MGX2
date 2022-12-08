# import packages
import plotly.express as px
import pandas as pd

# read kneaddata read counts
kd = pd.read_csv(str(snakemake.input), sep='\t')

# rename columns and prepare for plotting
kd.rename(columns = {'raw pair1': 'Before human read removal 1', 'final pair1': 'After human read removal 1'}, inplace=True)
kd['Before human read removal'] = kd['Before human read removal 1'] * 2
kd['After human read removal'] = kd['After human read removal 1'] * 2
kd_melt = kd.melt(id_vars=['Sample'], value_vars=['Before human read removal', 'After human read removal'])

# plot kneaddata read counts and save
fig = px.box(kd_melt, x='variable', y="value", points="all", color="variable", hover_name='Sample',
             labels={
                     "variable": "Preprocessing step",
                     "value": "Number of reads"     
                 })
fig.update_layout(title_text='Removal of human contamination')
fig.update_layout(showlegend=False)
fig.write_image(str(snakemake.output))