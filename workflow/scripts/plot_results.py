import pandas as pd
import seaborn as sns
import sys
import os

sns.set_context("paper")

results_filepath = sys.argv[1]
plot_filepath = sys.argv[2]

nice_name = os.path.split(plot_filepath)[1]
nice_name = nice_name.split(".")[0]  # Remove file ending
nicer_name = " ".join([x.capitalize() for x in nice_name.split('_')])

# Read in the data
data = pd.read_csv(results_filepath)


plot = sns.relplot(x="YEAR", y='VALUE', ci=None, kind="line", 
                   hue='TECHNOLOGY', data=data)
# plot.set_xticklabels(data.loc[:,'YEAR'].unique(), rotation=90)
plot.savefig(plot_filepath)
