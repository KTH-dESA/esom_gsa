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

index = data.columns[0]

# Remove rows with only zero values
data = data[(data.T != 0).any()]

reshaped = data.melt(id_vars=index, var_name='Years', value_name=nicer_name)
plot = sns.relplot(x="Years", y=nicer_name, ci=None, kind="line", 
                   hue=index, data=reshaped)
plot.set_xticklabels(reshaped.loc[:,'Years'].unique(), rotation=90)
plot.savefig(plot_filepath)
