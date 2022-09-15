import pathlib
import pandas as pd
import seaborn as sns
import sys
from pathlib import Path

def main(data:pd.DataFrame, plot_name:pathlib.Path):

    sns.set_context("paper")

    plot = sns.relplot(x="YEAR", y='VALUE', kind="line", hue='TECHNOLOGY', data=data)

    plot.savefig(plot_name)

if __name__ == '__main__':

    results_filepath = sys.argv[1]
    save_directory = sys.argv[2]

    data = pd.read_csv(results_filepath)
    data = data[~data['TECHNOLOGY'].str.startswith('MINE_')]
    nice_name = Path(results_filepath).stem.capitalize()
    plot_name = Path(save_directory, nice_name)

    main(data, plot_name)
