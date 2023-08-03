"""Creates a time series heat map from a model run. 

Arguments
---------
path_to_parameters : str
    File containing the parameters for generated sample
model_inputs : str
    File path to sample model inputs
model_outputs : str
    File path to model outputs
location_to_save : str
    File path to save results

Usage
-----
To run the script on the command line, type::

    python analyze_results.py path/to/parameters.csv path/to/inputs.txt 
        path/to/model/results.csv path/to/save/SA/results.csv

The ``parameters.csv`` CSV file should be formatted as follows::

    name,group,indexes,min_value,max_value,dist,interpolation_index,action
    CapitalCost,pvcapex,"GLOBAL,GCPSOUT0N",500,1900,unif,YEAR,interpolate
    DiscountRate,discountrate,"GLOBAL,GCIELEX0N",0.05,0.20,unif,None,fixed

The ``inputs.txt`` should be the output from SALib.sample.morris.sample

The ``model/results.csv`` must have an 'OBJECTIVE' column holding results OR
be a formatted output of an OSeMOSYS parameter 

"""

from pathlib import Path
from SALib.analyze import morris as analyze_morris
from SALib.plotting import morris as plot_morris
import numpy as np
import pandas as pd
import csv
import sys
import utils
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List

from logging import getLogger

logger = getLogger(__name__)

def add_zeros(df : pd.DataFrame, model_runs: List[str]) -> pd.DataFrame:
    """Add zeros to perserve shape of dataframes"""
    
    actual_runs = df["MODELRUN"].unique()
    data = []
    
    for run in model_runs:
        if run not in actual_runs:
            data.append([run, 0])
    
    df_to_add = pd.DataFrame(data, columns=["MODELRUN", "VALUE"])
    
    return pd.concat([df, df_to_add]).sort_values("MODELRUN").reset_index(drop=True)

def sort_results(df : pd.DataFrame, year : int) -> np.array:
    """Organizes a model variable results file for a morris analysis
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe grouped on model number and year
    year : int
        Year to sort results for 

    Returns
    -------
    Y : np.array
        Results for morris analysis
    """
    model_runs = df.index.get_level_values("MODELRUN").unique()
    try:
        df = df.xs(year, level=('YEAR')).reset_index()
        df = add_zeros(df, model_runs)
    except KeyError: 
        data = [[x,0] for x in model_runs]
        df = pd.DataFrame(data, columns=["MODELRUN", "VALUE"])
    df['NUM'] = df['MODELRUN'].map(lambda x: int(x.split('_')[1]))
    df = df.sort_values(by='NUM').reset_index(drop=True).drop(('NUM'), axis=1).set_index('MODELRUN')
    Y = df.to_numpy(dtype=float)
    return Y

def main(parameters: dict, X: np.array, model_results: pd.DataFrame, save_file: str, scaled: bool = False):
    """Performs SA and plots results. 

    Parameters
    ----------
    parameters : Dict
        Parameters for generated sample
    X : np.array
        Input Sample
    model_results : pd.DataFrame
        Model results for a OSeMOSYS variable 
    save_file : str
        File path to save results
    scaled : bool = False
        If the input sample is scaled
    """

    problem = utils.create_salib_problem(parameters)
    model_results = model_results.groupby(['MODELRUN','YEAR']).sum(numeric_only=True)

    years = model_results.index.unique(level='YEAR')
    SA_result_data = []

    for year in years:
        Y = sort_results(model_results, year)
        Si = analyze_morris.analyze(problem, X, Y, print_to_console=False, scaled=scaled)
        SA_result_data.append(Si['mu_star'])
    
    SA_result_data = np.ma.concatenate([SA_result_data])
    columns = list(set([x['group'] for x in parameters]))
    SA_results = pd.DataFrame(np.ma.getdata(SA_result_data), columns=columns, index=years).T

    # Save figure results
    title = Path(save_file).stem.capitalize()
    height = len(columns) / 2 + 1.5
    width = len(years) / 5
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(SA_results, cmap="coolwarm", ax=ax).set_title(title)
    ax.set_yticklabels(ax.get_yticklabels(), rotation = 0, fontsize = 8)
    fig.savefig(f'{save_file}.png', bbox_inches='tight')

if __name__ == "__main__":

    parameters_file = sys.argv[1]
    sample = sys.argv[2]
    result_file = sys.argv[3]
    save_file = str(Path(sys.argv[4]).with_suffix(''))
    scaled = sys.argv[5]
    scaled = False if scaled == "False" else True
    
    with open(parameters_file, 'r') as csv_file:
        parameters = list(csv.DictReader(csv_file))

    X = np.loadtxt(sample, delimiter=',')
    results = pd.read_csv(result_file)

    main(parameters, X, results, save_file, scaled)
    
