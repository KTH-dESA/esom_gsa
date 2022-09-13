"""Analyzes results from model

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

The ``model/results.csv`` must have an 'OBJECTIVE' column holding results

"""

from math import ceil
from SALib.analyze import morris as analyze_morris
from SALib.plotting import morris as plot_morris
import os
import numpy as np
import pandas as pd
import csv
from typing import List
import sys
import utils
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from logging import getLogger

logger = getLogger(__name__)

def plot_histogram(problem: dict, X: np.array, fig: plt.figure):

    # chnage histogram labels to legend for clarity
    problem_hist = problem.copy()
    problem_hist['names'] = [f'X{x}' for x, _ in enumerate(problem_hist['names'])]
    legend_labels = [f"{problem_hist['names'][num]} = {problem['names'][num]}" for num, _ in enumerate(problem['names'])]
    legend_handles = [mlines.Line2D([],[], color='w', marker='.', linewidth=0, markersize=0, label=label) for label in legend_labels]

    # plot histogram 
    ax = fig.subplots(1)
    plot_morris.sample_histograms(fig, X, problem_hist)
    fig.patch.set_visible(False)
    ax.axis('off')
    ncols = 2 if len(legend_labels) < 3 else ceil(len(legend_labels)/2) 
    fig.legend(handles=legend_handles, ncol=ncols, frameon=False)
    fig.suptitle(' ', fontsize=(ncols * 20))


def main(parameters: dict, X: np.array, Y: np.array, save_file: str):

    problem = utils.create_salib_problem(parameters)

    Si = analyze_morris.analyze(problem, X, Y, print_to_console=True)

    # Save text based results
    Si.to_df().to_csv(f'{save_file}.csv')
    
    # save graphical resutls 
    fig = plt.figure(figsize=(16, 8), constrained_layout=True)
    subfigs = fig.subfigures(1, 2)
    plot_histogram(problem, X, subfigs[0])
    axs_right = subfigs[1].subplots(2)
    plot_morris.horizontal_bar_plot(axs_right[0], Si, unit="(\$)")
    plot_morris.covariance_plot(axs_right[1], Si, unit="(\$)")

    fig.savefig(f'{save_file}.png')

if __name__ == "__main__":

    parameters_file = sys.argv[1]
    sample = sys.argv[2]
    model_results = sys.argv[3]
    save_file = sys.argv[4][:-4] #remove '.csv' extension

    with open(parameters_file, 'r') as csv_file:
        parameters = list(csv.DictReader(csv_file))
    
    X = np.loadtxt(sample, delimiter=',')
    Y = pd.read_csv(model_results)['OBJECTIVE'].to_numpy()
    
    main(parameters, X, Y, save_file)
