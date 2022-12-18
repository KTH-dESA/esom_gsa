"""Plots interactions per trajectory. 

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
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from math import ceil
import matplotlib.ticker as ticker
import sys
from pathlib import Path

from logging import getLogger

logger = getLogger(__name__)

def main(parameters: pd.DataFrame, X: np.array, results: pd.DataFrame, save_file: str):
    """Plots interactions. 

    Parameters
    ----------
    parameters : pd.DataFrame
        Parameters for generated sample
    X : np.array
        Input Sample
    model_results : pd.DataFrame
        Model results for a OSeMOSYS variable 
    save_file : str
        File path to save results
    """

    # Save sample data as dataframe 
    X = pd.DataFrame(X, columns=parameters['group'].to_list())
    X = X.loc[:,~X.columns.duplicated()]

    # Add objective cost column 
    Y_max = results['OBJECTIVE'].max()
    Y = results['OBJECTIVE'] / Y_max
    X['Objective_Cost'] = Y

    # save individual model results 
    models = {}
    for row in range(len(X)):
        models[f'model_{row}'] = X.loc[row, :].to_list()
    categories = list(X)

    # plot 
    num_plots = int(len(X) / len(X.columns))
    num_rows = ceil(num_plots / 2)
    runs_per_trajectory = len(X.columns) # + 1 taken care of by adding in result column

    fig, axs = plt.subplots(nrows=num_rows, ncols=2, sharex=True, sharey=True, figsize=(10,10), gridspec_kw = {'wspace':0.075, 'hspace':0.2})

    counter = 0
    row = 0
    col = 0
    secondary_ax = []
    for model_num, model_data in models.items():
        axs[row, col].plot(categories, model_data, label = model_num)
        axs[row, col].yaxis.set_major_locator(ticker.MultipleLocator(0.3333)) # 4 levels in morris
        axs[row, col].yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
        counter += 1
        if (counter) == runs_per_trajectory: 
            secondary_ax.append(axs[row, col].twinx())
            counter = 0
            col = (col + 1) % 2
            row = row if col == 1 else row + 1
    axs[num_rows - 1, 0].tick_params('x',labelrotation=90)
    axs[num_rows - 1, 1].tick_params('x',labelrotation=90)

    row = -1 # start at -1 to put counter at start of loop
    for index, _ in enumerate(secondary_ax):
        row += 1
        if row % 2 != 0:
            continue
        secondary_ax[index].get_shared_y_axes().join(secondary_ax[index], secondary_ax[index + 1])
        secondary_ax[index].plot()
        secondary_ax[index].yaxis.set_tick_params(labelright=False)
        secondary_ax[index].set_ylim(0, Y_max)

    fig.supylabel('Sample Value', fontsize='x-large', x=0.05)
    fig.suptitle('Interactions on Objective Cost', fontsize='x-large', y=0.93)
    fig.savefig(f'{save_file}.png', bbox_inches='tight')


if __name__ == "__main__":

    parameters_file = sys.argv[1]
    sample = sys.argv[2]
    result_file = sys.argv[3]
    save_file = str(Path(sys.argv[4]).with_suffix(''))
    parameters = pd.read_csv(parameters_file)

    X = np.loadtxt(sample, delimiter=',')
    results = pd.read_csv(result_file)

    main(parameters, X, results, save_file)
