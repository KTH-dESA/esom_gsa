"""Generates a sample from a list of parameters

Arguments
---------
replicates : int
    The number of model runs to generate
path_to_parameters : str
    File containing the parameters to generate a sample for

Usage
-----
To run the script on the command line, type::

    python create_sample.py 10 path/to/parameters.csv

The ``parameters.csv`` CSV file should be formatted as follows::

    name,group,indexes,min_value,max_value,dist,interpolation_index,action
    CapitalCost,pvcapex,"GLOBAL,GCPSOUT0N",500,1900,unif,YEAR,interpolate
    DiscountRate,discountrate,"GLOBAL,GCIELEX0N",0.05,0.20,unif,None,fixed

"""
from SALib.sample import latin
import os
import csv
from typing import List
import sys

from logging import getLogger

logger = getLogger(__name__)

def main(parameters: List, output_files: List):

    replicates = len(output_files)

    problem = {}
    problem['num_vars'] = len(parameters)

    names = []
    bounds = []
    groups = []
    for parameter in parameters:
        names.append(parameter['name'] + ";" + parameter['indexes'])

        min_value = float(parameter['min_value'])
        max_value = float(parameter['max_value'])

        groups.append(parameter['group'])

        bounds.append([min_value, max_value])

    problem['names'] = names
    problem['bounds'] = bounds
    problem['groups'] = groups

    sample = latin.sample(problem, replicates, seed=42)

    for model_run, row in enumerate(sample):
        filepath = output_files[model_run]
        with open(filepath, 'w') as csvfile:

            fieldnames = ['name', 'indexes', 'value', 'action', 'interpolation_index']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()

            for column, param in zip(row, parameters):
                data = {'name': param['name'],
                        'indexes': param['indexes'],
                        'value': column,
                        'action': param['action'],
                        'interpolation_index': param['interpolation_index']}
                writer.writerow(data)


if __name__ == "__main__":

    output_files = sys.argv[2:]
    parameters_file = sys.argv[1]
    with open(parameters_file, 'r') as csv_file:
        reader = list(csv.DictReader(csv_file))
    main(reader, output_files)
