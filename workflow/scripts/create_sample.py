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
from SALib.sample import morris
import os
import numpy as np
import csv
from typing import List
import sys

from logging import getLogger

logger = getLogger(__name__)

def main(parameters: List, sample_file: str, replicates: int):

    problem = {}
    problem['num_vars'] = len(parameters)

    names = []
    bounds = []
    groups = []
    for parameter in parameters:
        names.append(parameter['name'] + ";" + parameter['indexes'])
        groups.append(parameter['group'])
        min_value = 0
        max_value = 1
        bounds.append([min_value, max_value])

    problem['names'] = names
    problem['bounds'] = bounds
    problem['groups'] = groups

    sample = morris.sample(problem, replicates, seed=42)
    np.savetxt(sample_file, sample, delimiter=',')


if __name__ == "__main__":

    parameters_file = sys.argv[1]
    sample_file = sys.argv[2]
    replicates = int(sys.argv[3])
    with open(parameters_file, 'r') as csv_file:
        reader = list(csv.DictReader(csv_file))
    main(reader, sample_file, replicates)
