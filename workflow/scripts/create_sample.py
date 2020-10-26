from SALib.sample import latin
import os
import csv
from typing import List
import sys

PARAMETERS = os.path.join('config', 'parameters.csv')

def main(parameters: List, replicates: int):

    problem = {}
    problem['num_vars'] = len(parameters)

    names = []
    bounds = []
    for parameter in parameters:
        names.append(parameter['name'] + ";" + parameter['indexes'])
        
        min_value = float(parameter['min_value'])
        max_value = float(parameter['max_value'])
        
        bounds.append([min_value, max_value])

    problem['names'] = names
    problem['bounds'] = bounds

    sample = latin.sample(problem, replicates)

    for model_run, row in enumerate(sample):
        filename = "{}_sample.txt".format(model_run)
        filepath = os.path.join('modelruns', filename)
        with open(filepath, 'w') as csvfile:

            fieldnames = ['name', 'indexes', 'value']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()

            for column, param in zip(row, parameters):
                data = {'name': param['name'],
                        'indexes': param['indexes'],
                        'value': column}
                writer.writerow(data)


if __name__ == "__main__":

    replicates = sys.argv[1]
    with open(PARAMETERS, 'r') as csv_file:
        reader = list(csv.DictReader(csv_file))
    main(reader, int(replicates))