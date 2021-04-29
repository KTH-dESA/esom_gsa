"""Example data processing script

Run using the command

python run_status.py <output_file> <modelruns: int> <scenarios: int>
"""
import sys
from typing import List
from csv import DictReader, DictWriter

import pandas as pd

def read_data(filepath):

    with open(filepath, r) as csvfile:
        reader = DictReader(csvfile)
        data = list(reader.readlines())

    return df

def write_data(df, list_of_filepaths):
    """Write out data to disk
    """
    pass

def process_data(df):
    """
    """


    return processed_df

def read_input_data(paths):
    data = []
    for path in paths:
        data.append(pd.read_csv(path))


STATUS_MAP = {
    1: 'LOADED',
    2: "OPTIMAL",
    3: "INFEASIBLE",
    4: "INF_OR_UNBD",
    5: "UNBOUNDED",
    6: "CUTOFF",
    7: "ITERATION_LIMIT",
    8: "NODE_LIMIT",
    9: "TIME_LIMIT",
    10: "SOLUTION_LIMIT",
    11: "INTERRUPTED",
    12: "NUMERIC",
    13: "SUBOPTIMAL",
    14: "INPROGRESS",
    15: "USER_OBJ_LIMIT",
}

def main(output_file: str, modelruns: int, scenarios: int):

    input_data_paths = ["modelruns/{}/{}_sample.txt".format(x, y) for x in scenarios for y in modelruns]

    asset = []

    for x in range(1, int(scenarios) + 1):
        for y in range(0, int(modelruns)):
            path = "results/{}/{}.json".format(x, y)
            data = pd.read_json(path).reset_index()
            status = data[data['index'] == 'Status']['SolutionInfo'].values[0]

            status = STATUS_MAP[int(status)]

            # input_path = "modelruns/{}/{}_sample.txt".format(x, y)
            inputs = {}
            inputs['scenario'] = x
            inputs['modelrun'] = y
            inputs['status'] = status

            asset.append(inputs)

    with open(output_file, 'w') as csv_file:
        headers = ['scenario', 'modelrun', 'status']
        writer = DictWriter(csv_file, headers)
        writer.writeheader()
        writer.writerows(asset)

if __name__ == "__main__":

    main(sys.argv[1], sys.argv[2], sys.argv[3])
