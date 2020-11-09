"""Extracts results from the list of input files

Notes
-----
Note that this script is directly implemented into the Snakemake
workflow and pulls the arguments for the ``main()`` function directly from
``snakemake.input`` and ``snakemake.output[0]`` attributes on the snakemake
object passed into this module at run time.
"""
import pandas as pd
from typing import List, Tuple, Dict
from otoole.input import Strategy
import os
import sys

from logging import getLogger

logger = getLogger(__name__)

def add_dtypes(config: Dict):
    for name, details in config.items():
        if details["type"] == "param":
            dtypes = {}
            for column in details["indices"] + ["VALUE"]:
                if column == "VALUE":
                    dtypes["VALUE"] = details["dtype"]
                else:
                    dtypes[column] = config[column]["dtype"]
            details["index_dtypes"] = dtypes
    return config

strategy = Strategy()
config = strategy.results_config.copy()
config.update(strategy.input_config)
logger.debug(config)
config = add_dtypes(config)

def get_model_run_scenario_from_filepath(filename: str):
    """Parses filepath to extract useful bits

    "results/{{scenario}}/{modelrun}/{input_file}.csv"
    """
    filepath, name = os.path.split(filename)
    result_param = os.path.splitext(name)[0]
    scenario_path, model_run = os.path.split(filepath)
    scenario = os.path.split(scenario_path)[1]
    return {'model_run': model_run, 'scenario': scenario,
            'result_param': result_param, 'filepath': filepath}

def main(input_files: List, output_file: str, parameter: Tuple):
    """Iterate over list of CSV files, extract defined results, write to output file
    """
    aggregated_results = []
    for filename in input_files:

        bits = get_model_run_scenario_from_filepath(filename)

        column_dtypes = config[bits['result_param']]['index_dtypes']
        df_index = config[bits['result_param']]['indices']

        df = pd.read_csv(filename
            ).astype(column_dtypes
            ).set_index(df_index)

        try:
            interconnector = df.xs(parameter, drop_level=False)
            interconnector['SCENARIO'] = bits['scenario']
            interconnector['MODELRUN'] = bits['model_run']
            interconnector = interconnector.reset_index(
                ).set_index(['SCENARIO', 'MODELRUN'] + df_index)

            aggregated_results.append(interconnector)
        except KeyError:
            logger.warning("No results found for %s in %s", parameter, bits['filepath'])

    results = pd.concat(aggregated_results)

    results.to_csv(output_file)

# if __name__ == "__main__":
#     main(sys.argv[2:], sys.argv[1])
# else:
input_files : List = snakemake.input
output_file = snakemake.output[0]
parameter = tuple(snakemake.params['parameter'].split(","))
main(input_files, output_file, parameter)