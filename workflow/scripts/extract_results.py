"""Extracts results from the list of input files

Notes
-----
Note that this script is directly implemented into the Snakemake
workflow and pulls the arguments for the ``main()`` function directly from
``snakemake.input`` and ``snakemake.output[0]`` attributes on the snakemake
object passed into this module at run time.
"""
import pandas as pd
import pyarrow
from typing import List, Tuple, Dict
from otoole.input import Strategy
from utils import get_model_run_scenario_from_filepath, get_model_run_scenario_from_input_filepath
import os
import sys
from utils import write_results, read_results

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


def main(input_files: List, output_file: str, parameter: Tuple, config: Dict):
    """Iterate over list of CSV files, extract defined results, write to output file
    """
    aggregated_results = []
    for filename in input_files:

        bits = get_model_run_scenario_from_filepath(filename)

        column_dtypes = config[bits['param']]['index_dtypes']
        df_index = config[bits['param']]['indices']

        df = read_results(filename)
        df = df.astype(column_dtypes).set_index(df_index)

        try:
            results = df.xs(parameter, drop_level=False)
            results['SCENARIO'] = bits['scenario']
            results['MODELRUN'] = bits['model_run']
            results = results.reset_index(
                ).set_index(['SCENARIO', 'MODELRUN'] + df_index)

            aggregated_results.append(results)
        except KeyError:
            logger.warning("No results found for %s in %s", parameter, bits['filepath'])

    results = pd.concat(aggregated_results)

    write_results(results, output_file, True)


strategy = Strategy()
config = strategy.results_config.copy()
config.update(strategy.input_config)
logger.debug(config)
config = add_dtypes(config)

input_files : List = snakemake.input
output_file = snakemake.output[0]
parameter = tuple(snakemake.params['parameter'].split(","))
main(input_files, output_file, parameter, config)