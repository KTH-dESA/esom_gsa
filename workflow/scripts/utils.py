import os
import pandas as pd
from typing import List

from logging import getLogger

logger = getLogger(__name__)

def get_model_run_scenario_from_input_filepath(filename: str):
    """Parses filepath to extract useful bits

    "results/{{scenario}}/model_{modelrun}/data/{input_file}.csv"
    """
    filepath, name = os.path.split(filename)
    param = os.path.splitext(name)[0]
    scenario_path, model_run = os.path.split(filepath)
    resultsscenario, _ = os.path.split(scenario_path)
    scenario = os.path.split(resultsscenario)[1]
    return {'model_run': model_run, 'scenario': scenario,
            'param': param, 'filepath': filepath}

def get_model_run_scenario_from_filepath(filename: str):
    """Parses filepath to extract useful bits

    "results/{{scenario}}/{modelrun}/{input_file}.csv"
    """
    filepath, name = os.path.split(filename)
    param = os.path.splitext(name)[0]
    scenario_path, model_run = os.path.split(filepath)
    scenario = os.path.split(scenario_path)[1]
    return {'model_run': model_run, 'scenario': scenario,
            'param': param, 'filepath': filepath}

def read_results(input_filepath: str) -> pd.DataFrame:
    extension = os.path.splitext(input_filepath)[1]
    if extension == '.parquet':
        df = pd.read_parquet(input_filepath)
    elif extension == '.csv':
        df = pd.read_csv(input_filepath)
    elif extension == '.feather':
        df = pd.read_feather(input_filepath)
    return df


def write_results(df: pd.DataFrame, output_filepath: str, index=None) -> None:
    """Write out aggregated results to disk

    Arguments
    ---------
    df: pd.DataFrame
        Dataframe to write out
    output_filepath: str
        Path to the output file
    index=None
        Whether to write out the index or not
    """
    extension = os.path.splitext(output_filepath)[1]
    if extension == '.parquet':
        df.to_parquet(output_filepath, index=index)
    elif extension == '.csv':
        df.to_csv(output_filepath, index=index)
    elif extension == '.feather':
        if index:
            df = df.reset_index()
        df.to_feather(output_filepath)

def create_salib_problem(parameters: List) -> dict:
    """Creates SALib problem from scenario configuration file.
    
    Arguments
    ---------
    parameters: List
        List of dictionaries describing problem. Each dictionary must have
        'name', 'indexes', 'group' keys

    Returns
    -------
    problem: dict
        SALib formatted problem dictionary

    Raises
    ------
    ValueError
        If only one variable is givin, OR 
        If only one group is given
    """

    problem = {}
    problem['num_vars'] = len(parameters)
    if problem['num_vars'] <= 1:
        logger.error(f"Must define at least two variables in problem. User defined {problem['num_vars']} variable(s).")
        raise ValueError

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
    num_groups = len(set(groups))
    if num_groups <= 1:
        logger.error(f"Must define at least two groups in problem. User defined {num_groups} group(s).")
        raise ValueError

    return problem