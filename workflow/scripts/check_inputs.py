"""Performs checks on user inputs.

Arguments
---------
<config_file>
    Path to configuration file

To run this script on the command line, use the following::

    python ckeck_inputs.py <config_file>
"""

import sys
from pathlib import Path
from typing import Dict, List, Union, Any
from otoole.read_strategies import ReadCsv
from otoole.utils import _read_file
import csv
import pandas as pd
import yaml
import math
import os

from logging import getLogger

logger = getLogger(__name__)

def check_csv(scenario : int, path: str):
    """Checks for the existance of the scenario data.

    Parameters
    ----------
    scenario : int
        name of scenario
    path: str
        path to csv folder

    Raises
    ------
    FileNotFoundError
        If csv data does not exist
    """
    if not Path(path).is_dir():
        raise FileNotFoundError(
            f"csv directory {path} does not exist for scenario {scenario}."
        )

def check_scenario_file(path : str): 
    """Checks the validity of the scenario file.
    
    Parameters
    ----------
    path: str
        path to scenario.csv file

    Raises
    ------
    FileNotFoundError
        If the scenario file does not exist
    ValueError
        If the headings of the scenario file are incorrect OR
        If the scenario csv file is empty
    TypeError
        If scenario name is not an int
    """

    expected_headings = ['name','description','csv','config']

    if not Path(path).is_file():
        raise FileNotFoundError(
            f"Scenario file {path} does not exist. Create the file "
            f"'resources/scenarios.csv' with the headings {expected_headings}"
        )

    df = pd.read_csv(path)
    if list(df) != expected_headings:
        raise ValueError(
            f"Scenario file {path} not formatted correctly. Expected the "
            f"headings {expected_headings}. Recieved headings {list(df)}"
        )

    if df.empty:
        raise ValueError(f"Scenario file {path} is empty. Must input scenarios.")

    for name in df['name']:
        if not type(name) is int:
            raise TypeError(
                f"All scenario names must be of type int. Scenario {name} is "
                f"not an int"
            )

    for name, csv in zip(df['name'], df['csv']):
        check_csv(name, csv)

def check_model_file(path: str):
    """Checks for existance of model file. 
    
    Parameters
    ----------
    path: str
        path to model file

    Raises
    ------
    FileNotFoundError
        If the `resources/` folder does not exist OR
        If the model file does not exis 
    """
    model = Path(path)
    directory = model.parents[0]
    if not directory.is_dir():
        raise FileNotFoundError(
            f"The directory for the model {model.name} does not exist. "
            f"User provided the directory {directory}. "
            "Ensure the model is nested in the resources/ directory"
        )
    if not model.is_file():
        raise FileNotFoundError(f"The model file {path} does not exist.")

def check_solver(solver: str):
    """Checks for valid solver type. 
    
    Parameters
    ----------
    solver: str
        input solver selection

    Raises
    ------
    ValueError
        If solver selection is not valid
    """
    valid_solvers = ['cbc', 'gurobi', 'cplex']
    if solver not in valid_solvers:
        raise ValueError(
            f"Solver selection of {solver} not valid. Select from {valid_solvers}"
        )

def check_replicates(replicates: int):
    """Checks number of replicates. 
    
    Parameters
    ----------
    replicates: int
        number of replicates

    Raises
    ------
    ValueError
        If number of replicates to large (>50) OR
        If number of replicates to small (<2)
    """

    if replicates < 2:
        raise ValueError(
            f"Must select at least 2 replicates."
            f"User selected {replicates} replicates."
        )
    if replicates < 5: 
        logger.warning(f"User selection of {replicates} replicates is low. "
            "If unsure, 10 replicates is usually a good staring point. ")
    if replicates > 25:
        logger.warning(f"User selection of {replicates} replicates is large "
            "and will result in many model runs. If unsure, 10 replicates is "
            "usually a good staring point.")
    if replicates > 50: 
        raise ValueError(
            f"Must select less then 50 replicates."
            f"User selected {replicates} replicates"
        )

def check_parameter_index_data(
    param : Dict[str, Union[str, int, float]],
    user_config : Dict[str, Union[str, int, float]],
    model_params : Dict[str, pd.DataFrame]
):
    """Checks indices of parameter data for validity. 

    Parameters
    ----------
    param: Dict[str, Union[str, int, float]]
        Parameter data
    user_config : Dict[str, Union[str, int, float]]
        otoole configuration file 
    model_params : Dict[str, pd.DataFrame]
        Input data for model 
    
    Raises
    ------
    ValueError
        If parameter indices do not match config indices
    """

    config_indices = user_config[param['name']]['indices'].copy()
    param_index = param['indexes'].split(',')
    param_interp_index = param['interpolation_index']

    # While this will only search one of the scenarios data, all parameters that 
    # you are performing SA on must be in all the scenarios
    set_data = {}
    for index in config_indices:
        set_data[index] = model_params[index]['VALUE'].astype(str).to_list()

    if param_interp_index: # YEAR
        try:
            config_indices.remove(param_interp_index)
        except ValueError:
            raise ValueError(
                f"Interpolation index {param_interp_index} for {param['name']} "
                f"does not match user config indices for {param['name']}, which "
                f"include {param_index}"
            )

    if len(config_indices) != len(param_index):
        raise ValueError(
            f"The indices provided for {param['name']} of {param_index} do "
            f"not correspond to the following indices specified in the "
            f"datapackge: {config_indices}"
        )

    for num, index in enumerate(config_indices):
        param_data = param_index[num]
        if not param_data in set_data[index]:
            raise ValueError(
                f"Parameter {param['name']} has invalid input index data. "
                f"The value {param_data} does not exist in the models set data "
                f"for {index}, which includes {set_data[index]}"
            )

def check_parameter_interpolation_data(param : Dict[str, Union[str, int, float]]):
    """Checks input parameter interpolation options. 

    Parameters
    ----------
    param: Dict[str, Union[str, int, float]]
        Parameter data
    
    Raises
    ------
    RuntimeError
        If both start year values and end year values are the same OR
        If both start year values are the same and the action is fixed OR
        If action is interpolate and interpolation index is None
    """

    start_values_same = math.isclose(
        float(param['min_value_base_year']), float(param['max_value_base_year'])
    )
    end_values_same = math.isclose(
        float(param['min_value_end_year']), float(param['max_value_end_year'])
    )

    if start_values_same and end_values_same:
        raise RuntimeError(
            f"No Sensitivity analysis will be performed on {param['name']}. "
            f"Both start base year values of {param['min_value_base_year']} and "
            f"{param['max_value_base_year']} are the same, AND both end year "
            f"values of {param['min_value_end_year']} and "
            f"{param['max_value_end_year']} are the same. There must be a range "
            f"applied to at least one of the pair of points."
        )

    if start_values_same and param['action'] == 'fixed':
        raise RuntimeError(
            f"No Sensitivity analysis will be performed on {param['name']}. "
            f"Both start base year values of {param['min_value_base_year']} and "
            f"{param['max_value_base_year']} are the same, AND the interpolation "
            f"action is set to fixed. If you do not include a range in the start "
            f"year values, then the interpolation action must be set to 'interpolate'"
        )

    if param['action'] == 'interpolate' and param['interpolation_index'] != 'YEAR':
        raise RuntimeError(
            f"The interpolation index of {param['name']} is set to {param['action']}, "
            f"while the interpolation_index is set to {param['interpolation_index']}. "
            f"Only indexing over 'YEAR' is currently supported with interpolation. "
            f"Change the interpolation_index for {param['name']} to 'YEAR'"
        )

def is_number(s):
    """Checks if a string is a number.

    Parameters
    ----------
    s : str
        input string to check 

    Returns
    -------
    bool
        True if number, false if not number 
    """
    try:
        float(s)
        return True
    except ValueError:
        return False

def check_parameter_data(
    param : Dict[str, Union[str, int, float]],
    user_config : Dict[str, Union[str, int, float]]
):
    """Checks baisc input parameter data. 

    Parameters
    ----------
    param: Dict[str, Union[str, int, float]]
        Parameter data
    user_config : Dict[str, Union[str, int, float]]
        otoole configuration file 
    
    Raises
    ------
    ValueError
        If parameter name not a valid parameter name OR
        If dist is not 'unif' OR
        If 'interpolation_index' is not 'YEAR' or 'None'
    TypeError
        If group is not a string
        If start/end values are not ints or floats
    """

    try:
        user_config[param['name']]
    except ValueError:
        raise ValueError(
            f"{param['name']} is in parameter file, but not in scenario user_config"
        )

    if param['dist'] != 'unif':
        raise ValueError(
            f"{param['name']} distribution must be 'unif'. "
            f"User inputted {param['dist']}"
        )

    if param['interpolation_index'] not in ['YEAR', 'None', 'none', None]:
        if not param['interpolation_index']:
            pass
        else:
            raise ValueError(
                f"{param['name']} interpolation index must be 'YEAR' or None. "
                f"User inputted {param['interpolation_index']}"
            )
    
    if type(param['group']) != str:
        raise TypeError(
            f"{param['name']} group name must be a string. "
            f"User inputted {param['group']} of type {type(param['group'])}"
        )

    end_points = [
        'min_value_base_year', 
        'max_value_base_year', 
        'min_value_end_year', 
        'max_value_end_year'
    ]
    for end_point in end_points:
        if not ((is_number(param[end_point])) or (is_number(param[end_point]))):
            raise TypeError(
                f"{param['name']} {end_point} must be an number. "
                f"User inputted {param[end_point]} of type {type(param[end_point])}"
            )

def check_parameters(
    csv_dir : str, 
    parameters : List[Dict[str, Union[str, int, float]]], 
    user_config : Dict[str, Union[str, int, float]]
): 
    """Checks parameter file. 

    Parameters
    ----------
    csv_dir : str
        path to master csv data directory 
    parameters: List[Dict[str, Union[str, int, float]]]
        Flattened parameter file 
    user_config : Dict[str, Union[str, int, float]]
        otoole configuration file 
    
    Raises
    ------
    ValueError
        If parameter file not populated OR
        Duplicate input data for parameter name and index
    """
    if not parameters:
        raise ValueError(
            f"Empty parameter file. Create config/parameters.csv formattted as:"
            "\n"
            "name,group,indexes,min_value_base_year,max_value_base_year,min_value_end_year,max_value_end_year,dist,interpolation_index,action"
            "DiscountRate,discountrate,REGION,0.05,0.15,0.05,0.15,unif,None,fixed"
            "CapitalCost,CapitalCost,REGION,HYD1,2100,3100,742,1800,unif,YEAR,interpolate"
        )

    duplicate_check_data = []
    for parameter in parameters:
        duplicate_check_data.append([parameter['name'], parameter['indexes']])
    duplicate_check = pd.DataFrame(duplicate_check_data, columns=['name', 'indexes'])
    duplicates = duplicate_check[duplicate_check.duplicated(keep='last')]
    if not duplicates.empty:
        raise ValueError(
            f"Multiple occurances of the parameter(s) {duplicates['name'].to_list()} "
            f"with the corresponding index values of {duplicates['indexes'].to_list()}"
        )

    # get model parameter definitions 
    model_params, _ = ReadCsv(user_config=user_config).read(csv_dir)

    for parameter in parameters:
        check_parameter_data(parameter, user_config)
        check_parameter_interpolation_data(parameter)
        check_parameter_index_data(parameter, user_config, model_params)

def read_parameters_file(path : str) -> List[Dict[str, Union[str, int, float]]]:
    """Reads in a flattened CSV file

    Parameters
    ----------
    path : str
        file path to parameters.csv

    Returns
    --------
    parameters : List[Dict[str, Union[str, int, float]]]
        Flattened parameters file
    """
    with open(path, 'r') as csv_file:
        parameters = list(csv.DictReader(csv_file))
    return parameters

def check_file_extension(file_name : str, extension : str):
    """Checks the file for the correct extension.

    Parameters
    ----------
    file_name : str
        Name of file
    extension : str
        Expected file extension

    Rasies
    ------
    ValueError
        If the actual file extension is not the expected. 
    """
    _, ending = os.path.splitext(file_name)
    if not extension.startswith('.'):
        extension = f".{extension}"
    if ending != extension:
        raise ValueError(
            f"Input configuration file must be a {extension} file. Recieved the" 
            f"file {file_name} with the extension {ending}"
        )

def check_otoole_inputs(
    actual_data : str, 
    scenario : int,
    user_config: Dict[str, Union[str, int, float]],
):
    """Checks that scenario data matches what snakemake will expect

    Parameters
    ----------
    input_data : str
        path the directory holding otoole csv data
    expected_data : str = 'otoole_files.csv'
        name of the file in /resources holding input CSV data
    scenario : int
        scenario number

    Raises
    ------
    FileNotFoundError
        If the input csvs do not match the csvs in resources/otoole_files.csv
    """
    missing_files = []
    expected_files = [x for x, y in user_config.items() if y["type"] in ["param", "set"]]
    for csv in Path(actual_data).iterdir():
        if csv.stem not in expected_files:
            missing_files.append(f"{csv.stem}.csv")
    if missing_files:
        raise FileNotFoundError(
            f"The following CSV files are missing in the input data for scenario "
            f"{scenario} : {missing_files}"
        )

def main(config : Dict[str, Any]):
    """Runs a series of checks on input data. 

    Parameters
    ----------
    config : Dict[str, Any]
        Parsed config.yaml file

    Rasies
    ------
    ValueError
        If the input config file is not a yaml file
    """

    # check config file
    check_scenario_file(config['scenarios'])
    check_model_file(config['model_file'])
    check_solver(config['solver'])
    check_replicates(config['replicates'])

    # read in csv path and user_config file
    scenario_df = pd.read_csv(config['scenarios'])
    scenarios = scenario_df['name'].to_list()
    csvs = scenario_df['csv'].to_list()
    configs = scenario_df['config'].to_list()
    parameters = read_parameters_file(config['parameters'])
    with open(configs[0], "r") as f:
        user_config = _read_file(f, ".yaml")
    
    # check parameter file 
    check_parameters(csvs[0], parameters, user_config)

    # check otoole inputs 
    for scenario, csv in zip(scenarios, csvs):
        data_dir = Path(csv)
        check_otoole_inputs(str(data_dir), scenario, user_config)

if __name__ == "__main__":

    if len(sys.argv) != 2:
        raise ValueError(
            "Usage: python create_modelrun.py <config_file>"
        )
    else:
        with open(sys.argv[1], 'r') as stream:
            try:
                parsed_yaml = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
                raise exc
        main(parsed_yaml)