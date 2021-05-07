"""Creates model runs from a sample file and master model datapackage

Arguments
---------
<input_filepath>
    Path to the master model datapackage
<output_filepath>
    Path to the new model file
<sample_filepath>
    Path the sample file

The expected format of the input sample file is a CSV file with the following structure::

    name,indexes,value_base_year,value_end_year,action,interpolation_index
    CapitalCost,"GLOBAL,GCPSOUT0N",1024.3561663863075,2949.23939,interpolate,YEAR

It is very similar to the overall ``parameter.csv`` configuration file, except holds a sample value
rather than a range

To run this script on the command line, use the following::

    python create_modelrun.py <input_filepath> <output_filepath> <sample_filepath>

"""
import sys

from typing import Dict, List, Union, Tuple

import csv
import pandas as pd
import numpy as np
from otoole.read_strategies import ReadDatapackage
from otoole.write_strategies import WriteDatapackage

from logging import getLogger

logger = getLogger(__name__)

def process_data(df: pd.DataFrame, index: List,
                 start_year_value: float,end_year_value: float,
                 first_year: int, last_year: int
                 ) -> pd.DataFrame:
    """Interpolate data between min and max years

    Arguments
    ---------
    df: pd.DataFrame
    index: List[str]
        List of index names e.g. ``['GLOBAL', 'GCPSOUT0N']``
    start_year_value: float
        Value of the parameter in the start year
    end_year_value: float
        Value of the parameter in the end year
    first_year: int
        First year of the range to interpolate
    last_year: int
        Last year of the range to interpolate
    """
    # df.index = df.index.sortlevel(level=0)[0]

    values = np.interp([range(int(first_year), int(last_year) + 1)],
                       np.array([int(first_year), int(last_year)]),
                       np.array([float(start_year_value), float(end_year_value)])).T

    df.loc[tuple(index + [first_year]):tuple(index + [last_year])] = values
    return df.loc[tuple(index + [first_year]):tuple(index + [last_year])]


class TestInterpolate:

    def test_interpolate_parameter(self):

        index = "GLOBAL,GCPSOUT0N".split(",")
        start_year_value = 1000.0
        end_year_value = 6000.0
        data = [
            ['GLOBAL', 'GCPSOUT0N', 2015, 1000.0],
            ['GLOBAL', 'GCPSOUT0N', 2016, 1000.0],
            ['GLOBAL', 'GCPSOUT0N', 2017, 1000.0],
            ['GLOBAL', 'GCPSOUT0N', 2018, 1000.0],
            ['GLOBAL', 'GCPSOUT0N', 2019, 1000.0],
            ['GLOBAL', 'GCPSOUT0N', 2020, 1000.0],
            ['GLOBAL', 'GCPXXXT0N', 2020, 1000.0],
        ]

        df = pd.DataFrame(data, columns=['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']
                          ).set_index(['REGION', 'TECHNOLOGY', 'YEAR']).astype({'VALUE': float})

        actual = process_data(df, index, start_year_value, end_year_value, 2015, 2020)

        data = [
            ['GLOBAL', 'GCPSOUT0N', 2015, 1000.0],
            ['GLOBAL', 'GCPSOUT0N', 2016, 2000.0],
            ['GLOBAL', 'GCPSOUT0N', 2017, 3000.0],
            ['GLOBAL', 'GCPSOUT0N', 2018, 4000.0],
            ['GLOBAL', 'GCPSOUT0N', 2019, 5000.0],
            ['GLOBAL', 'GCPSOUT0N', 2020, 6000.0],
        ]

        expected = pd.DataFrame(data, columns=['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']
                          ).set_index(['REGION', 'TECHNOLOGY', 'YEAR'])

        pd.testing.assert_frame_equal(actual, expected)


def get_types_from_tuple(index: list, param: str, config: Dict) -> Tuple:
    depth = len(index)
    names = config[param]['indices'][0:depth]
    typed_index = []
    dtypes = config[param]['index_dtypes']
    for name, element in zip(names, index):
        this_type = dtypes[name]
        if this_type == 'str':
            typed_index.append(str(element))
        elif this_type == 'float':
            typed_index.append(float(element))
        elif this_type == 'int':
            typed_index.append(int(element))

    return typed_index


def modify_parameters(
        model_params: Dict[str, pd.DataFrame],
        parameters: List[Dict[str, Union[str, int, float]]],
        config: Dict):
    """
    """
    first_year = model_params['YEAR'].min().values[0]
    end_year = model_params['YEAR'].max().values[0]

    for parameter in parameters:

        name = parameter['name']
        df = model_params[name]
        untyped_index = parameter['indexes'].split(",")
        index = get_types_from_tuple(untyped_index, name, config)
        start_year_value = parameter['value_base_year']
        end_year_value = parameter['value_end_year']
        action = parameter['action']
        inter_index = parameter['interpolation_index']
        if action == 'interpolate':
            df2 = df.sort_index()
            snippet = df2.xs(tuple(index), drop_level=False)
            new_values = process_data(snippet, index, start_year_value, end_year_value, first_year, end_year)
        elif action == 'fixed':
            if inter_index == 'None':
                # Create new object and update inplace
                data = [index + [start_year_value]]
            elif inter_index == 'YEAR':
                # Create new object and update inplace
                data = [index + [x] + [start_year_value] for x in range(first_year, end_year + 1, 1)]
            columns = config[name]['indices']
            new_values = pd.DataFrame(data, columns=columns + ['VALUE']).astype(config[name]['index_dtypes']).set_index(columns)
        if all(new_values.index.isin(model_params[name].index)):
            logger.info("Updating values for {} in {}".format(index, name))
            model_params[name].update(new_values)
        else:
            logger.info("Appending values for {} in {}".format(index, name))
            model_params[name] = model_params[name].append(new_values)

    return model_params


class TestModifyParameters:

    def test_get_types_from_tuple(self):
        index = ["GLOBAL","AEIELEI0H","1"]
        param = 'VariableCost'
        var_cost_index = ['REGION','TECHNOLOGY','MODE_OF_OPERATION','YEAR']
        index_dtypes = {'REGION': 'str','TECHNOLOGY': 'str','MODE_OF_OPERATION': 'int','YEAR': 'int','VALUE': 'float'}
        config = {'VariableCost': {'indices': var_cost_index,
                                   'index_dtypes': index_dtypes}}
        actual = get_types_from_tuple(index, param, config)
        expected = ["GLOBAL","AEIELEI0H",1]
        assert actual == expected

    def test_fixed_year(self):
        name = 'VariableCost'
        var_cost_cols = ['REGION','TECHNOLOGY','MODE_OF_OPERATION','YEAR','VALUE']
        var_cost_index = ['REGION','TECHNOLOGY','MODE_OF_OPERATION','YEAR']
        index_dtypes = {'REGION': 'str','TECHNOLOGY': 'str','MODE_OF_OPERATION': 'int','YEAR': 'int','VALUE': 'float'}
        data = [
            ["GLOBAL","AEIELEI0H",1,2015,100.0],
            ["GLOBAL","AEIELEI0H",1,2016,100.0],
            ["GLOBAL","AEIELEI0H",1,2017,100.0],
            ["GLOBAL","AEIELEI0H",1,2018,100.0],
        ]

        model_params = {
            'VariableCost': pd.DataFrame(data=data, columns=var_cost_cols).astype(
                index_dtypes
            ).set_index(var_cost_index),
            'YEAR': pd.DataFrame(data=[2015, 2016, 2017, 2018], columns=['VALUE'])}
        parameters = [{'name': 'VariableCost',
                       'indexes': "GLOBAL,AEIELEI0H,1",
                       'value_base_year':  1,
                       'value_end_year': 1,
                       'dist': 'unif',
                       'interpolation_index': 'YEAR',
                       'action': 'fixed'}]

        config = {'VariableCost': {'indices': var_cost_index,
                                   'index_dtypes': index_dtypes}}
        actual = modify_parameters(model_params, parameters, config)
        expected = pd.DataFrame(data=[
            ["GLOBAL","AEIELEI0H",1,2015,1.0],
            ["GLOBAL","AEIELEI0H",1,2016,1.0],
            ["GLOBAL","AEIELEI0H",1,2017,1.0],
            ["GLOBAL","AEIELEI0H",1,2018,1.0],
        ], columns=var_cost_cols).astype(
                {'REGION': 'str','TECHNOLOGY': 'str','MODE_OF_OPERATION': 'int','YEAR': 'int','VALUE': 'float'}
            ).set_index(var_cost_index)
        pd.testing.assert_frame_equal(actual[name], expected)

    def test_interpolate(self):
        name = 'VariableCost'
        var_cost_cols = ['REGION','TECHNOLOGY','MODE_OF_OPERATION','YEAR','VALUE']
        var_cost_index = ['REGION','TECHNOLOGY','MODE_OF_OPERATION','YEAR']
        index_dtypes = {'REGION': 'str','TECHNOLOGY': 'str','MODE_OF_OPERATION': 'int','YEAR': 'int','VALUE': 'float'}
        data = [
            ["GLOBAL","AEIELEI0H",1,2015,100.0],
            ["GLOBAL","AEIELEI0H",1,2016,100.0],
            ["GLOBAL","AEIELEI0H",1,2017,100.0],
            ["GLOBAL","AEIELEI0H",1,2018,100.0],
        ]

        model_params = {
            'VariableCost': pd.DataFrame(data=data,
                                         columns=var_cost_cols
                ).astype(index_dtypes
                ).set_index(var_cost_index),
            'YEAR': pd.DataFrame(data=[2015, 2016, 2017, 2018], columns=['VALUE'])}
        parameters = [{'name': 'VariableCost',
                       'indexes': "GLOBAL,AEIELEI0H,1",
                       'value_base_year':  100.0,
                       'value_end_year': 1.0,
                       'dist': 'unif',
                       'interpolation_index': 'YEAR',
                       'action': 'interpolate'}]

        config = {'VariableCost': {'indices': var_cost_index,
                                   'index_dtypes': index_dtypes}}
        actual = modify_parameters(model_params, parameters, config)
        expected = pd.DataFrame(data=[
            ["GLOBAL","AEIELEI0H",1,2015,100.0],
            ["GLOBAL","AEIELEI0H",1,2016,67.0],
            ["GLOBAL","AEIELEI0H",1,2017,34.0],
            ["GLOBAL","AEIELEI0H",1,2018,1.0],
        ], columns=var_cost_cols).astype(
                index_dtypes
            ).set_index(var_cost_index)
        pd.testing.assert_frame_equal(actual[name], expected)


def main(input_filepath, output_filepath, parameters: List[Dict[str, Union[str, int, float]]]):

    reader = ReadDatapackage()

    writer = WriteDatapackage()

    logger.info("Reading datapackage {}".format(input_filepath))
    model_params, default_values = reader.read(input_filepath)
    config = reader.input_config
    model_params = modify_parameters(model_params, parameters, config)
    writer.write(model_params, output_filepath, default_values)


if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage: python create_modelrun.py <input_filepath> <output_filepath> <sample_filepath>")
    else:
        with open(sys.argv[3], 'r') as csv_file:
            sample = list(csv.DictReader(csv_file))
        main(sys.argv[1], sys.argv[2], sample)
