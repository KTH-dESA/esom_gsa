"""Example data processing script
"""
import sys

from typing import Dict, List, Union

import csv
import pandas as pd
import numpy as np
from otoole.read_strategies import ReadDatapackage
from otoole.write_strategies import WriteDatapackage


def process_data(df: pd.DataFrame, index: List, value: float, 
                 first_year: int, last_year: int) -> pd.DataFrame:
    """Return the interpolated data between min and max years
    """
    # df.index = df.index.sortlevel(level=0)[0]
    df = df.loc[tuple(index + [first_year]):tuple(index + [last_year])]
    df = df.reset_index().set_index('YEAR')
    df.loc[last_year, 'VALUE'] = value
    df.loc[first_year + 1:last_year - 1, 'VALUE'] = np.nan
    result = df.astype({'VALUE':'float'}).interpolate(method='values')
    return result.reset_index().set_index(['REGION', 'TECHNOLOGY', 'YEAR'])


class TestInterpolate:

    def test_interpolate_parameter(self):

        index = "GLOBAL,GCPSOUT0N".split(",")
        value = 6000.0
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
        
        actual = process_data(df, index, value, 2015, 2020)

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


def main(input_filepath, output_filepath, parameters: List[Dict[str, Union[str, int, float]]]):

    reader = ReadDatapackage()
    writer = WriteDatapackage()

    with open(parameters, 'r') as csv_file:
        parameters = list(csv.DictReader(csv_file))

    model_params, default_values = reader.read(input_filepath)

    first_year = model_params['YEAR'].min().values[0]
    end_year = model_params['YEAR'].max().values[0]

    for parameter in parameters:
        name = parameter['name']
        index = parameter['indexes'].split(",")
        value = parameter['value']
        df = model_params[name]
        df.index = df.index.sortlevel(level=0)[0]
        snippet = df.xs(tuple(index), drop_level=False)
        new_values = process_data(snippet, index, value, first_year, end_year)
        model_params[name].update(new_values)

    writer.write(model_params, output_filepath, default_values)


if __name__ == "__main__":

    main(sys.argv[1], sys.argv[2], sys.argv[3])