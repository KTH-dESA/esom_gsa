from extract_results import default_value_xs
import pandas as pd
from pandas.testing import assert_frame_equal
from pytest import fixture

@fixture
def data_in():
    return pd.DataFrame(
        data = {
            "REGION": ["R1"] * 5,
            "TECHNOLOGY": ["FOSSIL"] * 5,
            "YEAR": [2020, 2021, 2022, 2023, 2024],
            "VALUE": [5, 10, 15, 20, 25]
        }
        ).set_index(['REGION', 'TECHNOLOGY', 'YEAR'])

class TestDefaultValueXs:

    def test_default_value_xs(self, data_in):
        
        index = ('REGION', 'TECHNOLOGY')
        parameters = ('R1', 'HYDRO')
        
        actual = default_value_xs(df=data_in, index=index, parameters=parameters)

        expected_data = {
            "REGION": ["R1"] * 5,
            "TECHNOLOGY": ["HYDRO"] * 5,
            "YEAR": [2020, 2021, 2022, 2023, 2024],
            "VALUE": [0, 0, 0, 0, 0]
        }
        expected = pd.DataFrame(expected_data).set_index(['REGION', 'TECHNOLOGY', 'YEAR'])
        
        assert_frame_equal(actual, expected)
