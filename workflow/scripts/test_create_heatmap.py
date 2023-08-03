from create_heatmap import sort_results, add_zeros
import pandas as pd
from pandas.testing import assert_frame_equal
from numpy.testing import assert_allclose
from pytest import fixture

@fixture
def data_full():
    return pd.DataFrame(
        data = {
            "MODELRUN": ["model_0"] * 5,
            "YEAR": [2020, 2021, 2022, 2023, 2024],
            "VALUE": [1, 10, 20, 30, 40]
        }
    ).set_index(['MODELRUN', 'YEAR'])

@fixture
def data_short():
    return pd.DataFrame(
        data = {
            "MODELRUN": ["model_1"] * 3,
            "YEAR": [2022, 2023, 2024],
            "VALUE": [1, 2, 3]
        }
    ).set_index(['MODELRUN', 'YEAR'])
    
@fixture
def data_all_zeros():
    return pd.DataFrame(
        data = {
            "MODELRUN": ["model_2"] * 5,
            "YEAR": [2020, 2021, 2022, 2023, 2024],
            "VALUE": [0, 0, 0, 0, 0]
        }
    ).set_index(['MODELRUN', 'YEAR'])

@fixture
def data_add_zeros():
    return pd.DataFrame(
        data = {
            "MODELRUN": ["model_0", "model_1", "model_3", "model_5"],
            "VALUE": [1, 3, 10, 2]
        }
    )

class TestSortResults:

    def test_sort_results_full(self, data_full):
        
        year = 2020
        actual = sort_results(df=data_full, year=year)

        data = {
            "MODELRUN": ["model_0"],
            "VALUE": [1]
        }
        expected = pd.DataFrame(data).set_index(['MODELRUN']).to_numpy(dtype=float)
        assert_allclose(actual, expected, rtol=1e-5)

    def test_sort_results_short(self, data_short):
        
        year = 2020
        actual = sort_results(df=data_short, year=year)

        data = {
            "MODELRUN": ["model_1"],
            "VALUE": [0]
        }
        expected = pd.DataFrame(data).set_index(['MODELRUN']).to_numpy(dtype=float)
        assert_allclose(actual, expected, rtol=1e-5)
        
    def test_sort_results_zeros(self, data_all_zeros):
        
        year = 2020
        actual = sort_results(df=data_all_zeros, year=year)
        data = {
            "MODELRUN": ["model_2"],
            "VALUE": [0]
        }
        expected = pd.DataFrame(data).set_index(['MODELRUN']).to_numpy(dtype=float)
        assert_allclose(actual, expected, rtol=1e-5)
        
class TestAddZeros:
    
    def test_add_zeros(self, data_add_zeros):
        
        models = ["model_0", "model_1", "model_2", "model_3", "model_4", "model_5"]
        actual = add_zeros(data_add_zeros, models)
        
        expected = pd.DataFrame(
            data = {
                "MODELRUN": ["model_0", "model_1", "model_2", "model_3", "model_4", "model_5"],
                "VALUE": [1, 3, 0, 10, 0, 2]
            }
        )
        print(actual)
        print(expected)
        
        assert_frame_equal(actual, expected)