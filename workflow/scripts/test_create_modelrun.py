from create_modelrun import get_types_from_tuple, modify_parameters, process_data
import pandas as pd
import numpy as np

class TestModifyParameters:

    def test_get_types_from_tuple_1(self):
        index = ["GLOBAL","SOLPV"]
        param = 'DiscountRate'
        var_cost_index = ['REGION','TECHNOLOGY']
        index_dtypes = {'REGION': 'str','TECHNOLOGY': 'str','VALUE': 'float'}
        config = {'DiscountRate': {'indices': var_cost_index,
                                   'index_dtypes': index_dtypes}}
        actual = get_types_from_tuple(index, param, config)
        expected = ["GLOBAL","SOLPV"]
        assert actual == expected

    def test_get_types_from_tuple_2(self):
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
            ["GLOBAL","AEIELEXXX",1,2015,100.0],
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
            ["GLOBAL","AEIELEXXX",1,2015,100.0],
        ], columns=var_cost_cols).astype(
                {'REGION': 'str','TECHNOLOGY': 'str','MODE_OF_OPERATION': 'int','YEAR': 'int','VALUE': 'float'}
            ).set_index(var_cost_index)
        pd.testing.assert_frame_equal(actual[name], expected)

    def test_fixed_none(self):
        name = 'DiscountRate'
        var_cost_cols = ['REGION','VALUE']
        var_cost_index = ['REGION']
        index_dtypes = {'REGION': 'str','VALUE': 'float'}
        data = [
            ["GLOBAL", 0.05],
        ]
        model_params = {
            'DiscountRate': pd.DataFrame(data=data, columns=var_cost_cols).astype(
                index_dtypes
            ).set_index(var_cost_index),
            'YEAR': pd.DataFrame(data=[2015, 2016, 2017, 2018], columns=['VALUE'])}
        parameters = [{'name': 'DiscountRate',
                       'indexes': "GLOBAL",
                       'value_base_year':  0.1,
                       'value_end_year': 0.1,
                       'dist': 'unif',
                       'interpolation_index': 'None',
                       'action': 'fixed'}]

        config = {'DiscountRate': {'indices': var_cost_index,
                                   'index_dtypes': index_dtypes}}
        actual = modify_parameters(model_params, parameters, config)
        expected = pd.DataFrame(data=[
            ["GLOBAL", 0.1],
        ], columns=var_cost_cols).astype(
                {'REGION': 'str', 'VALUE': 'float'}
            ).set_index(var_cost_index)
        pd.testing.assert_frame_equal(actual[name], expected)

    def test_fixed_2_none(self):
        name = 'DiscountRate'
        var_cost_cols = ['REGION', 'TECHNOLOGY', 'VALUE']
        var_cost_index = ['REGION', 'TECHNOLOGY']
        index_dtypes = {'REGION': 'str', 'TECHNOLOGY': 'str', 'VALUE': 'float'}
        data = [
            ["GLOBAL", "SOLPV", 0.05],
        ]
        model_params = {
            'DiscountRate': pd.DataFrame(data=data, columns=var_cost_cols).astype(
                index_dtypes
            ).set_index(var_cost_index),
            'YEAR': pd.DataFrame(data=[2015, 2016, 2017, 2018], columns=['VALUE'])}
        parameters = [{'name': 'DiscountRate',
                       'indexes': "GLOBAL,SOLPV",
                       'value_base_year':  0.1,
                       'value_end_year': 0.1,
                       'dist': 'unif',
                       'interpolation_index': 'None',
                       'action': 'fixed'}]

        config = {'DiscountRate': {'indices': var_cost_index,
                                   'index_dtypes': index_dtypes}}
        actual = modify_parameters(model_params, parameters, config)
        expected = pd.DataFrame(data=[
            ["GLOBAL", "SOLPV", 0.1],
        ], columns=var_cost_cols).astype(index_dtypes
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


    def test_interpolate_negative(self):
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
                       'value_base_year':  -100.0,
                       'value_end_year': -1.0,
                       'dist': 'unif',
                       'interpolation_index': 'YEAR',
                       'action': 'interpolate'}]

        config = {'VariableCost': {'indices': var_cost_index,
                                   'index_dtypes': index_dtypes}}
        actual = modify_parameters(model_params, parameters, config)
        expected = pd.DataFrame(data=[
            ["GLOBAL","AEIELEI0H",1,2015,-100.0],
            ["GLOBAL","AEIELEI0H",1,2016,-67.0],
            ["GLOBAL","AEIELEI0H",1,2017,-34.0],
            ["GLOBAL","AEIELEI0H",1,2018,-1.0],
        ], columns=var_cost_cols).astype(
                index_dtypes
            ).set_index(var_cost_index)
        pd.testing.assert_frame_equal(actual[name], expected)

class TestInterpolate:

    def test_interpolate_parameter(self):

        start_year_value = 1000.0
        end_year_value = 6000.0

        actual = process_data(start_year_value, end_year_value, 2015, 2020)

        expected = np.array([[
            1000.0,
            2000.0,
            3000.0,
            4000.0,
            5000.0,
            6000.0,
        ]]).T

        np.testing.assert_equal(actual, expected)

