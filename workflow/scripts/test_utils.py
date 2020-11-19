from utils import get_model_run_scenario_from_input_filepath

def test_get_model_run_scenario_from_input_filepath():

    data = "results/{{scenario}}/model_{modelrun}/data/{input_file}.csv"
    actual = get_model_run_scenario_from_input_filepath(data)
    expected = {'model_run': '{model_run}', 'scenario': "{{scenario}}",
            'param': '{input_file}', 'filepath': "results/{{scenario}}/model_{modelrun}/data"}
    assert actual == expected