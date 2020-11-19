import os

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