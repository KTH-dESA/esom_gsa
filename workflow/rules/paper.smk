"""Rules to retrieve and run results for the submitted paper"""

MODELS = ["1a", "1b", "2a", "2b", "3a", "3b"]
MODEL_NUMS = [1,2,3]
GSA_FILES = ["parameters", "results", "scenarios"]

wildcard_constraints:
    model = "[a-z0-9]+"

rule setup_paper_data:
    input:
        csvs = expand("config/model_{model}/data/{input_file}.csv", model = MODELS, input_file=INPUT_FILES),
        model = expand("resources/model_{model}/osemosys_fast.txt", model = MODELS),
        otoole_config = expand("config/model_{model}/otoole.yaml", model = MODELS),
        parameters = expand("config/model_{model}/parameters.csv", model = MODELS),
        results = expand("config/model_{model}/results.csv", model = MODELS),
        scenarios = expand("config/model_{model}/scenarios.csv", model = MODELS),
        config = expand("config/model_{model}/config.yaml", model = MODELS),
        notebooks = expand("workflow/notebooks/model_{model_num}.ipynb", model_num = MODEL_NUMS),

rule retrieve_zenodo_data:
    params: 
        # doi = "https://zenodo.org/records/10107849",
        doi = "10107849",
        save = "./" # root folder
    output:
        csvs = expand("paper/model_{model}/model/data/{input_file}.csv", model = MODELS, input_file=INPUT_FILES),
        model = expand("paper/model_{model}/model/osemosys_fast.txt", model = MODELS),
        otoole_config = expand("paper/model_{model}/model/config.yaml", model = MODELS),
        gsa_files = expand("paper/model_{model}/config/{gsa_file}.csv", model = MODELS, gsa_file=GSA_FILES),
        gsa_config = expand("paper/model_{model}/config/config.yaml", model = MODELS),
        notebooks = expand("paper/notebooks/model_{model_num}_gsa.ipynb", model_num = MODEL_NUMS)
    log:
        "logs/retrieve_zenodo_data_model.log",
    script:
        "../scripts/retrieve_zenodo_data.py"

rule create_paper_model_directory:
    output:
        directory("config/model_{model}")
    run:
        from pathlib import Path

        dir_path = Path(output)
        if not dir_path.exists():
            dir_path.mkdir()

rule copy_paper_model_data:
    input:
        csvs = expand("paper/model_{{model}}/model/data/{input_file}.csv", input_file=INPUT_FILES),
        otoole_config = "paper/model_{model}/model/config.yaml"
    output:
        csvs = expand("config/model_{{model}}/data/{input_file}.csv", input_file=INPUT_FILES),
        otoole_config = "config/model_{model}/otoole.yaml"
    run:
        import shutil 
        from pathlib import Path

        dir_path = Path("config", f"model_{wildcards.model}", "data")
        if not dir_path.exists():
            dir_path.mkdir(parents=True)

        for csv in input.csvs:
            csv_name = Path(csv).name
            shutil.copy(csv, str(Path(dir_path, csv_name)))

        shutil.copy(input.otoole_config, output.otoole_config)

rule copy_paper_model:
    input:
        model = "paper/model_{model}/model/osemosys_fast.txt",
    output:
        model = "resources/model_{model}/osemosys_fast.txt",
    run: 
        import shutil
        dir_path = Path("resources", f"model_{wildcards.model}")
        if not dir_path.exists():
            dir_path.mkdir()
        shutil.copy(input.model, output.model)

rule copy_paper_gsa_data:
    input:
        parameters = "paper/model_{model}/config/parameters.csv",
        results = "paper/model_{model}/config/results.csv",
        scenarios = "paper/model_{model}/config/scenarios.csv",
        config = "paper/model_{model}/config/config.yaml",
    output:
        parameters = "config/model_{model}/parameters.csv",
        results = "config/model_{model}/results.csv",
        scenarios = "config/model_{model}/scenarios.csv",
        config = "config/model_{model}/config.yaml",
    run:
        import shutil 
        from pathlib import Path

        dir_path = Path("config", f"model_{wildcards.model}")
        if not dir_path.exists():
            dir_path.mkdir()

        shutil.copy(input.parameters, output.parameters)
        shutil.copy(input.results, output.results)
        shutil.copy(input.scenarios, output.scenarios)
        shutil.copy(input.config, output.config)

rule copy_paper_notebooks:
    input:
        notebook = "paper/notebooks/model_{model_num}_gsa.ipynb"
    output:
        notebook = "workflow/notebooks/model_{model_num}.ipynb"
    wildcard_constraints:
        model_num = "[0-9]"
    run:
        import shutil 
        shutil.copy(input.notebook, output.notebook)