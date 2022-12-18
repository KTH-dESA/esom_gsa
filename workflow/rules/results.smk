def get_input(wildcards):
    # input_file = RESULTS.set_index('filename').loc[wildcards.result_file]['resultfile']
    # return ["results/{scenario}/model_{modelrun}/results/{input_file}.csv".format(
    #     modelrun=x, input_file=input_file, scenario=y) for x in MODELRUNS for y in SCENARIOS.index]
    input_file = RESULTS.set_index('filename').loc[wildcards.result_file]['resultfile']
    return ["results/{wildcards.scenario}/model_{modelrun}/results/{input_file}.csv".format(
        modelrun=x, input_file=input_file, wildcards=wildcards) for x in MODELRUNS]

def get_indices(wildcards):
    indices = RESULTS.set_index('filename').loc[wildcards.result_file].dropna().drop('resultfile').to_dict()
    return {x:str(indices[x]).split(',') for x in indices}

rule extract_results:
    input: 
        csvs=get_input,
        config=config_from_scenario,
        # csvs=expand("results/{{scenario}}/model_{{model_run}}/results/{csv}.csv", csv=OUTPUT_FILES)
    params:
        parameter = get_indices,
        folder=directory("results/{{scenario}}_summary/")
    log: "results/log/extract_scenarion_{scenario}_{result_file}.log"
    output: expand("results/{{scenario}}/{{result_file}}.{ext}", ext=config['filetype'])
    conda: "../envs/otoole.yaml"
    script: "../scripts/extract_results.py"

rule calculate_hourly_demand:
    input: expand("results/annual_demand.{ext}", ext=config['filetype'])
    output: expand("results/hourly_demand.{ext}", ext=config['filetype'])
    conda: "envs/pandas.yaml"
    script: "scripts/calculate_hourly_demand.py"

rule calculate_hourly_generation:
    input: expand("results/annual_generation.{ext}", ext=config['filetype'])
    output: expand("results/hourly_generation.{ext}", ext=config['filetype'])
    conda: "envs/pandas.yaml"
    script: "scripts/calculate_hourly_generation.py"

rule calculate_SA_objective:
    message:
        "Calcualting objective cost sensitivity measures"
    params: 
        parameters=config['parameters'],
        result_type='objective'
    input: 
        sample = "modelruns/{scenario}/morris_sample.txt",
        results = "results/{scenario}/objective_{scenario}.csv"
    output: 
        expand("results/{{scenario}}_summary/SA_objective.{ext}",ext=['csv','png'])
    conda: "../envs/sample.yaml"
    shell: "python workflow/scripts/calculate_SA_results.py {params.parameters} {input.sample} {input.results} {output[0]} {params.result_type}"

rule calculate_SA_user_defined:
    message:
        "Calcualting user defined sensitivity measures"
    params: 
        parameters=config['parameters'],
        result_type='variable'
    input: 
        sample = "modelruns/{scenario}/morris_sample.txt",
        results=expand("results/{{scenario}}/{{result_file}}.{ext}", ext=config['filetype'])
    output: 
        expand("results/{{scenario}}_summary/SA_{{result_file}}.{ext}",ext=['csv','png'])
    conda: "../envs/sample.yaml"
    shell: "python workflow/scripts/calculate_SA_results.py {params.parameters} {input.sample} {input.results} {output[0]} {params.result_type}"

rule create_heatmap:
    message: 
        "Calculating user defined sensitivity measures"
    params: 
        parameters=config['parameters']
    input:
        sample="modelruns/{scenario}/morris_sample.txt",
        results=expand("results/{{scenario}}/{{result_file}}.{ext}", ext=config['filetype'])
    output:
        "results/{scenario}_summary/{result_file}_heatmap.png"
    shell: "python workflow/scripts/create_heatmap.py {params.parameters} {input.sample} {input.results} {output}"

    
