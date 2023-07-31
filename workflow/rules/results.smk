def get_input(wildcards):
    input_file = RESULTS.set_index('filename').loc[wildcards.result_file]['resultfile']
    return ["results/{wildcards.scenario}/model_{modelrun}/results/{input_file}.csv".format(
        modelrun=x, input_file=input_file, wildcards=wildcards) for x in MODELRUNS]

def get_sample_name(wildcards):
    if config['scale']:
        return f"modelruns/{wildcards.scenario}/morris_sample_scaled.txt"
    else: 
        return f"modelruns/{wildcards.scenario}/morris_sample.txt"
        
def get_indices(wildcards):
    indices = RESULTS.set_index('filename').loc[wildcards.result_file].dropna().drop('resultfile').to_dict()
    return {x:str(indices[x]).split(',') for x in indices}

rule extract_results:
    input: 
        csvs=get_input,
        config=config_from_scenario,
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
        result_type='objective',
        scaled = config['scale']
    input: 
        sample = get_sample_name,
        results = "results/{scenario}/objective_{scenario}.csv"
    output: 
        expand("results/{{scenario}}_summary/SA_objective.{ext}",ext=['csv','png'])
    conda: "../envs/sample.yaml"
    shell: "python workflow/scripts/calculate_SA_results.py {params.parameters} {input.sample} {input.results} {output[0]} {params.result_type} {params.scaled}"

rule calculate_SA_user_defined:
    message:
        "Calcualting user defined sensitivity measures"
    params: 
        parameters=config['parameters'],
        result_type='variable',
        scaled = config['scale']
    input: 
        sample = get_sample_name,
        results=expand("results/{{scenario}}/{{result_file}}.{ext}", ext=config['filetype'])
    output: 
        expand("results/{{scenario}}_summary/SA_{{result_file}}.{ext}",ext=['csv','png'])
    conda: "../envs/sample.yaml"
    shell: "python workflow/scripts/calculate_SA_results.py {params.parameters} {input.sample} {input.results} {output[0]} {params.result_type} {params.scaled}"

rule create_heatmap:
    message: 
        "Calculating user defined sensitivity measures"
    params: 
        parameters=config['parameters'],
        scaled = config['scale']
    input:
        sample="modelruns/{scenario}/morris_sample.txt",
        results=expand("results/{{scenario}}/{{result_file}}.{ext}", ext=config['filetype'])
    output:
        "results/{scenario}_summary/{result_file}_heatmap.png"
    shell: 
        "python workflow/scripts/create_heatmap.py {params.parameters} {input.sample} {input.results} {output} {params.scaled}"

rule plot_interactions:
    message:
        "Creating interaction plots"
    params:
        parameters=config['parameters']
    input:
        sample = "modelruns/{scenario}/morris_sample.txt",
        results = "results/{scenario}/objective_{scenario}.csv"
    output:
        "results/{scenario}_summary/SA_interactions.png"
    shell:
        "python workflow/scripts/plot_interactions.py {params.parameters} {input.sample} {input.results} {output}"

rule get_status:
    params:
        modelruns=len(MODELRUNS),
        scenarios=len(SCENARIOS)
    input: expand("temp/{scenario}/model_{model_run}.sol", model_run=MODELRUNS, scenario=SCENARIOS.index)
    output: "results/status.csv"
    log: "results/log/status.log"
    shell:
        """
        if [ {config[solver]} = gurobi ]
        then
            python workflow/scripts/run_status.py {output} {params.modelruns} {params.scenarios}
        elif [ {config[solver]} = cplex ]
            python workflow/scripts/run_status.py {output} {params.modelruns} {params.scenarios}
        else
            #! /bin/bash -x
            echo "SCENARIO,FILE,OBJECTIVE,STATUS" > {output}
            for FILE in ({input})
            do
            OBJ=$(head $FILE | grep -e 'objectiveValue' | cut -f 2 -d '=')
            STATUS=$(head $FILE | grep -e 'solutionStatusString' | cut -f 2 -d '=')
            JOB=$(echo $FILE | cut -f 3 -d '/' | cut -f 1 -d '.')
            echo "1,$JOB,$OBJ,$STATUS" >> {output}
            done
        fi
        """