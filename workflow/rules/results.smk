def get_input(wildcards):
    input_file = AGG_RESULTS.set_index('filename').loc[wildcards.agg_result_file]['resultfile']
    return ["results/{scenario}/{modelrun}/{input_file}.csv".format(
        modelrun=x, input_file=input_file, scenario=y) for x in MODELRUNS for y in SCENARIOS.index]

def get_indices(wildcards):
    return AGG_RESULTS.set_index('filename').loc[wildcards.agg_result_file]['indices']

rule extract_results:
    input: get_input
    params:
        parameter = get_indices
    log: "results/log/extract_{agg_result_file}.log"
    output: expand("results/{{agg_result_file}}.{ext}", ext=config['filetype'])
    conda: "envs/otoole.yaml"
    script: "scripts/extract_results.py"

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

rule calculate_SA_results:
    params: 
        parameters = config['parameters']
    input: 
        sample = "modelruns/{scenario}/morris_sample.txt",
        results = "results/{scenario}/objective_{scenario}.csv"
    output: 
        SA_csv = "results/SA_{scenario}.csv",
        SA_png = "results/SA_{scenario}.png"
    conda: "../envs/sample.yaml"
    shell: "python workflow/scripts/analyze_results.py {params.parameters} {input.sample} {input.results} {output.SA_csv}"