import json

wildcard_constraints:
    modelrun="\d+",
    scenarios="\d+"

def datapackage_from_scenario(wildcards):
    return SCENARIOS.loc[int(wildcards.scenario), 'path']

rule copy_datapackage:
    message: "Copying and modifying datapackage for '{output.folder}'"
    input:
        datapackage=datapackage_from_scenario,
        sample="modelruns/{scenario}/{model_run}_sample.txt"
    log: "results/copy_datapackage_{scenario}_{model_run}.log"
    conda: "../envs/otoole.yaml"
    output:
        folder=directory("results/{scenario}/model_{model_run, \d+}"),
        dummy="results/{scenario}/model_{model_run, \d+}/datapackage.json",
    shell:
        "python workflow/scripts/create_modelrun.py {input.datapackage} {output.folder} {input.sample}"

rule generate_datafile:
    message: "Generating datafile for '{output}'"
    input:
        datapackage="results/{scenario}/model_{model_run}/datapackage.json"
    output:
        "modelruns/{scenario}/model_{model_run}.txt"
    conda: "../envs/otoole.yaml"
    log:
        "results/otoole_{scenario}_{model_run}.log"
    shell:
        "otoole -v convert datapackage datafile {input} {output} 2> {log}"

rule generate_lp_file:
    message: "Generating the LP file for '{output}'"
    input:
        data=expand("modelruns/{{scenario}}/model_{{model_run}}.txt"),
        model=config['model_file']
    resources:
        mem_mb=7000,
        disk_mb=1300
    output:
        temporary("results/{scenario}/{model_run}.lp")
    benchmark:
        "benchmarks/gen_lp/{scenario}_{model_run}.tsv"
    log:
        "results/glpsol_{scenario}_{model_run}.log"
    conda: "../envs/osemosys.yaml"
    threads:
        1
    shell:
        "glpsol -m {input.model} -d {input.data} --wlp {output} --check > {log}"

rule solve_lp:
    message: "Solving the LP for '{output}'"
    input:
        "results/{scenario}/{model_run}.lp"
    output:
        protected("results/{scenario}/{model_run}.sol")
    log:
        "results/cbc_{scenario}_{model_run}.log"
    benchmark:
        "benchmarks/cbc/{scenario}_{model_run}.tsv"
    resources:
        mem_mb=3000,
        disk_mb=33
    threads:
        1
    shell:
        "cbc {input} solve -sec 1500 -solu {output} > {log}"

rule process_solution:
    message: "Processing CBC solution for '{output}'"
    input:
        solution="results/{scenario}/{model_run}.sol",
        data="results/{scenario}/model_{model_run}/datapackage.json"
    output: ["results/{{scenario}}/{{model_run, \d+}}/{}.csv".format(x) for x in RESULTS.index]
    log: "results/process_solution_{scenario}_{model_run}.log"
    params:
        folder="results/{scenario}/{model_run}"
    shell:
        "mkdir -p {params.folder} && otoole -v results cbc csv {input.solution} {params.folder} --input_datapackage {input.data} 2> {log}"

# rule solve_gurobi:
#     message: "Solving the LP for '{output}' using Gurobi"
#     input:
#         "results/{scenario}/{model_run}.lp"
#     output:
#         json="results/{scenario}/{model_run}.json",
#         solution="results/{scenario}/{model_run}.sol",
#     log:
#         "results/solver_{scenario}_{model_run}.log"
#     params:
#         ilp="results/{scenario}/{model_run}.ilp"
#     benchmark:
#         "benchmarks/gurobi/{scenario}_{model_run}.tsv"
#     resources:
#         mem_mb=3000,
#         disk_mb=33
#     threads:
#         1
#     shell:
#         "gurobi_cl OutputFlag=0 Method=2 Threads={threads} ResultFile={output.solution} ResultFile={output.json} ResultFile={params.ilp} {input} > {log}"

# rule process_solution:
#     message: "Processing Gurobi solution for '{output}'"
#     input:
#         solution="results/{scenario}/{model_run}.sol",
#         data="results/{scenario}/model_{model_run}/datapackage.json"
#     output: ["results/{{scenario}}/{{model_run, \d+}}/{}.csv".format(x) for x in RESULTS.index]
#     conda: "../envs/otoole.yaml"
#     log: "results/process_solution_{scenario}_{model_run}.log"
#     params:
#         folder="results/{scenario}/{model_run}"
#     shell:
#         "mkdir -p {params.folder} && otoole -v results gurobi csv {input.solution} {params.folder} --input_datapackage {input.data} 2> {log}"
