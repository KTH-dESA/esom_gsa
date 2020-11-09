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
        folder=directory("results/{scenario}/gcc_india_{model_run, \d+}"),
        dummy="results/{scenario}/gcc_india_{model_run, \d+}/datapackage.json",
    shell:
        "python workflow/scripts/create_modelrun.py {input.datapackage} {output.folder} {input.sample}"

rule generate_datafile:
    message: "Modifying data and generating datafile for '{output}'"
    input:
        datapackage="results/{scenario}/gcc_india_{model_run}/datapackage.json"
    output:
        "modelruns/{scenario}/gcc_india_{model_run}.txt"
    log:
        "results/otoole_{scenario}_{model_run}.log"
    shell:
        "otoole convert datapackage datafile {input} {output} 2> {log}"

rule generate_lp_file:
    message: "Generating the LP file for '{output}'"
    input:
        data=expand("modelruns/{{scenario}}/gcc_india_{{model_run}}.txt"),
        model=config['model_file']
    resources:
        mem_mb=5000,
        disk_mb=1300
    output:
        temporary("results/{scenario}/{model_run}.lp")
    log:
        "results/glpsol_{scenario}_{model_run}.log"
    threads:
        1
    shell:
        "glpsol -m {input.model} -d {input.data} --wlp {output} --check --log {log} 2> {log}"

rule solve_lp:
    message: "Solving the LP for '{output}'"
    input:
        "results/{scenario}/{model_run}.lp"
    output:
        protected("results/{scenario}/{model_run}.sol")
    log:
        "results/cbc_{scenario}_{model_run}.log"
    resources:
        mem_mb=3000,
        disk_mb=33
    threads:
        1
    shell:
        "cbc {input} solve -solu {output} > {log}"

rule process_solution:
    message: "Processing CBC solution for '{output}'"
    input:
        solution="results/{scenario}/{model_run}.sol",
        data="results/{scenario}/gcc_india_{model_run}/datapackage.json"
    output: ["results/{{scenario}}/{{model_run, \d+}}/{}.csv".format(x) for x in RESULTS.index]
    log: "results/process_solution_{scenario}_{model_run}.log"
    params:
        folder="results/{scenario}/{model_run}"
    shell:
        "mkdir -p {params.folder} && otoole -v results cbc csv {input.solution} {params.folder} --input_datapackage {input.data} 2> {log}"