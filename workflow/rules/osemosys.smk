wildcard_constraints:
    modelrun="\d+"

rule copy_datapackage:
    message: "Copying and modifying datapackage for '{output.folder}'"
    input:
        datapackage=config['datapackage'],
        sample="modelruns/{model_run, \d+}_sample.txt"
    output:
        folder=directory("processed_data/gcc_india_{model_run, \d+}"),
        dummy="processed_data/gcc_india_{model_run, \d+}/datapackage.json",
    shell:
        "python workflow/scripts/create_modelrun.py {input.datapackage} {output.folder} {input.sample}"

rule generate_datafile:
    message: "Modifying data and generating datafile for '{output}'"
    input:
        datapackage="processed_data/gcc_india_{model_run, \d+}/datapackage.json"
    output:
        "modelruns/gcc_india_{model_run}.txt"
    shell:
        "otoole convert datapackage datafile {input} {output}"

rule generate_lp_file:
    message: "Generating the LP file for '{output}'"
    input: 
        data=expand("modelruns/gcc_india_{{model_run}}.txt"),
        model=config['model_file']
    resources:
        mem_mb=8192
    output:
        protected("results/{model_run}.lp.gz")
    log:
        "results/glpsol_{model_run}.log"
    threads:
        1
    shell:
        "glpsol -m {input.model} -d {input.data} --wlp {output} --check --log {log}"

rule solve_lp:
    message: "Solving the LP for '{output}'"
    input:
        "results/{model_run}.lp.gz"
    output:
        protected("results/{model_run}.sol")
    log:
        "results/cbc_{model_run}.log"
    resources:
        mem_mb=8192
    threads:
        1
    shell:
        "cbc {input} -dualpivot pesteep -psi 1.0 -pertv 52 -duals solve -solu {output} > {log}"

rule process_solution:
    message: "Processing CBC solution for '{output}'"
    input:
        solution="results/{model_run}.sol",
        data="processed_data/gcc_india_{model_run}/datapackage.json"
    output: expand("results/{{model_run}}/{result}.csv", result=RESULTS.index)
    params:
        folder="results/{model_run}"
    shell:
        "mkdir -p {params.folder} && otoole results cbc csv {input.solution} {params.folder} --input_datapackage {input.data}"