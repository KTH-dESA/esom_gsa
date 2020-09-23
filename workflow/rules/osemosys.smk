rule generate_lp_file:
    message: "Generating the LP file for '{wildcards.model_run}'"
    input: data="resources/simplicity_{model_run}.txt", model=config['model_file']
    output:
        protected("results/{model_run}.lp.gz")
    log:
        "results/glpsol_{model_run}.log"
    threads:
        1
    shell:
        "glpsol -m {input.model} -d {input.data} --wlp {output} --check --log {log}"

rule solve_lp:
    message: "Solving the LP for '{wildcards.model_run}'"
    input:
        "results/{model_run}.lp.gz"
    output:
        protected("results/{model_run}.sol")
    log:
        "results/cbc_{model_run}.log"
    threads:
        1
    shell:
        "cbc {input} -dualpivot pesteep -psi 1.0 -pertv 52 -duals solve -solu {output} > {log}"

rule process_solution:
    message: "Processing CBC solution for '{wildcards.model_run}'"
    input:
        solution="results/{model_run}.sol",
        data="resources/simplicity_{model_run}.txt"
    output: expand("results/{{model_run}}/{result}.csv", result=RESULTS.index)
    params:
        folder="results/{model_run}"
    shell:
        "mkdir -p {params.folder} && otoole results cbc csv {input.solution} {params.folder} --input_datafile {input.data}"