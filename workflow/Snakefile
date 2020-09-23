import pandas as pd
configfile: "config/config.yaml"

RESULTS = pd.read_csv(config["result_params"]).set_index('name')
RUNS = pd.read_csv(config["modelruns"]).set_index('name')

rule all:
    input: expand("results/{modelrun}/{results}.pdf", modelrun=RUNS.index, results=RESULTS.index)
    message: "Running pipeline to generate the files '{input}'"

rule generate_lp_file:
    message: "Generating the LP file for '{wildcards.model_run}'"
    input: data="data/simplicity_{model_run}.txt", model=config['model_file']
    output:
        protected("processed_data/{model_run}.lp.gz")
    log:
        "processed_data/glpsol_{model_run}.log"
    threads:
        1
    shell:
        "glpsol -m {input.model} -d {input.data} --wlp {output} --check --log {log}"

rule solve_lp:
    message: "Solving the LP for '{wildcards.model_run}'"
    input:
        "processed_data/{model_run}.lp.gz"
    output:
        protected("processed_data/{model_run}.sol")
    log:
        "processed_data/cbc_{model_run}.log"
    threads:
        1
    shell:
        "cbc {input} -dualpivot pesteep -psi 1.0 -pertv 52 -duals solve -solu {output} > {log}"

rule process_solution:
    message: "Processing CBC solution for '{wildcards.model_run}'"
    input:
        solution="processed_data/{model_run}.sol",
        data="data/simplicity_{model_run}.txt"
    output: expand("processed_data/{{model_run}}/{result}.csv", result=RESULTS.index)
    params:
        folder="processed_data/{model_run}"
    shell:
        "mkdir -p {params.folder} && otoole results cbc csv {input.solution} {params.folder} --input_datafile {input.data}"

rule clean:
    shell:
        "rm -rf results/* && rm -rf processed_data/*"

rule clean_plots:
    shell:
        "rm -f results/{modelrun}/*.pdf"

rule plot:
    input: "processed_data/{modelrun}/{result}.csv"
    output: "results/{modelrun}/{result}.pdf"
    conda: "env/plot.yaml"
    message: "Generating plot using '{input}' and writing to '{output}'"
    shell:
        "python scripts/plot_results.py {input} {output}"

rule make_dag:
    output: pipe("dag.txt")
    shell:
        "snakemake --dag > {output}"

rule plot_dag:
    input: "dag.txt"
    output: "dag.png"
    conda: "env/dag.yaml"
    shell:
        "dot -Tpng {input} > dag.png && open dag.png"
