RESULTS = ['AnnualTechnologyEmission', 'TotalCapacityAnnual', 'RateOfActivity']
RUNS = ['low', 'central', 'high']

rule all:
    input: expand("results/{modelrun}/{results}.pdf", modelrun=RUNS, results=RESULTS)
    message: "Running pipeline to generate the files '{input}'"

rule generate_lp_file:
    input: data="data/simplicity_{model_run}.txt", model="model/osemosys_short.txt"
    output:
        protected("processed_data/{model_run}.lp.gz")
    log:
        "processed_data/glpsol_{model_run}.log"
    threads:
        1
    shell:
        "glpsol -m {input.model} -d {input.data} --wlp {output} --check --log {log}"

rule solve_lp:
    input:
        "processed_data/{model_run}.lp.gz"
    output:
        protected("processed_data/{model_run}.sol")
    log:
        "processed_data/cbc_{model_run}.log"
    threads:
        2
    shell:
        "cbc {input} -dualpivot pesteep -psi 1.0 -pertv 52 -duals solve -solu {output} > {log}"

rule process_solution:
    input:
        solution="processed_data/{model_run}.sol",
        data="data/simplicity_{model_run}.txt"
    output: expand("processed_data/{{model_run}}/{result}.csv", result=RESULTS)
    params:
        folder="processed_data/{model_run}"
    shell:
        "mkdir -p {params.folder} && otoole results cbc csv {input.solution} {params.folder} --input_datafile {input.data}"

# rule solve:
#     input: data="data/simplicity_{modelrun}.txt", model="model/osemosys_short.txt"
#     output: 
#         expand("processed_data/{{model_run}}/{result}.csv", result=RESULTS)
#     log: "processed_data/{modelrun}/glpsol.log"
#     conda: "env/osemosys.yaml"
#     shell:
#         "glpsol -d {input.data} -m {input.model} > {log} && echo {output}"

rule clean:
    shell:
        "rm -r processed_data/*"

rule clean_plots:
    shell:
        "rm -f processed_data/{modelrun}/*.pdf"

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