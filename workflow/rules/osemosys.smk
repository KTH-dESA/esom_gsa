import json

wildcard_constraints:
    modelrun="\d+",
    scenarios="\d+"

rule add_export:
    message: "Adds matching export parameters"
    params:
        parameters=config['parameters']
    group: "gen_lp"
    input: "modelruns/{scenario}/{model_run}_sample_import.txt"
    output: "modelruns/{scenario}/{model_run}_sample_export.txt"
    log: "results/log/add_export_{scenario}_{model_run}_sample_export.log"
    shell:
        """
        grep -v -e 'ELGX' {input} > {output} 2> {log}
        grep -e 'ILGX' {input} | sed -e 's/ILGX/ELGX/' -e 's/-//g' >> {output} 2> {log}
        """

def datapackage_from_scenario(wildcards):
    return SCENARIOS.loc[int(wildcards.scenario), 'path']

rule copy_datapackage:
    message: "Copying and modifying datapackage for '{output.folder}'"
    input:
        datapackage=datapackage_from_scenario,
        sample="modelruns/{scenario}/{model_run}_sample_export.txt"
    log: "results/log/copy_datapackage_{scenario}_{model_run}.log"
    conda: "../envs/otoole.yaml"
    group: "gen_lp"
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
        temp("results/{scenario}/model_{model_run}.txt")
    conda: "../envs/otoole.yaml"
    group: "gen_lp"
    log:
        "results/log/otoole_{scenario}_{model_run}.log"
    shell:
        "otoole -v convert datapackage datafile {input} {output} 2> {log}"

rule modify_model_file:
    message: "Adding MODEX sets to model file"
    input:
        "results/{scenario}/model_{model_run}.txt"
    output:
        temp("results/{scenario}/model_{model_run}_modex.txt")
    group: "gen_lp"
    threads:
        1
    conda: "../envs/otoole.yaml"
    shell:
        "python workflow/scripts/add_modex_sets.py otoole {input} {output}"

rule generate_lp_file:
    message: "Generating the LP file for '{output}'"
    input:
        data="results/{scenario}/model_{model_run}_modex.txt",
        model=config['model_file']
    resources:
        mem_mb=64000,
        disk_mb=16000,
        time=180
    output:
        temp("{config[scratch]}/results/{scenario}/{model_run}.lp.gz")
    benchmark:
        "benchmarks/gen_lp/{scenario}_{model_run}.tsv"
    log:
        "results/log/glpsol_{scenario}_{model_run}.log"
    conda: "../envs/osemosys.yaml"
    group: "gen_lp"
    threads:
        1
    shell:
        "glpsol -m {input.model} -d {input.data} --wlp {output} --check > {log}"

rule unzip:
    message: "Unzipping LP file"
    input:
        "{config[scratch]}/results/{scenario}/{model_run}.lp.gz"
    group:
        "solve"
    output:
        temp("{config[scratch]}/results/{scenario}/{model_run}.lp")
    shell:
        "gunzip -fq {input}"

rule solve_lp:
    message: "Solving the LP for '{output}' using {config[solver]}"
    input:
        "{config[scratch]}/results/{scenario}/{model_run}.lp"
    output:
        json="results/{scenario}/{model_run}.json",
        solution="results/{scenario}/{model_run}.sol",
    log:
        "results/log/solver_{scenario}_{model_run}.log"
    params:
        ilp="results/{scenario}/{model_run}.ilp",
        cplex="results/{scenario}/{model_run}.cplex",
    benchmark:
        "benchmarks/solver/{scenario}_{model_run}.tsv"
    resources:
        mem_mb=30000,
        disk_mb=20000,
        time=1080
    group: "solve"
    threads:
        3
    shell:
        """
        if [ {config[solver]} = gurobi ]
        then
          gurobi_cl Method=2 Threads={threads} LogFile={log} LogToConsole=0 ScaleFlag=2 NumericFocus=3 ResultFile={output.solution} ResultFile={output.json} ResultFile={params.ilp} {input}
        elif [ {config[solver]} = cplex ]
        then
          echo "set threads {threads}"   > {params.cplex}
          echo "read {input}" 	     >> {params.cplex}
          echo "baropt"                  >> {params.cplex}
          echo "write {output.solution}" >> {params.cplex}
          echo "quit"                    >> {params.cplex}
        cplex < {params.cplex} > {log} && touch {output.json}
        else
          cbc {input} solve -sec 1500 -solu {output.solution} 2> {log} && touch {output.json}
        fi
        """

rule process_solution:
    message: "Processing {config[solver]} solution for '{output}'"
    input:
        solution="results/{scenario}/{model_run}.sol",
        data="results/{scenario}/model_{model_run}/datapackage.json"
    output: ["results/{{scenario}}/{{model_run, \d+}}/{}.csv".format(x) for x in RESULTS.index]
    conda: "../envs/otoole.yaml"
    log: "results/process_solution_{scenario}_{model_run}.log"
    params:
        folder="results/{scenario}/{model_run}"
    shell:
        "mkdir -p {params.folder} && otoole -v results {config[solver]} csv {input.solution} {params.folder} --input_datapackage {input.data} 2> {log}"