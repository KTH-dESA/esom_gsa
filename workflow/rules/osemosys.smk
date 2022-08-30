import json

wildcard_constraints:
    modelrun="\d+",
    scenarios="\d+"

ruleorder: unzip_solution > solve_lp

def solver_output(wildcards):
    if config['solver'] == 'cplex':
        return rules.sort_transformed_solution.output
    else: 
        return rules.unzip_solution.output

rule add_export:
    message: "Adds matching export parameters for TEMBA"
    params:
        parameters=config['parameters']
    group: "gen_lp"
    input: "modelruns/{scenario}/{model_run}_sample_import.txt"
    output: "modelruns/{scenario}/{model_run}_sample_export.txt"
    log: "results/log/add_export_{scenario}_{model_run}_sample_export.log"
    shell:
        """
        set +o pipefail
        grep -v -e 'ILGX' {input} > {output} 2> {log}
        grep -e 'ELGX' {input} | sed -e 's/ELGX/ILGX/' -e 's/-//g' >> {output} 2> {log}
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
        expand("{scratch}/results/{{scenario}}/{{model_run}}.lp.gz", scratch=SCRATCH)
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
        expand("{scratch}/results/{{scenario}}/{{model_run}}.lp.gz", scratch=SCRATCH)
    group:
        "solve"
    output:
        temp(expand("{scratch}/results/{{scenario}}/{{model_run}}.lp", scratch=SCRATCH))
    shell:
        "gunzip -fcq {input} > {output}"

rule solve_lp:
    message: "Solving the LP for '{output}' using {config[solver]}"
    input:
        expand("{scratch}/results/{{scenario}}/{{model_run}}.lp", scratch=SCRATCH)
    output:
        json="results/{scenario}/{model_run}.json",
        solution=temp(expand("{scratch}/results/{{scenario}}/{{model_run}}.sol", scratch=SCRATCH))
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
        time=720
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
          echo "set timelimit 43200"     >> {params.cplex}
          echo "read {input}" 	         >> {params.cplex}
          echo "baropt"                  >> {params.cplex}
          echo "write {output.solution}" >> {params.cplex}
          echo "quit"                    >> {params.cplex}
        cplex < {params.cplex} > {log} && touch {output.json}
        else
          cbc {input} solve -sec 1500 -solu {output.solution} 2> {log} && touch {output.json}
        fi
        """

rule zip_solution:
    message: "Zip up solution file {input}"
    group: "solve"
    input: expand("{scratch}/results/{{scenario}}/{{model_run}}.sol", scratch=SCRATCH)
    output: expand("{scratch}/results/{{scenario}}/{{model_run}}.sol.gz", scratch=SCRATCH)
    shell: "gzip -fcq {input} > {output}"

rule unzip_solution:
    message: "Unzip solution file {input}"
    group: "results"
    input: expand("{scratch}/results/{{scenario}}/{{model_run}}.sol.gz", scratch=SCRATCH)
    output: temp(expand("{scratch}/results/{{scenario}}/{{model_run}}.sol", scratch=SCRATCH))
    shell: "gunzip -fcq {input} > {output}"

rule transform_file:
    message: "Transforming CPLEX sol file '{input}'"
    group: 'results'
    input: rules.unzip_solution.output
    conda: "../envs/otoole.yaml"
    output:
        temp(expand("{scratch}/results/{{scenario}}/{{model_run}}_trans.sol", scratch=SCRATCH))
    shell:
        "python workflow/scripts/transform_31072013.py {input} {output}"

rule sort_transformed_solution:
    message: "Sorting transformed CPLEX sol file '{input}'"
    group: 'results'
    input:
        expand("{scratch}/results/{{scenario}}/{{model_run}}_trans.sol", scratch=SCRATCH)
    output:
        temp(expand("{scratch}/results/{{scenario}}/{{model_run}}_sorted.sol", scratch=SCRATCH))
    shell:
        "sort {input} > {output}"

rule process_solution:
    message: "Processing {config[solver]} solution for '{output}'"
    group: 'results'
    input:
        solution=solver_output,
        data="results/{scenario}/model_{model_run}/datapackage.json"
    output: ["results/{{scenario}}/{{model_run, \d+}}/{}.csv".format(x) for x in RESULTS.index]
    conda: "../envs/otoole.yaml"
    log: "results/log/process_solution_{scenario}_{model_run, \d+}.log"
    params:
        folder="results/{scenario}/{model_run, \d+}"
    shell: """
        mkdir -p {params.folder}
        otoole -v results {config[solver]} csv {input.solution} {params.folder} --input_datapackage {input.data} &> {log}
        """

rule get_statistics:
    message: "Extract the {config[solver]} statistics from the sol file"
    input: rules.solve_lp.output.solution
    output: "results/{scenario}/{model_run}.stats"
    group: "solve"
    shell: 
        """
        if [ {config[solver]} = cplex ]
        then
          head -n 27 {input} | tail -n 25 > {output}
        else
          head -n 1 {input} > {output}
        fi
        """

rule get_objective_value:
    input: expand("results/{{scenario}}/{model_run}.stats", model_run=MODELRUNS)
    output: "results/objective_{scenario}.csv"
    shell:
        """
        echo "FILE,OBJECTIVE,STATUS" > {output}
        if [ {config[solver]} = cplex ]
        then
          for FILE in {input}
          do
          OBJ=$(head $FILE | grep -e 'objectiveValue' | cut -f 2 -d '=')
          STATUS=$(head $FILE | grep -e 'solutionStatusString' | cut -f 2 -d '=')
          JOB=$(echo $FILE | cut -f 3 -d '/' | cut -f 1 -d '.')
          echo "$JOB,$OBJ,$STATUS" >> {output}
          done
        elif [ {config[solver]} = cbc ]
        then
          i=0
          for FILE in {input}
          do
          OBJ=$(head $FILE | cut -f 5 -d ' ')
          STATUS=$(head $FILE | cut -f 1 -d ' ')
          JOB=$FILE
          echo "$JOB,$OBJ,$STATUS" >> {output}
          ((i=i+1))
          done
        else
          echo "To be done"
        fi
        """