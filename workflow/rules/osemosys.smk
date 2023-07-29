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

def csv_from_scenario(wildcards):
    return SCENARIOS.loc[int(wildcards.scenario), 'csv']

def config_from_scenario(wildcards):
    return SCENARIOS.loc[int(wildcards.scenario), 'config']

rule create_model_data:
    message: "Copying and modifying data for '{params.folder}'"
    input:
        csv=csv_from_scenario,
        sample="modelruns/{scenario}/model_{model_run}/sample_{model_run}.txt",
        config=config_from_scenario
    log: "results/log/copy_csv_{scenario}_{model_run}.log"
    params:
        folder=directory("results/{scenario}/model_{model_run}/data"),
    conda: "../envs/otoole.yaml"
    group: "gen_lp"
    output:
        csvs = expand("results/{{scenario}}/model_{{model_run}}/data/{x}.csv", x=INPUT_FILES)
    shell:
        "python workflow/scripts/create_modelrun.py {input.csv} {params.folder} {input.sample} {input.config}"

rule copy_otoole_config:
    message: "Copying otoole configuration file for '{params.folder}'"
    input:
        yaml=config_from_scenario,
        csvs=rules.create_model_data.output,
    log: "results/log/copy_csv_{scenario}_{model_run}.log"
    params:
        folder="results/{scenario}/model_{model_run}",
    output:
        yaml="results/{scenario}/model_{model_run}/config.yaml",
    log:
        "results/log/copy_csv_{scenario}_{model_run}.log"
    shell:
        """
        cp {input.yaml} {output.yaml}
        """

rule generate_datafile:
    message: "Generating datafile for '{output}'"
    input:
        csv = rules.create_model_data.output,
        config="results/{scenario}/model_{model_run}/config.yaml"
    output:
        datafile = temp("temp/{scenario}/model_{model_run}.txt")
    params:
        csv_dir = "results/{scenario}/model_{model_run}/data/"
    conda: "../envs/otoole.yaml"
    group: "gen_lp"
    log:
        "results/log/otoole_{scenario}_{model_run}.log"
    shell:
        "otoole -v convert csv datafile {params.csv_dir} {output.datafile} {input.config} 2> {log}"

rule modify_model_file:
    message: "Adding MODEX sets to model file"
    input:
        "temp/{scenario}/model_{model_run}.txt"
    output:
        temp("temp/{scenario}/model_{model_run}_modex.txt")
    group: "gen_lp"
    threads:
        1
    conda: "../envs/otoole.yaml"
    shell:
        "python workflow/scripts/add_modex_sets.py otoole {input} {output}"

rule generate_lp_file:
    message: "Generating the LP file for '{output}'"
    input:
        data="temp/{scenario}/model_{model_run}_modex.txt",
        model=config['model_file']
    resources:
        mem_mb=64000,
        disk_mb=16000,
        time=180
    output:
        temp(expand("temp/{{scenario}}/model_{{model_run}}.lp{zip_extension}", zip_extension=ZIP))
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
        "temp/{scenario}/model_{model_run}.lp.gz"
    group:
        "solve"
    output:
        temp("temp/{scenario}/model_{model_run}.lp")
    shell:
        "gunzip -fcq {input} > {output}"

rule solve_lp:
    message: "Solving the LP for '{output}' using {config[solver]}"
    input:
        "temp/{scenario}/model_{model_run}.lp"
    output:
        json="modelruns/{scenario}/model_{model_run}/{model_run}.json",
        solution=temp("temp/{scenario}/model_{model_run}.sol")
    log:
        "results/log/solver_{scenario}_{model_run}.log"
    params:
        ilp="results/{scenario}/model_{model_run}/solve.ilp",
        cplex="results/{scenario}/model_{model_run}/solve.cplex",
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
    input: "temp/{scenario}/model_{model_run}.sol"
    output: expand("temp/{{scenario}}/{{model_run}}.sol{zip_extension}", zip_extension=ZIP)
    shell: "gzip -fcq {input} > {output}"

rule unzip_solution:
    message: "Unzip solution file {input}"
    group: "results"
    input: "temp/{scenario}/model_{model_run}.sol.gz"
    output: temp("temp/{scenario}/model_{model_run}.sol")
    shell: "gunzip -fcq {input} > {output}"

rule transform_file:
    message: "Transforming CPLEX sol file '{input}'"
    group: 'results'
    input: rules.unzip_solution.output
    conda: "../envs/otoole.yaml"
    output:
        temp("temp/{scenario}/model_{model_run}_trans.sol")
    shell:
        "python workflow/scripts/transform_31072013.py {input} {output}"

rule sort_transformed_solution:
    message: "Sorting transformed CPLEX sol file '{input}'"
    group: 'results'
    input:
        "temp/{scenario}/model_{model_run}_trans.sol"
    output:
        temp("temp/{scenario}/model_{model_run}_sorted.sol")
    shell:
        "sort {input} > {output}"

rule process_solution:
    message: "Processing {config[solver]} solution for '{output}'"
    group: 'results'
    input:
        solution=solver_output,
        datafile="temp/{scenario}/model_{model_run}.txt",
        config="results/{scenario}/model_{model_run}/config.yaml",
    output: 
        expand("results/{{scenario}}/model_{{model_run}}/results/{csv}.csv", csv=OUTPUT_FILES)
    conda: "../envs/otoole.yaml"
    log: "results/log/process_solution_{scenario}_{model_run}.log"
    params:
        folder="results/{scenario}/model_{model_run}/results"
    shell: 
        """
        mkdir -p {params.folder}
        otoole -v results {config[solver]} csv {input.solution} {params.folder} --input_datafile {input.datafile} {input.config} &> {log}
        """

rule get_statistics:
    message: "Extract the {config[solver]} statistics from the sol file"
    input: rules.solve_lp.output.solution
    output: "modelruns/{scenario}/model_{model_run}/{model_run}.stats"
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
    input: expand("modelruns/{{scenario}}/model_{model_run}/{model_run}.stats",  model_run=MODELRUNS)
    output: "results/{scenario}/objective_{scenario}.csv"
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