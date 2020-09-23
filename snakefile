RESULTS = ['tid_demand', 'total_annual_capacity']
RUNS = ['low', 'central', 'high', '2degrees', '1.5degrees']

rule all:
	input: expand("processed_data/{modelrun}/{x}.pdf", x=RESULTS, modelrun=RUNS)
	message: "Running pipeline to generate the files '{input}'"

rule prepare_model:
	output: "processed_data/{modelrun}/osemosys.txt"
	params: modelrun="{modelrun}"
	shell:
		"bash scripts/osemosys_run.sh {params.modelrun}"

rule solve:
	input: data="data/simplicity.txt", model="processed_data/{modelrun}/osemosys.txt"
	output: results="processed_data/{modelrun}/results.sol", default="processed_data/{modelrun}/SelectedResults.csv"
	log: "processed_data/{modelrun}/glpsol.log"
	conda: "env/osemosys.yaml"
	shell:
		"glpsol -d {input.data} -m {input.model} -o {output.results} > {log}"

rule clean:
	shell:
		"rm -f processed_data/*/*.pdf processed_data/*/*.sol processed_data/*/*.csv *.png processed_data/*.csv processed_data/*.sol processed_data/*.log"

rule clean_plots:
	shell:
		"rm -f processed_data/{modelrun}/*.pdf"

rule extract_tid_demand:
	input: "processed_data/{modelrun}/SelectedResults.csv"
	output: "processed_data/{modelrun}/tid_demand.csv"
	shell:
		"head -n 33 {input} | tail -n 22 > {output}"

rule extract_total_annual_capacity:
	input: "processed_data/{modelrun}/SelectedResults.csv"
	output: "processed_data/{modelrun}/total_annual_capacity.csv"
	shell:
		"head -n 326 {input} | tail -n 29 > {output}"

rule plot:
	input: "processed_data/{modelrun}/{result}.csv"
	output: "processed_data/{modelrun}/{result}.pdf"
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