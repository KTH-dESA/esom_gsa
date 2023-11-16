rule create_sample:
    message: "Creating sample for '{params.replicates}' trajectories and '{params.parameters}' parameters"
    params:
        replicates=config['replicates'],
        parameters=config['parameters']
    output: "modelruns/{scenario}/morris_sample.txt"
    conda: "../envs/sample.yaml"
    log: "results/log/create_{scenario}_sample.log"
    shell:
        "python workflow/scripts/create_sample.py {params.parameters} {output} {params.replicates}"

rule expand_sample:
    params:
        parameters=config['parameters']
    input: "modelruns/{scenario}/morris_sample.txt"
    output: expand("modelruns/{{scenario}}/model_{model_run}/sample_{model_run}.txt", model_run=MODELRUNS)
    conda: "../envs/sample.yaml"
    log: "results/log/expand_{scenario}_sample.log"
    shell:
        "python workflow/scripts/expand_sample.py {input} {params.parameters} {output}"

rule scale_sample:
    params:
        parameters=config['parameters']
    input:
        "modelruns/{scenario}/morris_sample.txt"
    output:
        "modelruns/{scenario}/morris_sample_scaled.txt"
    conda: "../envs/sample.yaml"
    log: "results/log/expand_{scenario}_sample.log"
    shell:
        "python workflow/scripts/scale_sample.py {input} {params.parameters} {output}"