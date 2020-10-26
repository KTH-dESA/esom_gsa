# gui_workflow
The snakemake workflow for the Gulf UnderSea Interconnector feasibility study

## Installation

Install snakemake using conda into a new environment called `snakemake`:

```bash
conda install -c conda-forge mamba
mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal SALib pandas numpy graphviz glpk seaborn
```

Then, activate the environment using `source activate snakemake` on Mac and Linux, or `activate snakemake` on Windows.

Now install the other dependencies using pip:

```python
pip install -r requirements.txt
```

## Running the workflow

To run the workflow, using the command `snakemake --use-conda --cores all --resources mem_mb=16000`

## Plotting the workflow

To visualise the workflow, run the following rule: `snakemake plot_dag --use-conda  --cores 2`

## Folder structure

This repository follows the snakemake guidelines for reproducibility:

    ├── .gitignore
    ├── README.md
    ├── LICENSE.md
    ├── workflow
    │   ├── rules
    |   │   ├── module1.smk
    |   │   └── module2.smk
    │   ├── envs
    |   │   ├── tool1.yaml
    |   │   └── tool2.yaml
    │   ├── scripts
    |   │   ├── script1.py
    |   │   └── script2.R
    │   ├── notebooks
    |   │   ├── notebook1.py.ipynb
    |   │   └── notebook2.r.ipynb
    │   ├── report
    |   │   ├── plot1.rst
    |   │   └── plot2.rst
    |   └── Snakefile
    ├── config
    │   ├── config.yaml
    │   └── some-sheet.tsv
    ├── results
    └── resources