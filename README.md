# gui_workflow
The snakemake workflow for the Gulf UnderSea Interconnector feasibility study

## Installation

Install snakemake using conda into a new environment called `snakemake`:

```bash
conda install -c conda-forge mamba
mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal
```

Then, activate the environment using `source activate snakemake` on Mac and Linux, or `activate snakemake` on Windows.

Now install the other dependencies using pip:

```python
pip install -r requirements.txt
```

## Configure the workflow

Update the paths in the `config/config.yaml` file for `datapackage:` and `model_file:` keys. These should both hold relative paths to the combined GUI model `datapackage.json` and osemosys model file (e.g. `osemosys_fast.txt`).

Each of these paths can then point to the repositories outside of the gui_workflow. So deployment to an HPC will involve:

1. Copying a release package of OSeMOSYS to the server and unzipping
2. Cloning the gui_osemosys repository to the server
3. Cloning this repository to the server
4. Updating the config so that the paths point to the correct locations

## Running the workflow

To run the workflow, using the command `snakemake --use-conda --cores 4`

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