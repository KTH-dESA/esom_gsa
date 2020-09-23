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

## Running the workflow

To run the workflow, using the command `snakemake --use-conda --cores 4`

## Plotting the workflow

To visualise the workflow, run the following rule: `snakemake plot_dag --use-conda  --cores 2`
