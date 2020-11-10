sudo apt-get update
sudo apt-get install coinor-cbc unzip

# 1. Install environment for snakemake
wget 

conda install -c conda-forge mamba
mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal SALib pandas numpy graphviz glpk seaborn

conda activate snakemake
pip install -r requirements.txt

# 2. Get the OSeMOSYS model
wget https://github.com/OSeMOSYS/OSeMOSYS_GNU_MathProg/releases/download/v1.0/osemosys_gnu_mathprog_v1.0.zip
unzip osemosys_gnu_mathprog_v1.0.zip

# 3. Get the workflow
git clone git@github.com:ClimateCompatibleGrowth/gui_workflow.git

# 4. Get the model files
git clone git@github.com:ClimateCompatibleGrowth/gui_osemosys.git
cd gui_osemosys
git fetch origin
git checkout osemosys_v1