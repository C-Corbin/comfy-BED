# comfy_BED
A program to produce BED files from LRG files. 


## Usage
comfy_BED is designed to be integrated into bioinformatics pipelines. It runs from the command line.
The user can either provide a local copy of a LRG file (offline mode), or have the LRG XML file fetched from the LRG website by its LRG ID (online mode). The user specifies a transcript (e.g. t1, t2) and genome build (GRCh37 or GRCh38)

## Setup

### Download code

`git clone https://github.com/C-Corbin/comfy_BED.git`

### Set up either a conda or virtual environment:

comfy_BED has been designed to work with both Python 2 and Python 3

#### Conda environment

- Set up conda environment: `conda env create -f requirements_conda.yaml`

- Activate the conda environment: `conda activate comfy-env`

- To exit conda environment: `conda deactivate`

#### Virtual environment

- Set up virtual env with python2 as Python interpreter: `virtualenv comfy-env`

- Activate the virtual environment: `source comfy-env/bin/activate`

- Install requirements: `pip install -r requirements_pip.txt`

- To exit virtual environment: `deactivate`
