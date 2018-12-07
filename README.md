# comfy_BED

comfy_BED is a program which produces BED files from LRG files. comfy_BED is designed to be integrated into bioinformatics pipelines. It runs from the command line.
To run comfy_BED, the user can either provide a local copy of a LRG file (offline mode), or have the LRG XML file fetched from the LRG website by its LRG ID (online mode). The user specifies a transcript (e.g. t1, t2) and genome build (GRCh37 or GRCh38). comfy_BED makes a BED file, and appends run information to a daily log file.


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


## Flags

-w OR -l: run comfy_BED in 'web' (-w) or 'local' (-l) mode. Provide an LRG ID for web mode. Provide a filepath for local mode.
-t: transcript option, e.g. t1.
-g: genome option, either GRCh37 or GRCh38.


## Usage examples

python comfy_BED.py -w LRG_1 -t t1 -g GRCh37
          Pulls LRG_1 from the web and outputs a BED file of transcript 1 in GRCh37
        
python comfy_BED.py -l ~/Documents/LRG_1.xml -t t1,t2 -g GRCh38
          Loads a local copy of LRG_1 and outputs a BED file in GRCh38 for each of 
          transcript 1 and transcript 2
