# comfy_BED

comfy_BED is a program which produces BED files from LRG files. comfy_BED is designed to be integrated into bioinformatics pipelines. It runs from the command line.  
To run comfy_BED, the user can either provide a local copy of a LRG file (offline mode), or have the LRG XML file fetched from the LRG website by its LRG ID (online mode). The user specifies a transcript (e.g. t1, t2) and genome build (GRCh37 or GRCh38). comfy_BED makes a BED file, and appends run information to a daily log file.


## Setup

### Download code

`git clone https://github.com/C-Corbin/comfy_BED.git`

### Set up either a conda or virtual environment:

comfy_BED has been designed to work with Python 2 only

#### Conda environment

- Set up conda environment: `conda env create -f requirements_conda.yaml`

- Activate the conda environment: `conda activate comfy-env`

- To exit conda environment: `conda deactivate`

#### Virtual environment

- Set up virtual env with python2 as Python interpreter: `virtualenv --python=/usr/bin/python2.7 comfy-env`

- Activate the virtual environment: `source comfy-env/bin/activate`

- Install requirements: `pip install -r requirements_pip.txt`

- To exit virtual environment: `deactivate`


## Flags

`-w` **OR** `-l`: Run comfy_BED in 'web' (-w) or 'local' (-l) mode. **Required** (only one of these options is required):
- Web mode `-w`: Pulls LRG data from the web. Provide an LRG ID, HGNS gene name or RefSeq/Ensembl ID. 
- Local mode `-l`: Loads LRG data from a local file. Provide a filepath to an LRG XML file.  

`-t`: Choice of transcript(s) to make BED file for. **Required**. Must match the transcript ID in the LRG, e.g. t1. Multiple transcripts can be processed by separating each transcript with a comma (no spaces), e.g. t1,t2.  

`-g`: Genome build option, either GRCh37 or GRCh38. **Optional**, defaults to GRCh37 if empty.  


## Usage examples

`python comfy_BED.py -w LRG_1 -t t1 -g GRCh37`  
Pulls LRG_1 from the web and outputs a BED file of transcript 1 in GRCh37
        
`python comfy_BED.py -l ~/Documents/LRG_1.xml -t t1,t2 -g GRCh38`  
Loads a local copy of LRG_1 and outputs a BED file in GRCh38 for each of transcript 1 and transcript 2
