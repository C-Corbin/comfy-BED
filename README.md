# comfy_BED

comfy_BED is a program which produces BED files from LRG files. comfy_BED is designed to be integrated into bioinformatics pipelines. It runs from the command line.  

To run comfy_BED, the user can either provide a local copy of a LRG file (offline mode), or have the LRG XML file fetched from the LRG website by its LRG ID (online mode). The user specifies a transcript (e.g. t1, t2) and genome build (GRCh37 or GRCh38). comfy_BED makes a BED file, and appends run information to a daily log file.


## comfy_BED setup

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


## Running comfy_BED

### Flags

`-h`: Show help text and exit

`-w` **OR** `-l`: Run comfy_BED in 'web' (-w) or 'local' (-l) mode. **Required** (only one of these options is required):
- Web mode `-w`: Pulls LRG data from the web. Provide an LRG ID, HGNC gene name or RefSeq/Ensembl ID. 
- Local mode `-l`: Loads LRG data from a local file. Provide a filepath to an LRG XML file.  

`-t`: Choice of transcript(s) to make BED file for. **Required**. Must match the transcript ID in the LRG, e.g. t1. Multiple transcripts can be processed by separating each transcript with a comma (no spaces), e.g. t1,t2.  

`-g`: Genome build option, either GRCh37 or GRCh38. **Optional**, defaults to GRCh37 if empty.  

### Usage examples

`python comfy_BED.py -w LRG_1 -t t1 -g GRCh37`  
Pulls LRG_1 from the web and outputs a BED file of transcript 1 in GRCh37
        
`python comfy_BED.py -l ~/Documents/LRG_1.xml -t t1,t2 -g GRCh38`  
Loads a local copy of LRG_1 and outputs a BED file in GRCh38 for each of transcript 1 and transcript 2

### Output

comfy_BED will output a BED file of the genomic co-ordinates of the LRG transcript selected. The output is in the standard BED format, with the exon number also included in the 4th column.  

The file will be named `<LRG_ID>_<transcript_ID>.bed`, where `<LRG_ID>` is the LRG number and `<transcript_ID>` is the transcript number, e.g. LRG_1_t1.bed.  

**Warning**: If the script is run when there is already a file of the same name in the directory, the old file will be overwritten.

### Logging

comfy_BED automatically produces a log file detailing the steps that it has carried out. This is useful if you require an audit trail in the case of an error, or to verify that comfy_BED has performed as expected.

The log file will be saved in the current directory with the name `<date>_comfy_BED.log`, where `<date>` is the current date. If comfy_BED is run multiple times in the same day, the logs from each run will be appended onto the same log file.

### Testing

comfy_BED is unit tested using the pytest package. To run the tests, navigate to the comfy_BED directory and run `pytest`.
