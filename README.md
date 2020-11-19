# RNA FISH probes designer

The RNA FISH probes designer is a Nextflow workflow designed to simplify and automate designing oligonucleotide probes for RNA FISH single molecule imaging.
It is based primarily on OligoMiner python scripts (REF), as well as common bioinformatics tool.  

***

## Requirements
* Linux or OS operating system
* Bash &GreaterEqual; 3.2
* Java &GreaterEqual; 8


## Installation and configuration

1. Install conda (Anaconda3 or Miniconda3).
Check if conda is correctly installed with `conda info `

2. Navigate to a directory where you wish your workflow to be created.
Clone the RNA FISH probe designer repository from github with git
```
git clone https://github.com/episkadlo/RNAFISHProbeDesigner.git  
```
or manually download the zip file containing the necessary files by navigating to the repository listed above and clicking on Code > Downaload ZIP, then unpacking it into your current directory.

3. Move inside the RNAFISHProbeDesigner-main directory (containing main.nf), then clone or download files from the OligoMiner repository, either manually or with git
```
git clone https://github.com/beliveau-lab/OligoMiner.git
```  

4. Inside the RNAFISHProbeDesigner-main directory, install Nextflow workflow manager
```
curl -fsSL get.nextflow.io | bash
```

Your directory structure should now look as follows:

```
RNAFISHProbeDesigner-main
|-- genomes
|   `-- raw
|-- helperScripts
|   |-- buildJellyfishIndexes.sh
|   `-- customBed2Fasta.py
|-- env_python2.yml
|-- env.python3.yml
|-- main.nf
|-- nextflow
|-- nextflow.config
|-- OligoMiner-master
|-- README.md
|-- Results
|-- UPLOAD_FASTA_HERE
```

5. Inside the RNAFISHProbeDesigner-main directory, where the .yml files are, create two conda environments ProbeMakerEnv_python2 and ProbeMakerEnv_python3:
```
conda env create -f env_python2.yml
conda env create -f env_python3.yml
```
Note the full paths of the conda environments you just created. You can check those by running:
```
conda env lists
```
6. Open the nextflow.config file in any text editor. In the indicated places fill in the local paths to conda environments you just created.

7. To test if the workflow is ready and to see the help message with instructions of designing probes, run in the RNAFISHProbeDesigner-main directory:
```
./nextflow main.nf --help
```

## Prepare genomic indexes for specificy filtering

## Run the workflow to generate RNA FISH probes
### General instructions


### Examples

##### endogenous RNA target

##### exogenous RNA target


## References

## License
