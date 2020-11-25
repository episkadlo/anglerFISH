# RNA FISH probes designer

The RNA FISH probes designer is a Nextflow [1] workflow designed to simplify and automate designing oligonucleotide probes for RNA FISH single molecule imaging.  
It is based primarily on OligoMiner [2], as well as several common bioinformatics tool to provide a streamlined manner of designing probes for endogenous and exogenous target RNAs.


## Contents
* [Installation](#installation)
* [Usage](#usage)
* [Example](#example)
* [Output Visualization](#output-visualization)
* [Workflow Scheme](#workflow-scheme)
* [References](#references)
* [Citation and Licence](#citation-and-license)


## Installation
This workflow requires a couple of things to be installed before we can get started. Make sure you have a **Unix** operating system running:
* Bash &GreaterEqual; 3.2
    * Both Linux and MacOS ship with bash, check your version with: `bash --version`
    * Installation on Linux: `sudo apt-get update && sudo apt-get install --only-upgrade Bash```
    * Installation on MacOS: `brew install bash`
* Java &GreaterEqual; 8
    * Check your installation with: `java -version`
    * Installation on Linux: `sudo apt install default-jdk`
    * Installation on MacOS: `brew install java11`
* Miniconda3 (or Anaconda3)
    * Check your installation with: `conda info`
    * Installation on Linux / MacOS: [https://docs.conda.io/en/latest/miniconda.html]()

If all requirements were met you can go ahead by installing and configuring this workflow:

1. Navigate to a directory where you wish your workflow to be created.
Clone the RNA FISH probe designer repository from GitHub with git
    ```
    git clone https://github.com/episkadlo/RNAFISHProbeDesigner.git  
    ```
    or manually download the zip file containing the necessary files by navigating to the repository listed above and clicking on Code > Download ZIP, then unpacking it into your current directory.

1. Move inside the RNAFISHProbeDesigner directory (containing main.nf), then clone or download files from the OligoMiner [2] repository, either manually or with git
    ```
    git clone https://github.com/beliveau-lab/OligoMiner.git
    ```  

1. Inside the RNAFISHProbeDesigner directory, install Nextflow workflow manager:
    ```
    wget -qO- https://get.nextflow.io | bash
    ```
    or
    ```
    curl -fsSL get.nextflow.io | bash
    ```
    An executable file nextflow will appear in the RNAFISHProbeDesigner directory.  
    Now, your directory structure should look as follows:
    ```
    RNAFISHProbeDesigner
    |-- genomes
    |   |-- indexes
    |   `-- raw
    |-- helperScripts
    |   |-- buildJellyfishIndexes.sh
    |   `-- customBed2Fasta.py
    |-- OligoMiner-master
    |-- Results
    |-- UPLOAD_FASTA_HERE
    |-- ProbeMakerEnv_python2.yml
    |-- ProbeMakerEnv_python3.yml
    |-- main.nf
    |-- nextflow
    |-- nextflow.config
    |-- OligoMiner-master
    |-- README.md
    ```  
1. Inside the RNAFISHProbeDesigner directory, where the .yml files are, create two conda environments ProbeMakerEnv_python2 and ProbeMakerEnv_python3:
    ```
    conda env create -f ProbeMakerEnv_python2.yml
    conda env create -f ProbeMakerEnv_python3.yml
    ```
    Note the full paths of the conda environments you just created. You can check those by running:
    ```
    conda env list
    ```
1. Open the nextflow.config file in any text editor. In the indicated places fill in the local paths to conda environments you just created.

1. To test if the workflow is ready and to see the help message with instructions of designing probes, run in the RNAFISHProbeDesigner directory:
    ```
    ./nextflow main.nf --help
    ```

## Usage
### Prepare genomic indexes for specificity filtering
After the first set of possible probes are identified by OligoMiner, there are several steps that ensure maximal specificity and performance of probes.

Firstly, the probes are mapped to the genome of the target organism. This step is important to filter out unspecific probes.
First level of specificity filtering is based on mapping the probes to the reference genome of the target organism with HISAT2 aligner software [3] and filtering out non-specifically binding probes. If the RNA target of the probes is endogenous (it occurs in the reference genome, for example tubulin mRNA), the filters will keep only the probes that bind only once to the genome, presumably probe only bind to its specific genomic target. If the probe is targeting an exogenous sequence (not found in a reference genome, for example GFP mRNA from GFP reporter inserted into a mouse genome), only the probes that do not bind anywhere in the target genome will pass the filtering ,to ensure that the probes should not bind to any other targets.  

The second level of specificity filtering is based on removing probes which have the sequence of common k-mers occuring in the target genome with Jellyfish software [4], which could be missed by HISAT2 filtering.

Both steps described above require preparing reference genome of the target organism in form of indexes/dictionaries for HISAT2 software and Jellyfish. The main.nf workflow is able to generate the necessary indexes given the reference genome:

1. Download the reference genome in fasta format. One excellent source of genomes is [here](https://hgdownload.soe.ucsc.edu/downloads.html) If the genome is in .fa.gz format, unpack it first by:
    ```bash
    gzip <genome>.fa.gz
    ```
    Place the extracted `<genome>.fa` file into `RNAFISHProbeDesigner/genomes/raw` directory.
2. From `RNAFISHProbeDesigner/` directory, run the pipeline in mode of creating HISAT2 and Jellyfish indexes by adding a flag `--createIndexes`, supplying the genome path or name and provide length range of oligos to create Jellyfish dictionaries for these lengths range:
    ```bash
    ./nextflow main.nf\
        --createIndexes\
        --genome_index <genome name>\
        --l <min length of probes>\
        --L <max length of probes>
    ```
    Note: The order of the flags is irrelevant.   
    The indexes will be automatically generated and placed in the `RNAFISHProbeDesigner/genomes/indexes/<genome name> folder.

3. (**Optional**) Adding extra indexes for Jellyfish with different k-mer lengths.  
If the HISAT2 index has already been prepared and you wish to add additional Jellyfish indexes (for example you already have indexes for probes 18-23-nucleotides-long, and now you wish to design probes of lengths 24-28-nucleotides-long), you have an option to create only the Jellyfish indexes for given k-mer range, and not HISAT2 index. To do it run together flags `--createIndexes --jf_only`.
    ```
    ./nextflow main.nf\
        --createIndexes\
        --jf_only\
        --genome_index <genome name>\
        --l <min length of probes>\
        --L <max length of probes>
    ```

### Run the workflow to generate RNA FISH probes

1. Copy the RNA target sequence in .fa format into the RNAFISHProbeDesigner/UPLOAD_FASTA_HERE/ directory.
2. From RNAFISHProbeDesigner/ directory, run the workflow to design probes specifying the detailed parameters of the probes:
```bash
./nextflow main.nf\
    --genome_index <genome name>\
    --name <name of .fa sequence, without .fa> <probes parameters>\
    --outputName <basename of output files>
```
Note: The order of the flags is irrelevant.    
Read about most important probe parameters you can specify in the `./nextflow. main.nf --help`
Typically, you will want to specify: <default values>
* **--l** minimal probe length (nucleotides) <18>
* **--L** maxinal probe length (nucleotides) <23>
* **--spacing** minimal spacing between probes (nucleotides) <2>
* **--t** minimal probe melting temperature (&deg;C) <47>
* **--T** maximal probe melting temperature (&deg;C) <60>
* **--F** formamide concentration in buffers (%) <10>
* **--mode** type of sequence: endogenous or exogenous (endo/exo) <endo>
* **--strand** strandness of the target RNA (-/+)(for endogenous only) <->

After a successful workflow execution, a zipped file will appear in RNAFISHProbeDesigner/Results directory, containing the following output files:
* .fa file with sequences of the probes that passed through all filters
* .fa file wih reverse complementary sequences of probes
* .bam and .bam.bai alignment files of probes against genome for endogenous sequences, or .doc alignment file of probes against target RNA sequence for exogenous sequences
* .txt file with melting temperatures of the final probes
* .txt log file with stored basic parameters of the run


## Example
We included two RNA sequences (`smc2_dm.fa`, `renilla.fa`) that you can use to test your local workflow.

1. Download the Drosophila melanogaster dm6.fa.gz genome (size: 43M) by clicking [here](https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz).
2. Unpack it and place it in the `RNAFISHProbeDesigner/genomes/raw/` directory.
3. Prepare indexes for HISAT2 and Jellyfish (for probes 18-23nt-long):
    ```bash
    ./nextflow main.nf --genome_index dm6 --l 18 --L 23 --createIndexes
    ```
4. Follow the following steps for either an [endogenous](#endogenous-rna-target) or [exogenous](#exogenous-rna-target) RNA target.

#### Endogenous RNA target
Design probes for an endogenous Structure Maintenance of chromosomes protein 2 SMC2 mRNA (RefSeq NM_137151) from *D.melanogaster*, located on minus (-) strand. The sequence is stored in `smc2_dm.fa` file in the repository.
To create 18-20nt long probes run:
```bash
./nextflow main.nf\
    --genome_index dm6\
    --name smc2_dm\
    --l 18\
    --L 20\
    --mode endo\
    --strand -\
    --outputName smc2_dm_test
```
The output name should be smc2_dm_test. Once run, unpack and check the log file. You should get 165 final probes (out of 174 probes before filtering).

#### Endogenous RNA target
Design probes for an exogenous RNA coding renilla insterted into *D.melanogaster* cells. Sequence is stored in renilla.fa file. Probes are supposed to be 18-23nt-long. The output name should be renilla_test
```bash
./nextflow main.nf\
    --genome_index dm6\
    --name renilla\
    --l 18\
    --L 23\
    --mode exo\
    --outputName renilla_test
```
After the run, unpack and check the log file. You should get 40 final probes (out of 44 probes before filtering).


## Output Visualization
For endogenous alignments, go to the [Integrative Genomic Viewer](https://igv.org/app/) and click on the *Genome > selector* and select **the same genome version** that was used to filter and map your probes. You can also upload your own genome, if a non-standard genome was used.

Then click on *Tracks > Local File* to upload the probe alignment files created by the pipeline. Select **both** the `.bam` file and `.bam.bai` file at the same time. Search for the name of your gene in the input box next to the magnifying glass :magifying-glass-tilted-right: to zoom in to the gene and view the aligned probes.


## Workflow Scheme
The workflow utilizes the following tools to automatically design RNA FISH probes:  
* Nextflow [1]
* OligoMiner [2]
* HISAT2 [3]
* Jellyfish [4]
* NUPACK [5-7]
* SAMtools [8]
* FASTX-Toolkit [9]
* EMBOSS/needle [10]   

![workflow_scheme](https://github.com/episkadlo/RNAFISHProbeDesigner/blob/main/workflow_scheme.jpg)


## References
[1]  P. Di Tommaso, M. Chatzou, E. W. Floden, P. P. Barja, E. Palumbo, and C. Notredame, “Nextflow enables reproducible computational workflows,” Nat. Biotechnol., vol. 35, no. 4, pp. 316–319, 2017.  
[2]  B. J. Beliveau et al., “OligoMiner provides a rapid, flexible environment for the design of genome-scale oligonucleotide in situ hybridization probes,” Proc. Natl. Acad. Sci. U. S. A., vol. 115, no. 10, pp. E2183–E2192, Mar. 2018.  
[3]  D. Kim, B. Langmead, and S. L. Salzberg, “HISAT: a fast spliced aligner with low memory requirements,” Nat. Methods, vol. 12, no. 4, pp. 357–360, 2015.   
[4]  G. Marçais and C. Kingsford, “A fast, lock-free approach for efficient parallel counting of occurrences of k-mers.,” Bioinformatics, vol. 27, no. 6, pp. 764–770, Mar. 2011.    
[5]  R. M. Dirks and N. A. Pierce, “A partition function algorithm for nucleic acid secondary structure including  pseudoknots.,” J. Comput. Chem., vol. 24, no. 13, pp. 1664–1677, Oct. 2003.   
[6]  R. M. Dirks and N. A. Pierce, “An algorithm for computing nucleic acid base-pairing probabilities including  pseudoknots.,” J. Comput. Chem., vol. 25, no. 10, pp. 1295–1304, Jul. 2004.  
[7]  R. M. Dirks, J. S. Bois, J. M. Schaeffer, E. Winfree, and N. A. Pierce, “Thermodynamic Analysis of Interacting Nucleic Acid Strands,” SIAM Rev., vol. 49, no. 1, pp. 65–88, Jan. 2007.   
[8]  H. Li et al., “The Sequence Alignment/Map format and SAMtools.,” Bioinformatics, vol. 25, no. 16, pp. 2078–2079, Aug. 2009.   
[9] Hannon, G.J. (2010) FASTX-Toolkit. http:\//hannonlab.cshl.edu fastx_toolkit  
[10] P. Rice, I. Longden, and A. Bleasby, “EMBOSS: the European Molecular Biology Open Software Suite.,” Trends Genet., vol. 16, no. 6, pp. 276–277, Jun. 2000.   

## Citation and License
This workflow comes without any warranty can be freely used under the [MIT license](opensourse.org/licenses/MIT).
