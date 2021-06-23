#!/usr/bin/env nextflow

/*
======================================================
                    anglerFISH
======================================================
simple pipeline to design highly specific and
customizable RNA FISH probes against endo- and
exogenous target sequences

version 1.1
------------------------------------------------------

To be run with python 2.x
	createJellyfishIndex
  blockParse
  filteringEndo
  filteringExo
  checkStrand
  kmerFilter
  structureCheck
  bed2fasta
  bed2fastq
  sam2bam
  sortBam
  sortIndexBam
  revComplement
  probeTm
  alignExo

To be run with python 3.x
	createHisat2index
  initial_mapping
  mappingFinalProbes
*/

//collect and configure input parameters
params.name = ""
params.outputName = "${params.name}_output"
params.inFilePath = "$projectDir/UPLOAD_HERE/fasta/${params.name}.fa"
params.l = 18
params.L = 23
params.spacing = 2
params.F = 10
params.s = 390
params.g = 20
params.G = 80
params.t = 47
params.T = 60
params.hybrTemp = 37
params.genome_index = ""
params.genomeIndexPath = "$projectDir/UPLOAD_HERE/genome_index/"
params.mode = "endo"
params.strand = "-"
params.overlapMode = "no"
params.createIndexes = false
params.rawGenomePath = "$projectDir/UPLOAD_HERE/genome_raw/${params.genome_index}.fa"
params.jf_only = false
params.help = ""

// help message text
def helpMessage() {
  log.info """

  =========================================================================================================================
                                                      anglerFISH
  =========================================================================================================================
                  Simple nextflow pipeline to design highly specific and customizable RNA FISH probes.
  -------------------------------------------------------------------------------------------------------------------------

  Usage example:
    To design RNA FISH probes against sequence in example.fa of length 45-50 nucleotides, working in 50% formamide solution
    with check against specificity in mm10 genome:

    nextflow run main.nf --name 'example.fa' --genome_index mm10 --mode exo --l 45 --L 50 --F 50

  Manditory arguments:
      --name                    name of the input .fa file for which the probes are to be designer; if the file is not
                                located in the default UPLOAD_HERE/fasta folder, specify the path via --inFIlePath instead

      --genome_index            the name of the genome and it's indexes against which the probes will be aligned against

      --mode                    accepted values: endo / exo; use endo when target sequence occurs naturally in the genome
                                (i.e. probes against actin mRNA); use exo for targets that cannot be found in the genome
                                (i. e. probes against GFP)

      --strand                  strand of the target genome for the endogenous (endo) type of sequence; accepted values
                                are ["-", "+", "minus", "plus"], default: "-"

      --outputName              name of the output file; by default <name>_output

      --inFilePath              path to the .fa file; by default <projectDir>/UPLOAD_HERE/fasta/<name.fa>

      --genomeIndexPath         path to directory where the indexes for the genome are located;
                                by default <projectDir>/genomes/indexes/

  Specyfing parameters of the probes:
      --l                       minimal length of probe (nucleotides ), default: 18

      --L                       maximal length of probe (nucleotides), default: 23

      --spacing                 minimal allowed spacing between two probes, defualt: 2

      --F                       formamide concentration of formamide in the buffers (%), default: 10

      --s                       Na+ concentration in the buffers (mM), default: 390

      --g                       minimal GC content in a probe (%), default: 20

      --G                       maxinmal GC content in a probe (%), default: 80

      --t                       minimal Tm (melting temperature) of a probe (degC), default: 47

      --T                       minimal Tm (melting temperature) of a probe (degC), default: 60

      --hybrTemp                hybridization temperature (degC), default: 37

      --overlapMode             allow for overlaping probes (true/false), default: false

    Creating genome indexes for HISAT2 and Jellyfish:
      --createIndexes           trigger creating the HISAT2 and Jellyfish indexes of name specified
                                via --genome_index tag using the default path to the .fa file; the --genome_
                                index must match the fasta file name with the intup genome sequence; by deafult,
                                the createIndex will look for genome fasta file in <projectDir>/genomes/raw/,
                                unless it is overridden by specifying path of the genomic fasta
                                with --rawGenomePath tag; it will create Jellyfish dictionaries for lengths
                                specified via --l and --L (maximum and minimum probe length); it will
                                override existing the genome_index for HISAT2 and existing Jellyfish dictionaries;
                                with this flag, only specifying --genome_index (or --rawGenomePath if not in default
                                location) and minimal and maximal desired probes length (--l and --L) are required

      --jf_only                 only create Jellyfish indexes (dictionaries) for given --l and --L (maximum and minimum
                                probe length) and do NOT create HISAT2 idex; useful if --createIndex was
                                used before to create HISAT2 index, and the lengths of designed probes is now changing
                                and the Jellyfish indexes must be expanded to acommodate extra probe lengths

      --rawGenomePath           specify path of the fasta file with genomic sequence, which will be used for
                                creating HISAT2 and/or Jellyfish dictionaries, used with --createIndex flag;
                                default: <projectDir>/genomes/raw/<genome_index>
  """.stripIndent()
}

// show help message when the --help frag is called
if (params.help) {
    helpMessage()
    exit 0
}

// create HISAT2 index of genome using genome_index name matching the fasta file with genomic
// sequence, and optionally the path to the fasta file if not in the default location
process createHisat2index {

    input:
    path rawGenomePath from "$params.rawGenomePath"

    when:
    params.createIndexes == true && params.jf_only == false

    script:
    """
    hisat2-build ${rawGenomePath} ${params.genome_index}
    mkdir -p "${projectDir}/UPLOAD_HERE/genome_index/${params.genome_index}/"
    cp *.ht2 ${projectDir}/UPLOAD_HERE/genome_index/${params.genome_index}/
    """
}

// if the HISAT2 index of the genome index exists, collect the index files and create the channel for aignment
process collectHisat2index {

    input:
    path hisat2index from Channel.fromPath("${params.genomeIndexPath}/${params.genome_index}/*.ht2").collect()

    output:
    val "done" into hiast2IndexProcess_collected
    path hisat2index into hs2_indices_collected, hs2_indices_collected_duplicate

    when:
    params.createIndexes == false

    script:
    """
    echo $hisat2index
    """
}

// create Jellyfish index/dictionary of genome using genome_index name matching the fasta file with genomic
// sequence, and optionally the path to the fasta file if not in the default location
process createJellyfishIndex {

    input:
    path rawGenomePath from "$params.rawGenomePath"
    path buildJellyfishIndexes from "$projectDir/helperScripts/buildJellyfishIndexes.sh"

    when:
    params.createIndexes == true

    script:
    """
    bash ${buildJellyfishIndexes} ${params.l} ${params.L} ${params.genome_index} ${rawGenomePath}
    mkdir -p "${projectDir}/UPLOAD_HERE/genome_index/${params.genome_index}/"
    cp *.jf ${projectDir}/UPLOAD_HERE/genome_index/${params.genome_index}/
    """
}

// if the Jellyfish index/dictionary of the genome index exists for given length, collect the index files and create the channel for aignment
process collectJellyfishindex {

    input:
    path jellyfishIndex from Channel.fromPath("${params.genomeIndexPath}/${params.genome_index}/*.jf").collect()

    output:
    val "done" into jellyfishIndexProcess_collected
    path jellyfishIndex into jellyfish_indices_collected, jellyfish_indices_collected_duplicate

    when:
    params.createIndexes == false

    script:
    """
    echo $jellyfishIndex
    """
}

// call OligoMiner's script blockParse to create list of candidate probes with given parameters
process blockParse {

    input:
    val hisat2indexDone from hiast2IndexProcess_collected
    val jellyfishIndexDone from jellyfishIndexProcess_collected
    path blockParseScript from "$projectDir/OligoMiner/blockParse.py"
    path inFile from "${params.inFilePath}"

    output:
    path "${params.name}_blockparse.fastq" into blockparseProcess, blockparseProcess_count

    script:
    if( params.overlapMode == "no" )
      """
      python ${blockParseScript}\
        -f ${inFile}\
        -o ${params.name}_blockparse\
        -l ${params.l}\
        -L ${params.L}\
        --Spacing ${params.spacing}\
        -F ${params.F}\
        -s ${params.s}\
        -g ${params.g}\
        -G ${params.G}\
        -t ${params.t}\
        -T ${params.T}\
        --header chr1:10000000-10050000
        """

    else if( params.overlapMode == "yes" )
      """
      python ${blockParseScript}\
        -f ${inFile}\
        -o ${params.name}_blockparse\
        -l ${params.l}\
        -L ${params.L}\
        --Spacing ${params.spacing}\
        -F ${params.F}\
        -s ${params.s}\
        -g ${params.g}\
        -G ${params.G}\
        -t ${params.t}\
        -T ${params.T}\
        --header chr1:10000000-10050000\
        --OverlapMode
      """

    else
      exit 1, "Invalid alignment mode: ${params.overlapMode}"

}

// count the number of probes before filtering for reporting
process countInitialProbes {

    input:
    path "${params.name}_blockparse.fastq" from blockparseProcess_count

    output:
    stdout result into InitialNprobesProcess

    script:
    """
    initNprobes=\$(grep -c "@" ${params.name}_blockparse.fastq)
    echo "Number of probes found: \${initNprobes}"
    """
}

// use HISAT2 to align the candidate probes to genome specified in genome_index
process initial_mapping {

    input:
    path "${params.name}_blockparse.fastq" from blockparseProcess
    path hs2_indices from hs2_indices_collected.collect()


    output:
    path "${params.name}_hisat2.sam" into initialMappingProcess_endo, initialMappingProcess_exo

    script:
    """
    hisat2\
      -x ${params.genome_index}\
      -U ${params.name}_blockparse.fastq\
      -S ${params.name}_hisat2.sam
    """
}

// filter out unspecific cprobes using Oligominer's outputClean script for endogenous
// sequence, keeping only the probes that align exacly once (flag -u)
process filteringEndo {

    input:
    path outputCleanScript from "$projectDir/OligoMiner/outputClean.py"
    path "${params.name}_hisat2.sam" from initialMappingProcess_endo

    output:
    path "${params.name}_cleaned.bed" into outputCleanProcessEndo

    when:
    params.mode == "endo"

    script:
    """
    samtools view -q 60 ${params.name}_hisat2.sam  > ${params.name}_cleaned_temp.sam

    python ${outputCleanScript}\
      -F ${params.F}\
      -f ${params.name}_cleaned_temp.sam\
      -o ${params.name}_cleaned\
      -u
    """
}

// filter out unspecific cprobes using Oligominer's outputClean script for ezogenous
// sequence, keeping only the probes that do not align to the genome(flag --zero)
process filteringExo {

    input:
    path outputCleanScript from "$projectDir/OligoMiner/outputClean.py"
    path "${params.name}_hisat2.sam" from initialMappingProcess_exo

    output:
    path "${params.name}_cleaned.bed" into outputCleanProcessExo

    when:
    params.mode == "exo"

    script:
    """
    samtools view -f 4 ${params.name}_hisat2.sam  > ${params.name}_cleaned_temp.sam

    python ${outputCleanScript}\
      -F ${params.F}\
      -f ${params.name}_cleaned_temp.sam\
      -o ${params.name}_cleaned\
      --zero
    """
}

// check the standedness of the sequence, reverse-transcribe probes using OligoMiner's probeRC scripts
// if the original fasta sequence is on + strand for endogenous or if it is an exogenous sequences
process checkStrand {

  input:
    path "${params.name}_cleaned.bed" from outputCleanProcessEndo.mix(outputCleanProcessExo)
    path probeRCScript from "$projectDir/OligoMiner/probeRC.py"

  output:
    path "${params.name}_strandChecked.bed" into checkStrandProcess

  script:
    if (((params.strand == "+" || params.strand == "plus") & params.mode == "endo") || params.mode == "exo")
      """
      python ${probeRCScript}\
        -f ${params.name}_cleaned.bed\
        -o ${params.name}_strandChecked
      """

    else
      """
      cp ${params.name}_cleaned.bed ${params.name}_strandChecked.bed
      """
}

// using OligoMiner's kmerFileter script to filter out common k-mers from probes (utilizing jellyfish)
process kmerFilter {

    input:
    path "${params.name}_strandChecked.bed" from checkStrandProcess
    path kmerFilterScript from "$projectDir/OligoMiner/kmerFilter.py"
    path jfDict from jellyfish_indices_collected.collect()

    output:
    path "${params.name}_kmerFilter.bed" into kmerFilterProcess

    shell:
  	'''
  	for i in `seq !{params.l} !{params.L}`; do \
      python !{kmerFilterScript}\
        -f !{params.name}_strandChecked.bed\
        -m $i\
        -j "!{projectDir}/genomes/indexes/!{params.genome_index}/!{params.genome_index}_$i.jf"\
        -k 4 -o !{params.name}_kmerFilter_$i; \
    done
    awk '{print}' *_kmerFilter_*.bed | sort -u | grep . > !{params.name}_kmerFilter.bed
    '''
}

// using OligoMiner's structureCheck script to filter out
// probable secondary structures-forming probes (utilizing NUPACK)
process structureCheck {

    input:
    path "${params.name}_kmerFilter.bed" from kmerFilterProcess
    path structureCheckScript from "$projectDir/OligoMiner/structureCheck.py"

    output:
    path "${params.name}_structureCheck.bed" into structureCheckProcessFastq, structureCheckProcessFasta

    script:
    """
    python ${structureCheckScript}\
      -f ${params.name}_kmerFilter.bed\
      -o ${params.name}_structureCheck\
      -t 0.05\
      -F ${params.F}\
      --hybTemp ${params.hybrTemp}
    """
}

// get filtered probes in form of a fasta file
process bed2fasta {

    input:
    path "${params.name}_structureCheck.bed" from structureCheckProcessFasta
    path bed2FastaScript from "$projectDir/helperScripts/customBed2Fasta.py"

    output:
    path "${params.name}_finalProbes.fasta" into bed2fastaProcessEndoZip, bed2fastaProcessExoZip, bed2fastaProcessExoAlign, bed2fastaProcessRevCompl, bed2fastaProcessTm, bed2fastaProcessTab

    script:
    """
    python ${bed2FastaScript}\
      -f "${params.name}_structureCheck.bed"\
      -o "${params.name}_finalProbes"\
      --header "${params.name}"
    """
}

// converst fasta file with the final probes into a tab format for easier oligos ordering
process fasta2tab {

    input:
    path "${params.name}_finalProbes.fasta" from bed2fastaProcessTab

    output:
    path "${params.name}_oder.tab" into fasta2tabEndoZip, fasta2tabExoZip

    script:
    """
    seqkit fx2tab ${params.name}_finalProbes.fasta > ${params.name}_oder.tab
    """
}


// get filtered probes in form of a fastq file for visualization purposes in endogenous sequences
process bed2fastq {

    input:
    path "${params.name}_structureCheck.bed" from structureCheckProcessFastq
    path bedToFastqScript from "$projectDir/OligoMiner/bedToFastq.py"

    output:
    path "${params.name}_finalProbes.fastq" into finalFastqProcessMapping, finalFastqProcessCount

    script:
    """
    python ${bedToFastqScript}\
      -f ${params.name}_structureCheck.bed\
      -o ${params.name}_finalProbes
    """
}

// count the number of probes after all filtering steps
process countFinalProbes {

    input:
    path "${params.name}_finalProbes.fastq" from finalFastqProcessCount

    output:
    stdout result into FinalNprobesProcess

    script:
    """
    final_Nprobes=\$(grep -c "@" ${params.name}_finalProbes.fastq)
    echo "Number of probes that passed all the filters: \${final_Nprobes}"
    """
}

// map the filtered probes to the genome with HISAT2 (for later visualization purposes, endogenous sequences only)
process mappingFinalProbes {

    input:
    path "${params.name}_finalProbes.fastq" from finalFastqProcessMapping
    path hs2_indices_duplicate from hs2_indices_collected_duplicate.collect()

    output:
    path "${params.name}_finalProbes.sam" into mappingFinalProcess

    when:
    params.mode == "endo"

    script:
    """
    hisat2\
      -x ${params.genome_index}\
      -U ${params.name}_finalProbes.fastq\
      -S ${params.name}_finalProbes.sam
    """
}

// convert the sam file to bam, using samtools
process sam2bam {

    input:
    path "${params.name}_finalProbes.sam" from mappingFinalProcess

    output:
    path "${params.name}_finalProbes.bam" into finalBamProcess

    script:
    """
    samtools view -Sb ${params.name}_finalProbes.sam > ${params.name}_finalProbes.bam
    """
}

// sort the bam file, using samtools
process sortBam {

    input:
    path "${params.name}_finalProbes.bam" from finalBamProcess

    output:
    path "${params.name}_finalProbes_sorted.bam" into alignment1ProcessSort, alignment1ProcessZip

    script:
    """
    samtools sort ${params.name}_finalProbes.bam > ${params.name}_finalProbes_sorted.bam
    """
}

// index the bam file, using samtools
process sortIndexBam {

    input:
    path "${params.name}_finalProbes_sorted.bam" from alignment1ProcessSort

    output:
    path "${params.name}_finalProbes_sorted.bam.bai" into alignment2Process

    script:
    """
    samtools index ${params.name}_finalProbes_sorted.bam
    """
}

// create a fasta with reverse complementary sequences to the filtered probes using fastx
process revComplement {

    input:
    path "${params.name}_finalProbes.fasta" from bed2fastaProcessRevCompl

    output:
    path "${params.name}_finalProbes_revComplement.fasta" into revComplEndo, revComplExo

    script:
    """
    fastx_reverse_complement\
      -i ${params.name}_finalProbes.fasta\
      -o ${params.name}_finalProbes_revComplement.fasta
    """
}

// estimate the melting temperatures of the filtered probes using OligoMiner's probeTm script
process probeTm {

    input:
    path "${params.name}_finalProbes.fasta" from bed2fastaProcessTm
    path probeTmScript from "$projectDir/OligoMiner/probeTm.py"

    output:
    path "${params.name}_finalProbes_Tm.txt" into probeTmEndo, probeTmExo

    script:
    """
    awk 'BEGIN{RS=">"}{print "#"\$1"\t"\$2;}' ${params.name}_finalProbes.fasta |\
      tail -n+2 > ${params.name}_finalProbes.txt

    python ${probeTmScript}\
      -f ${params.name}_finalProbes.txt\
      -F ${params.F}\
      -o ${params.name}_finalProbes_Tm
    """
}

// create a log file holding key informations from the parameters used in the rnn
process makeLog {

    input:
    val initNprobes from InitialNprobesProcess
    val finalNprobes from FinalNprobesProcess

    output:
    stdout result into makeLogProcessPrint
    path "${params.name}_log.txt" into makeLogProcessEndo, makeLogProcessExo

    script:
    """
    touch ${params.name}_log.txt

    run_date=\$(date)

    echo -e "anglerFISH version 1.1" >> ${params.name}_log.txt
    echo -e "Sequence name: ${params.name}\n" >> ${params.name}_log.txt
    echo -e "Date of run: \$run_date \n" >> ${params.name}_log.txt
    echo -e "\nParameters of the run:\n" >> "${params.name}_log.txt"
    echo -e "Type of sequence (endogenous / exogenous): ${params.mode}\n" >> ${params.name}_log.txt
    echo -e "Target organism genome: ${params.genome_index}\n" >> ${params.name}_log.txt
    echo -e "Min probe length, nt: ${params.l}\n" >> ${params.name}_log.txt
    echo -e "Max probe length, nt: ${params.L}\n" >> ${params.name}_log.txt
    echo -e "Probes spacing, nt: ${params.spacing}\n" >> ${params.name}_log.txt
    echo -e "Formamide concentration, %: ${params.F}\n" >> ${params.name}_log.txt
    echo -e "min GC content, %: ${params.g}\n" >> ${params.name}_log.txt
    echo -e "max GC content, %: ${params.G}\n" >> ${params.name}_log.txt
    echo -e "min Tm, degC: ${params.t}\n" >> ${params.name}_log.txt
    echo -e "max Tm, degC: ${params.T}\n" >> ${params.name}_log.txt
    echo -e "\n${initNprobes}\n" >> ${params.name}_log.txt
    echo -e "\n${finalNprobes}" >> ${params.name}_log.txt
    """
}

// create local alignment between filtered probes and the input sequence, using needle from EMBOSS
process alignExo{

    input:
    path "${params.name}_finalProbes.fasta" from bed2fastaProcessExoAlign
    path inFile from "${params.inFilePath}"


    output:
    path "${params.name}_alignment.doc" into alignExoProcess

    when:
    params.mode == "exo"

    script:
    """
    needle ${inFile} ${params.name}_finalProbes.fasta\
      -sreverse2\
      -outfile ${params.name}_alignment.doc\
      -gapopen 10\
      -gapextend 10
    """
}

// create and zip file containign all the output file, for endogenous sequeence
process zipOutFilesEndo {

    publishDir "${projectDir}/UPLOAD_HERE/results/"

    output:
    path "./${params.outputName}.zip" into zip_endo

    input:
    path "${params.name}_finalProbes_sorted.bam" from alignment1ProcessZip
    path "${params.name}_finalProbes_sorted.bam.bai" from alignment2Process
    path "${params.name}_finalProbes.fasta" from bed2fastaProcessEndoZip
    path "${params.name}_oder.tab" from fasta2tabEndoZip
    path "${params.name}_finalProbes_revComplement.fasta" from revComplEndo
    path "${params.name}_finalProbes_Tm.txt" from probeTmEndo
    path "${params.name}_log.txt" from makeLogProcessEndo

    when:
    params.mode == "endo"

    script:
    """
    zip ${params.outputName}.zip\
      ${params.name}_finalProbes_sorted.bam\
      ${params.name}_finalProbes_sorted.bam.bai\
      ${params.name}_finalProbes.fasta\
      ${params.name}_oder.tab\
      ${params.name}_finalProbes_revComplement.fasta\
      ${params.name}_finalProbes_Tm.txt\
      ${params.name}_log.txt
    """
}

// create and zip file containign all the output file, for exogenous sequeence
process zipOutFilesExo {

    publishDir "${projectDir}/UPLOAD_HERE/results/"

    input:
    path "${params.name}_finalProbes.fasta" from bed2fastaProcessExoZip
    path "${params.name}_oder.tab" from fasta2tabExoZip
    path "${params.name}_alignment.doc" from alignExoProcess
    path "${params.name}_finalProbes_revComplement.fasta" from revComplExo
    path "${params.name}_finalProbes_Tm.txt" from probeTmExo
    path "${params.name}_log.txt" from makeLogProcessExo

    output:
    path "./${params.outputName}.zip" into zip_exo

    when:
    params.mode == "exo"

    script:
    """
    zip ${params.outputName}.zip\
      ${params.name}_finalProbes.fasta\
      ${params.name}_oder.tab\
      ${params.name}_alignment.doc\
      ${params.name}_finalProbes_revComplement.fasta\
      ${params.name}_finalProbes_Tm.txt\
      ${params.name}_log.txt
    """
}
