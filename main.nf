#!/usr/bin/env nextflow

params.name = ""
params.outputName = "${params.name}_output"
inFilePath = "$projectDir/UPLOAD_FASTA_HERE/${params.name}.fa"
params.l = 45
params.L = params.l
params.spacing = 2
params.F = 50
params.s = 390
params.g = 20
params.G = 80
params.t = 37
params.T = 50
params.hybrTemp = 37
params.genome_index = ''
genomeBasePath = "$projectDir/genomes/indexes"
genomeIndexPath = "$genomeBasePath/${params.genome_index}"
params.m = params.l
params.mode = "endo"
params.strand = "-"
params.overlapMode = "no"
params.createIndexes = false
params.rawGenomePath = "$projectDir/genomes/raw/${params.genome_index}.fa"
params.jf_only = false


process createHisat2index {
  container = 'quay.io/biocontainers/hisat2:2.2.0--py36hf0b53f7_4'

  input:
    path rawGenomePath from "$params.rawGenomePath"

  output:
    val 'done' into hiast2IndexProcess_created
    path "${params.genome_index}.*.ht2*" into hs2_indices_created, hs2_indices_created_duplicate

  when:
    params.createIndexes == true & params.jf_only == false

  script:

      """
      hisat2-build ${rawGenomePath} ${params.genome_index}
      mkdir -p "${projectDir}/genomes/indexes/${params.genome_index}/"
      cp *.ht2 ${projectDir}/genomes/indexes/${params.genome_index}/
      """
}


process collectHisat2index {

  input:
    path hisat2index from Channel.fromPath("$projectDir/genomes/indexes/${params.genome_index}/*.ht2").collect()

  output:
    val 'done' into hiast2IndexProcess_collected
    path hisat2index into hs2_indices_collected, hs2_indices_collected_duplicate

  when:
    params.createIndexes == false || params.jf_only == true

  script:

      """
      echo $hisat2index
      """
}


process createJellyfishIndex {
  container = 'jellyfish_ps:1.2'

  input:
    path rawGenomePath from "$params.rawGenomePath"
    path buildJellyfishIndexes from "$projectDir/helperScripts/buildJellyfishIndexes.sh"

  output:
    val 'done' into jellyfishIndexProcess_created
    path "*.jf" into jellyfish_indices_created, jellyfish_indices_created_duplicate

  when:
    params.createIndexes == true

  script:

      """
      bash ${buildJellyfishIndexes} ${params.l} ${params.L} ${params.genome_index} ${rawGenomePath}
      cp *.jf ${projectDir}/genomes/indexes/${params.genome_index}/
      """
}


process collectJellyfishindex {

  input:
    path jellyfishIndex from Channel.fromPath("$projectDir/genomes/indexes/${params.genome_index}/*.jf").collect()

  output:
    val 'done' into jellyfishIndexProcess_collected
    path jellyfishIndex into jellyfish_indices_collected, jellyfish_indices_collected_duplicate

  when:
    params.createIndexes == false

  script:

      """
      echo $jellyfishIndex
      """
}



process blockParse {
  container = 'testdocker3'

  input:
  val hisat2indexDone from hiast2IndexProcess_created.mix(hiast2IndexProcess_collected)
  val jellyfishIndexDone from jellyfishIndexProcess_created.mix(jellyfishIndexProcess_collected)
  path blockParseScript from "$projectDir/OligoMiner-master/blockParse.py"
  path inFile from inFilePath

  output:
  path "${params.name}_blockparse.fastq" into blockparseProcess, blockparseProcess_count

	script:
  if( params.overlapMode == "no" )
    """
    python ${blockParseScript} \
      -f ${inFile} \
      -o ${params.name}_blockparse \
      -l ${params.l} \
      -L ${params.L} \
      --Spacing ${params.spacing} \
      -F ${params.F} \
      -s ${params.s} \
      -g ${params.g} \
      -G ${params.G} \
      -t ${params.t} \
      -T ${params.T} \
      --header chr1:10000000-10050000
	  """

  else if( params.overlapMode == "yes" )
    """
    python ${blockParseScript} \
      -f ${inFile} \
      -o ${params.name}_blockparse \
      -l ${params.l} \
      -L ${params.L} \
      --Spacing ${params.spacing} \
      -F ${params.F} \
      -s ${params.s} \
      -g ${params.g} \
      -G ${params.G} \
      -t ${params.t} \
      -T ${params.T} \
      --header chr1:10000000-10050000 \
      --OverlapMode
    """

  else
    error "Invalid alignment mode: ${params.overlapMode}"

}


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


process initial_mapping {

  container = 'quay.io/biocontainers/hisat2:2.2.0--py36hf0b53f7_4'

	input:
  path "${params.name}_blockparse.fastq" from blockparseProcess
  path hs2_indices from hs2_indices_created.mix(hs2_indices_collected) .collect()


	output:
	path "${params.name}_hisat2.sam" into initialMappingProcess_endo, initialMappingProcess_exo

	script:
	"""
  hisat2 -x ${params.genome_index}  -U ${params.name}_blockparse.fastq -S ${params.name}_hisat2.sam
	"""
}


process filteringEndo {

  container = 'testdocker3'

	input:
  path outputCleanScript from "$projectDir/OligoMiner-master/outputClean.py"
  path "${params.name}_hisat2.sam" from initialMappingProcess_endo

	output:
	path "${params.name}_cleaned.bed" into outputCleanProcessEndo

  when:
  params.mode == "endo"

	script:

	"""
	python ${outputCleanScript} -F ${params.F}  -f ${params.name}_hisat2.sam -o ${params.name}_cleaned -u
	"""
}


process filteringExo {

  container = 'testdocker3'

	input:
  path outputCleanScript from "$projectDir/OligoMiner-master/outputClean.py"
  path "${params.name}_hisat2.sam" from initialMappingProcess_exo

	output:
	path "${params.name}_cleaned.bed" into outputCleanProcessExo

  when:
  params.mode == "exo"

	script:

	"""
	python ${outputCleanScript} -F ${params.F} -f ${params.name}_hisat2.sam -o ${params.name}_cleaned --zero
	"""
}


process checkStrand {

  container = 'testdocker3'

  input:
  path "${params.name}_cleaned.bed" from outputCleanProcessEndo.mix(outputCleanProcessExo)
  path probeRCScript from "$projectDir/OligoMiner-master/probeRC.py"

  output:
  path "${params.name}_strandChecked.bed" into checkStrandProcess


  script:

  if ((params.strand == "+" & params.mode == "endo") || params.mode == "exo")
    """
    python ${probeRCScript} -f ${params.name}_cleaned.bed -o ${params.name}_strandChecked
    """

  else
    """
    cp ${params.name}_cleaned.bed ${params.name}_strandChecked.bed
    """
}


process kmerFilter {

  container = 'jellyfish_ps:1.2'

	input:
  path "${params.name}_strandChecked.bed" from checkStrandProcess
  path kmerFilterScript from "$projectDir/OligoMiner-master/kmerFilter.py"
  // path jfDict from "${projectDir}/genomes/indexes/$params.genome_index/${params.genome_index}_${params.m}.jf"
  path jfDict from jellyfish_indices_created.mix(jellyfish_indices_collected) .collect()

	output:
  path "${params.name}_kmerFilter.bed" into kmerFilterProcess


  shell:
	'''
	for i in `seq !{params.l} !{params.L}`; do \
    python !{kmerFilterScript} -f !{params.name}_strandChecked.bed -m $i -j "!{projectDir}/genomes/indexes/!{params.genome_index}/!{params.genome_index}_$i.jf" -k 4 -o !{params.name}_kmerFilter_$i; \
  done
  awk '{print}' *_kmerFilter_*.bed | sort -u | grep . > !{params.name}_kmerFilter.bed
  '''

}


process structureCheck {

  // container = 'python2-nupack:1.3'
  conda '/home/ewa/anaconda3/envs/envForStructureCheck'

	input:
  path "${params.name}_kmerFilter.bed" from kmerFilterProcess
  path structureCheckScript from "$projectDir/OligoMiner-master/structureCheck.py"

	output:
  path "${params.name}_structureCheck.bed" into structureCheckProcessFastq, structureCheckProcessFasta

	script:

	"""
	python ${structureCheckScript} -f ${params.name}_kmerFilter.bed -o ${params.name}_structureCheck -t 0.05 -F ${params.F} --hybTemp ${params.hybrTemp}
	"""
}


process bed2fasta {

  conda '/home/ewa/anaconda3/envs/envForStructureCheck'

	input:
  path "${params.name}_structureCheck.bed" from structureCheckProcessFasta
  path bed2FastaScript from "$projectDir/helperScripts/customBed2Fasta.py"

	output:
	path "${params.name}_finalProbes.fasta" into bed2fastaProcessEndoZip, bed2fastaProcessExoZip, bed2fastaProcessExoAlign, bed2fastaProcessRevCompl, bed2fastaProcessTab

	script:

	"""
  python ${bed2FastaScript} -f "${params.name}_structureCheck.bed" -o "${params.name}_finalProbes"
	"""
}


process bed2fastq {

  container = 'testdocker3'


	input:
  path "${params.name}_structureCheck.bed" from structureCheckProcessFastq
  path bedToFastqScript from "$projectDir/OligoMiner-master/bedToFastq.py"

	output:
  path "${params.name}_finalProbes.fastq" into finalFastqProcessMapping, finalFastqProcessCount

  script:
	"""
	python ${bedToFastqScript} -f ${params.name}_structureCheck.bed  -o ${params.name}_finalProbes
	"""
}


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


process mappingFinalProbes {

  container = 'quay.io/biocontainers/hisat2:2.2.0--py36hf0b53f7_4'

	input:
  path "${params.name}_finalProbes.fastq" from finalFastqProcessMapping
  path hs2_indices_duplicate from hs2_indices_created_duplicate.mix(hs2_indices_collected_duplicate) .collect()

	output:
  path "${params.name}_finalProbes.sam" into mappingFinalProcess

  when:
  params.mode == "endo"

	script:
	"""
	hisat2 -x ${params.genome_index} -U ${params.name}_finalProbes.fastq -S ${params.name}_finalProbes.sam
	"""
}

process sam2bam {

  container = 'samtools_ps:1.0'

	input:
  path "${params.name}_finalProbes.sam" from mappingFinalProcess

	output:
  path "${params.name}_finalProbes.bam" into finalBamProcess

	script:
	"""
	samtools view -Sb ${params.name}_finalProbes.sam > ${params.name}_finalProbes.bam
	"""
}


process sortBam {

  container = 'samtools_ps:1.0'

	input:
  path "${params.name}_finalProbes.bam" from finalBamProcess

	output:
  path "${params.name}_finalProbes_sorted.bam" into alignment1ProcessSort, alignment1ProcessZip


	script:
	"""
	samtools sort ${params.name}_finalProbes.bam > ${params.name}_finalProbes_sorted.bam
	"""
}


process sortIndexBam {

  container = 'samtools_ps:1.0'

	input:
  path "${params.name}_finalProbes_sorted.bam" from alignment1ProcessSort

	output:
  path "${params.name}_finalProbes_sorted.bam.bai" into alignment2Process

	script:
	"""
  samtools index ${params.name}_finalProbes_sorted.bam
	"""
}



process revComplement {

  container = 'biocontainers/fastxtools:v0.0.14_cv2'

  input:
  path "${params.name}_finalProbes.fasta" from bed2fastaProcessRevCompl

  output:
  path "${params.name}_finalProbes_revComplement.fasta" into revComplEndo, revComplExo

  script:
  """
  fastx_reverse_complement -i ${params.name}_finalProbes.fasta -o ${params.name}_finalProbes_revComplement.fasta
  """
}

process probeTm {

  container = 'testdocker3'

  input:
  path "${params.name}_finalProbes.fasta" from bed2fastaProcessTab
  path probeTmScript from "$projectDir/OligoMiner-master/probeTm.py"

  output:
  path "${params.name}_finalProbes_Tm.txt" into probeTmEndo, probeTmExo

  script:
  """
  awk 'BEGIN{RS=">"}{print "#"\$1"\t"\$2;}' ${params.name}_finalProbes.fasta | tail -n+2 > ${params.name}_finalProbes.txt
  python ${probeTmScript} -f ${params.name}_finalProbes.txt -F ${params.F} -o ${params.name}_finalProbes_Tm
  """
}


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


  echo -e "Sequence name: ${params.name}\n" >> ${params.name}_log.txt
  echo -e "Date of run: \$run_date \n" >> ${params.name}_log.txt
  echo -e "\nParameters of the run:\n" >> "${params.name}_log.txt"
  echo -e "Type of sequence (endogenous / exogenous): ${params.mode}\n" >> ${params.name}_log.txt
  echo -e "Target organism genome: ${params.genome_index}\n"  >> ${params.name}_log.txt
  echo -e "Min probe length, nt: ${params.l}\n"  >> ${params.name}_log.txt
  echo -e "Max probe length, nt: ${params.L}\n"  >> ${params.name}_log.txt
  echo -e "Probes spacing, nt: ${params.spacing}\n"  >> ${params.name}_log.txt
  echo -e "Formamide concentration, %: ${params.F}\n"  >> ${params.name}_log.txt
  echo -e "min GC content, %: ${params.g}\n"  >> ${params.name}_log.txt
  echo -e "max GC content, %: ${params.G}\n"  >> ${params.name}_log.txt
  echo -e "min Tm, degC: ${params.t}\n"  >> ${params.name}_log.txt
  echo -e "max Tm, degC: ${params.T}\n" >> ${params.name}_log.txt
  echo -e "\n${initNprobes}\n" >> ${params.name}_log.txt
  echo -e "\n${finalNprobes}" >> ${params.name}_log.txt

  """

}



process alignExo{

  container = 'biocontainers/emboss:v6.6.0dfsg-7b1-deb_cv1'

  input:
  path "${params.name}_finalProbes.fasta" from bed2fastaProcessExoAlign

  output:
  path "${params.name}_alignment.doc" into alignExoProcess

  when:
  params.mode == "exo"

  script:
  """
  needle ${params.name}_finalProbes.fasta -sreverse2 -outfile ${params.name}_alignment.doc -gapopen 10 -gapextend 10
  """

}

process zipOutFilesEndo {

  input:
  path "${params.name}_finalProbes_sorted.bam" from alignment1ProcessZip
  path "${params.name}_finalProbes_sorted.bam.bai" from alignment2Process
  path "${params.name}_finalProbes.fasta" from bed2fastaProcessEndoZip
  path "${params.name}_finalProbes_revComplement.fasta" from revComplEndo
  path "${params.name}_finalProbes_Tm.txt" from probeTmEndo

  path "${params.name}_log.txt" from makeLogProcessEndo

  when:
  params.mode == "endo"

  script:
  """
  zip ${params.outputName}.zip ${params.name}_finalProbes_sorted.bam ${params.name}_finalProbes_sorted.bam.bai ${params.name}_finalProbes.fasta ${params.name}_finalProbes_revComplement.fasta ${params.name}_finalProbes_Tm.txt ${params.name}_log.txt
  cp ./${params.outputName}.zip ${projectDir}/Results/
  """
}



process zipOutFilesExo {
  input:
  path "${params.name}_finalProbes.fasta" from bed2fastaProcessExoZip
  path "${params.name}_alignment.doc" from alignExoProcess
  path "${params.name}_finalProbes_revComplement.fasta" from revComplExo
  path "${params.name}_finalProbes_Tm.txt" from probeTmExo
  path "${params.name}_log.txt" from makeLogProcessExo

  when:
  params.mode == "exo"

  script:
  """
  zip ${params.outputName}.zip ${params.name}_finalProbes.fasta ${params.name}_alignment.doc ${params.name}_finalProbes_revComplement.fasta ${params.name}_finalProbes_Tm.txt ${params.name}_log.txt
  cp ./${params.outputName}.zip ${projectDir}/Results/

  """
}
