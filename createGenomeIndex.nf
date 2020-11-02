#!/usr/bin/env nextflow

params.l = ""
params.L = ""
params.genome = ""

genomePath = "$projectDir/genomes/raw/${params.genome}.fa"
params.indexName = "$params.genome"
outPath = "$projectDir/genomes/indexes/$params.indexName/"


process createOutDir {

  output:
    val created into dirCreated1, dirCreated2

  script:
    created=true
    """
    mkdir -p ${outPath}
    """
}

process buildHisat2Index {
  container = 'quay.io/biocontainers/hisat2:2.2.0--py36hf0b53f7_4'

  input:
    val created from dirCreated1
    path genomePath

  script:
    """
    hisat2-build ${genomePath} ${outPath}/${params.indexName}
    chmod -R 777 ${outPath}
    """
}



process buildJellyfishIndex {
  container = 'jellyfish_ps:1.0'


  input:
    val created from dirCreated2
    path buildJellyfishIndexes from "$projectDir/helperScripts/buildJellyfishIndexes.sh"
    path genomePath

  output:
    stdout results

  script:
    """
    bash ${buildJellyfishIndexes} ${params.l} ${params.L} ${params.genome} ${genomePath}
    mv *.jf ${outPath}
    """


}
