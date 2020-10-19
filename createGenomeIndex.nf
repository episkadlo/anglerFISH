#!/usr/bin/env nextflow

params.genome = ""
genomePath = "$projectDir/genomes/raw/${params.genome}.fa"
params.indexName = ""
outPath = "$projectDir/genomes/indexes"

process createIndex {
  container = 'quay.io/biocontainers/hisat2:2.2.0--py36hf0b53f7_4'

  input:
  path genomePath

  script:
  """

  hisat2-build $genomePath ${outPath}/${params.indexName}
  mkdir ${outPath}/${params.indexName}
  mv ${outPath}/${params.indexName}.*.ht2 ${outPath}/${params.indexName}/
  chmod -R 777 ${outPath}/${params.indexName}/

  """

}
