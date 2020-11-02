#!/usr/bin/env nextflow
params.l = ""
params.L = ""
genomePath = "$projectDir/genomes/raw/${params.genome}.fa"
outPath = "$projectDir/genomes/indexes"

process imATest {
  container = 'jellyfish_ps:1.0'

  input:
    path buildJellyfishIndexes from "$projectDir/helperScripts/buildJellyfishIndexes.sh"
    path genomePath

  output:
    stdout results

  script:
    """
    bash ${buildJellyfishIndexes} ${params.l} ${params.L} ${params.genome} ${genomePath}
    mv *.jf ${outPath}/mm10_ch10/
    """
}


results.view { it }
