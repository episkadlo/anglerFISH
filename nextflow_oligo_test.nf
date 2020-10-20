#!/usr/bin/env nextflow


###projectDir does not need to be explicitly defined, it's automatically created by nextflow (location of the pieplien)
process testBlockParse {

	input:
		path blockParseScript from "$projectDir/blockParse.py"
		path inFile from "$projectDir/example.fa"


	script:
	"""
		python $blockParseScript -f $inFile -o outTest
		cp outTest.fastq ${projectDir}
	"""

}
