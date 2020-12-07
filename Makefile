install:
	bash helperScripts/install.sh

clean:
	rm -rf ./work
	rm -rf .nextflow/
	rm .nextflow.l*

uninstall:
	conda env remove --name ProbeMakerEnv_python2
	conda env remove --name ProbeMakerEnv_python3
	rm -rf ./work
	rm -rf .nextflow/
	rm .nextflow.l*
