# Check Java installation
if type -p java; then
    echo "found java executable in PATH"
    _java=java
else
    echo "Please install Java"
    return
fi

if [[ "$_java" ]]; then
    version=$("$_java" -version 2>&1 | awk -F '"' '/version/ {print $2}')
    if [[ "$version" < "1.8" ]]; then
        echo "Please update Jave to version 8 / 1.8"
        return
    fi
fi

# Check git installation
if type -p git; then
    echo "found git executable in PATH"
else
    echo "Please install Git or follow the manual procedure described in the main README"
    return
fi

# Clone Oligominer and download nextflow
if [[ -d "./OligoMiner-master" ]]; then
    mv ./OligoMiner-master ./OligoMiner
fi
if [[ ! -d "./OligoMiner" ]]; then
    git clone https://github.com/beliveau-lab/OligoMiner.git
fi

# Install nextflow
if [[ ! -f "./nextflow" ]]; then
    if type -p wget; then
        wget -qO- https://get.nextflow.io | bash
    elif type -p curl; then
        curl -fsSL get.nextflow.io | bash
    else
        echo "Please install wget or curl to install nextflow"
        return
    fi
fi

# Check for conda and update nextflow.config
if type -p conda; then
    echo found conda executable in PATH

    if [[ ! -f "./ProbeMakerEnv_python2.yml" ]] || [[ ! -f "./ProbeMakerEnv_python3.yml" ]]; then
        echo "Environment.yml files not found"
        return
    fi

    # If exists, conda update
    # -eq 2 because both for name of environment and folder name
    if [[ $(conda env list | grep -o ProbeMakerEnv_python2 | wc -l) -eq 2 ]] &&
        [[ $(conda env list | grep -o ProbeMakerEnv_python3 | wc -l) -eq 2 ]]; then
        conda env update -f ProbeMakerEnv_python2.yml
        conda env update -f ProbeMakerEnv_python3.yml
    else
        conda env create -f ProbeMakerEnv_python2.yml
        conda env create -f ProbeMakerEnv_python3.yml
    fi

    # Fancy bash magic to replace the placeholders in nextflow.config
    # with the paths to the conda environments
    configFile="./nextflow.config"
    if [[ -f $configFile ]] &&
        [[ $(grep -oi "PATH/TO/ENV/" $configFile | wc -l) -eq 2 ]]; then
        condaEnvs=$(conda env list)
        prefix=" "
        suffix="ProbeMakerEnv_python2"
        regex="$prefix[a-zA-Z0-9\/\.]\+$suffix"
        condaPath=$(echo $condaEnvs | grep -oi $regex)
        condaPath=${condaPath#"$prefix"}
        condaPath=${condaPath%"$suffix"}
        echo "found conda path $condaPath"
        perl -i -pe "s+PATH/TO/ENV/+${condaPath}+g" $configFile
    else
        echo "nextflow.config file not found or not properly formated!"
        return
    fi
fi

# Check nextflow installation
./nextflow main.nf --help && echo Installation successful \:\) || echo Installation failed \:\(
