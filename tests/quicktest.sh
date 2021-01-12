#!/bin/bash

set -eEuo pipefail

eecho(){
    echo $@ 1>&2
}

CONDA_ACTIVATE="$(which activate)"
CONDA_ACTIVE=0

activate_conda(){
    echo $@
    set +ue
    source "$CONDA_ACTIVATE" --stack "$@"
    set -ue
    CONDA_ACTIVE=$((CONDA_ACTIVE+1))
    echo CONDA ENVIRONMENTS "$CONDA_ACTIVE"
}

deactivate_conda(){
    set +ue
    conda deactivate
    set -ue
    CONDA_ACTIVE=$((CONDA_ACTIVE-1))
    echo ACTIVE CONDA ENVIRONMENTS $CONDA_ACTIVE
}

stop_setup(){

    while [ $CONDA_ACTIVE -gt 0 ]; do
        deactivate_conda
    done

    unset CONDA_ACTIVE CONDA_ACTIVATE
    exit 0
}

main(){
    local PROXY=""
    snake_args=( "" ) 
	local covpipe_package_exists=0
	if conda list | grep covpipe ; then
		covpipe_package_exists=1
	fi  

    local skipconda=${covpipe_package_exists}
    while getopts "p:sn" flag; do
        case "$flag" in 
            p) PROXY=$OPTARG
               ;;
            s) skipconda=1
               ;;
            n) snake_args=( "${snake_args[@]}" "-n")
               ;;
        esac  
    done
    shift $(( $OPTIND - 1 ))  
    local SOURCE="$( realpath "${BASH_SOURCE[0]}" )"
    local SOURCE_DIR="$( dirname "$SOURCE" )"
    local TEST_DATA_DIR="${SOURCE_DIR}/TestData"
    local GIT_ROOT="$( realpath "${SOURCE_DIR}/..")"
    local SETUP_TESTENV="${SOURCE_DIR}/development_envs/setups/testDev.setup/setup_env.sh"
    local SETUP_PROGENV="${GIT_ROOT}/ncov_minipipe.conda.setup/setup_env.sh"
    local WORKDIR="${SOURCE_DIR}/QuickTest"
    mkdir -p "$WORKDIR"
    local TESTENVDIR=${TESTENVDIR-"${WORKDIR}/envs"}
    mkdir -p "$TESTENVDIR"
    # Build Conda Environments
    if [ $skipconda -eq 1 ]; then  
        echo Setting Up
        "$SETUP_PROGENV" "-s" "$TESTENVDIR"
    else 
        "$SETUP_PROGENV" "$TESTENVDIR"
        "$SETUP_TESTENV" "$TESTENVDIR"
    fi
    local PROGENV="${TESTENVDIR}/covpipe_environment"
    local TESTENV="${TESTENVDIR}/test_dev_env"
    cd "$WORKDIR"
    run_tests $PROXY
}

run_tests(){
    local http_proxy=${http_proxy:-""}
    local proxy=${1:-$http_proxy}
    local GEN_SCR="${TEST_DATA_DIR}/GEN_SAMPLES/gen_samples.sh"
    local VCF_SCR="${TEST_DATA_DIR}/VCF_FILES/get_vcf.sh"
    local REF_SCR="${TEST_DATA_DIR}/References/get_ref.sh"
    local ANN_SCR="${TEST_DATA_DIR}/References/get_gff.sh"
    activate_conda "$TESTENV"
    IFS=';' read -a READS  <<< "$("$GEN_SCR" Sample_7 gen)"
    IFS=';' read -a READS_TEST02  <<< "$("$GEN_SCR" Sample_6 gen)"
    local REF="$( http_proxy=$http_proxy $REF_SCR SarsCov2Wuhan1 gen)"
    local VCF="$($VCF_SCR VCF_ALL_VARS gen)"
    local adapter="${GIT_ROOT}/adapter.fasta"
    deactivate_conda
    local testid=test01
    mkdir -p $testid
    local testpath="$(realpath $testid)"
    local testsampleconf="${testpath}/sample.conf.yaml"
    local testconf="${testpath}/ncov_minipipe.config"
    printf "%s\n" \
            "output: \"$testpath\""\
            "samples: \"$testsampleconf\"" \
            "reference: \"$REF\"" \
            "annotation: 'false'" \
            "run_id:" \
            "amplification: 'false'" \
            "cutadapt_list_5prime: 'unknown'" \
            "cutadapt_list_3prime: 'unknown'" \
            "adapter_list: '$adapter'" \
            "kraken:" \
            "krakenTaxID: 2697049"\
            > "$testconf"
    printf "%s\n" \
            "Sample1:"\
            "    read1: \"${READS[0]}\""\
            "    read2: \"${READS[1]}\""\
            > "$testsampleconf"

    if [ $covpipe_package_exists -eq 0 ]; then  
    	activate_conda "$PROGENV"
	fi
    #echo ncov_minipipe --conf "$testconf" -o $testid \
    #    --conda_prefix "${TESTENVDIR}" --blame
    #ncov_minipipe --conf "$testconf" -o $testid \
    #    --conda_prefix "${TESTENVDIR}" --blame ${snake_args[@]}
    echo ncov_minipipe --input ${READS_TEST02[@]} -o test02 --reference $REF \
        --conda_prefix "${TESTENVDIR}" --wgs --blame
    ncov_minipipe --input ${READS_TEST02[@]} --reason -o test02 --reference $REF \
        --conda_prefix "${TESTENVDIR}" --wgs --blame  ${snake_args[@]} --debug
    
    ncov_minipipe --reason \
        --conf "$( echo test02/*ncov_minipipe.config | cut -d' ' -f1 )" \
        --conda_prefix "${TESTENVDIR}" --wgs --blame  ${snake_args[@]} --debug
    #  test3

    mkdir -p test03/input
    mkdir -p test03/input2
    cp -rf ${READS_TEST02[0]} test03/input/200102_20-01222_nk0_S1_L001_R1_001.fastq.gz 
    cp -rf ${READS_TEST02[1]} test03/input/200102_20-01222_nk0_S1_L001_R2_001.fastq.gz 
    cp -rf ${READS_TEST02[0]} test03/input/200102_20-01223_nk1_S1_L001_R1_001.fastq.gz 
    cp -rf ${READS_TEST02[1]} test03/input/200102_20-01223_nk1_S1_L001_R2_001.fastq.gz 
    cp -rf ${READS_TEST02[0]} test03/input/200102_20-01224_nk1_S1_L001_R1_001.fastq.gz 
    cp -rf ${READS_TEST02[1]} test03/input/200102_20-01224_nk1_S1_L001_R2_001.fastq.gz 
    cp -rf ${READS_TEST02[0]} test03/input2/200102_20-01223_nk1_S1_L001_R1_001.fastq.gz 
    cp -rf ${READS_TEST02[1]} test03/input2/200102_20-01223_nk1_S1_L001_R2_001.fastq.gz 

    ncov_minipipe --input test03/input --input test03/input2 -o test03/output --reference $REF \
        --conda_prefix "${TESTENVDIR}" --reason --run_id "run01-202001" --wgs \
        --blame ${snake_args[@]} --debug

    deactivate_conda
}


main $@
