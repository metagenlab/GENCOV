#!/bin/bash

set -eEuo pipefail
RET_VAL=0


shutdown(){
    local ret=$RET_VAL
    unset RET_VAL
    exit $ret
}


error_echo(){
    RET_VAL=1
    echo -e $@
}

usage(){
    local error_state=${1-"0"}
    [[ "$error_state" != 0 ]] && error_echo "$error_state"
    printf "%s\n" \
        "generate_modified_reference [variation.vcf] [reference.fa] [optional output file]"
    [[ "$error_state" != 0 ]] && shutdown
}


bgzip_step(){
    local varfile="$1"
    local comp="${varfile##*vcf}"
    [[  $comp == ".gz" ]] || [ -f "${varfile}.gz" ] || bgzip "$varfile"
}

tabix_step(){
    local varfile="$1"
    [ -f "${varfile}.tbi" ] || tabix "${varfile}"
}

generate_step(){
    local varfile="$1"
    local reffile="$2"
    local outfile="$3"
    if [[ "$outfile" == "0" ]]; then
        bcftools consensus -f ${reffile} "${varfile}"
    else
        bcftools consensus -f ${reffile} "${varfile}" > "$outfile"
    fi
}

main(){
    local num_args=$#
    [ $num_args -lt 2 ] &&  usage "ERROR: Need to provide input files!\n"
    local varfile="$1"
    local reffile="$2"
    bgzip_step "$varfile"
    varfile="${varfile%%vcf*}vcf.gz"
    tabix_step "$varfile"
    if [ $num_args -eq 3 ];  then local outfile="$3"; fi
    generate_step "$varfile" "$reffile" "${outfile-0}"

    shutdown
}



main $@
