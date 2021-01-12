#!/bin/bash

set +eEuo pipefail

eecho(){
    echo $@ 1>&2
}

main(){
    local scrname="GET MODREF SCRIPT"
    local OLD_IFS=$IFS
    local QUERY="${1-".*"}"
    local runmode="${2-"info"}"
    local SOURCE="$( realpath "${BASH_SOURCE[0]}" )"
    local SOURCE_DIR="$( dirname "$SOURCE" )"
    local REF_SCR="$( realpath "$SOURCE_DIR/../References/get_ref.sh" )"
    local VCF_SCR="$( realpath "$SOURCE_DIR/../VCF_FILES/get_vcf.sh")"
    if [ -f "$QUERY" ]; then
        echo "$QUERY"
    elif [ -f "${SOURCE_DIR}/$( basename "$QUERY" )" ]; then
        echo "${SOURCE_DIR}/$( basename "$QUERY" )"
    else
        local numRes="$(grep "^$QUERY" "${SOURCE_DIR}/modref.csv" | wc -l)"
        if [ $numRes -eq 0 ]; then
            echo "${QUERY}_not_found"
            return
        fi
        eecho "Found $numRes fitting generator targets for query"
        while read line; do
            getref "$line" $runmode
        done < <(grep "^$QUERY" "${SOURCE_DIR}/modref.csv")
    fi

}

getref(){
    IFS=';' read -a args <<< $1
    local gen_scr="${SOURCE_DIR}/generate_modified_reference.sh"
    local numargs=${#args[@]}
    local runmode="${2-"info"}"
    local alt_runmode="info"
    [[ $runmode =~ g.* ]] && alt_runmode=$runmode
    local id=${args[0]}
    local vcf=${args[1]}
    local temp_vcf="$( "$VCF_SCR" "$vcf" "$alt_runmode" )"
    if [[ $temp_vcf =~ "not_found|_error_" ]]; then
        eecho -e "${scrname}: Attempt to find vcf failed\n ID=$id\n Return: $temp_vcf"
        echo "${vcf}_vcf_input_not_found"
        return
    fi
    vcf="$temp_vcf"
    local ref=${args[2]}
    local temp_ref="$( "$REF_SCR" "$ref" "$alt_runmode" )"
    if [[ $temp_ref =~ "not_found|_error_" ]]; then
        eecho -e "${scrname}: Attempt to find reference failed\n ID=$id\n Return: $temp_ref"
        echo "${ref}_ref_input_not_found"
        return
    fi
    ref="$temp_ref"
    local outfile="${SOURCE_DIR}/${id}.fa"

    case $runmode in
        i*)
            echo "$outfile"
            ;;
        c*)
            printf "%s\n" \
                "$gen_scr $vcf $ref $outfile"
            ;;
        r*|g*)
            "$gen_scr" "$vcf" "$ref" "$outfile"
            echo "$outfile"
            ;;
    esac
}


main $@
