#!/bin/bash

set -eEuo pipefail

eecho(){
    echo $@ 1>&2
}

main(){
    local scrname="GET VCF SCRIPT"
    local OLD_IFS=$IFS
    local QUERY="${1-".*"}"
    local runmode="${2-"info"}"
    local SOURCE="$( realpath "${BASH_SOURCE[0]}" )"
    local SOURCE_DIR="$( dirname "$SOURCE" )"
    local REF_SCR="$( realpath "$SOURCE_DIR/../References/get_ref.sh" )"
    if [ -f "$QUERY" ]; then
        echo "$QUERY"
    elif [ -f "${SOURCE_DIR}/$( basename "$QUERY" )" ]; then
        echo "${SOURCE_DIR}/$( basename "$QUERY" )"
    else
        local numRes="$(grep "^$QUERY" "${SOURCE_DIR}/vcfs.csv" | wc -l)"
        if [ $numRes -eq 0 ]; then
            echo "${QUERY}_not_found"
            return
        fi
        eecho "Found $numRes fitting generator targets for query"
        while read line; do
            vcf_prog "$line" $runmode
        done < <(grep "^$QUERY" "${SOURCE_DIR}/vcfs.csv")
    fi
}


vcf_prog() {
    IFS=';' read -a args <<< $1
    local numargs=${#args[@]}
    local runmode="${2-"info"}"
    local alt_runmode="info"
    [[ $runmode =~ g.* ]] && alt_runmode=$runmode
    local id=${args[0]}
    local type=${args[1]}
    case $type in
        f*|file)
            get_vcf_file $1 $runmode
        ;;
        g*|gen)
            gen_vcf_file $1 $runmode
        ;;
        *)
            gen_vcf_file $1 $runmode
        ;;
    esac
}


gen_vcf_file(){
    local gen_scr="${SOURCE_DIR}/generate_vcf.py"
    IFS=';' read -a args <<< $1
    local numargs=${#args[@]}
    local runmode="${2-"info"}"
    local alt_runmode="info"
    [[ $runmode =~ g.* ]] && alt_runmode=$runmode
    local id=${args[0]}
    local type=${args[1]}
    local cmd=${args[2]}
    local mutfile=${args[3]}
    [ -f "$mutfile" ] || \
        mutfile="$(realpath "${SOURCE_DIR}/../mutations.yaml")" 

    local prelim_outfile="${SOURCE_DIR}/${id}.vcf"
    local outfile="${prelim_outfile}.gz"
    cmd=( $cmd )

    case $runmode in
        i*)
            echo "$outfile"
            ;;
        c*)
            printf "%s\n" \
                "$gen_scr $mutfile ${cmd[@]} -o $prelim_outfile" \
                "bgzip  $prelim_outfile"
            ;;
        r*|g*)

            if [ ! -f  $outfile ]; then
                "$gen_scr" "$mutfile" ${cmd[@]} -o "$prelim_outfile"
                bgzip "$prelim_outfile"
            fi
            echo "$outfile"
            ;;

    esac

}

get_vcf_file(){
   eecho File Reference Not implemented Yet
   echo "${id}_not_found"
}

main $@
