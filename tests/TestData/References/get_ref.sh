#!/bin/bash

set -eEuo pipefail

eecho(){
    echo $@ 1>&2
}

main(){
    local PROXY=""
    while getopts "p:" flag; do
        case "$flag" in 
            p) PROXY=$OPTARG;;
        esac  
    done
    local OLD_IFS=$IFS
    local QUERY="${1-".*"}"
    local runmode="${2-"info"}"
    local SOURCE="$( realpath "${BASH_SOURCE[0]}" )"
    local SOURCE_DIR="$( dirname "$SOURCE" )"
    if [ -f "$QUERY" ]; then
        echo "$QUERY"
    elif [ -f "${SOURCE_DIR}/$( basename "$QUERY" )" ]; then
        echo "${SOURCE_DIR}/$( basename "$QUERY" )"
    else
        local numRes="$(grep "^$QUERY" "${SOURCE_DIR}/references.csv" | wc -l)"
        if [ $numRes -eq 0 ]; then
            echo "${QUERY}_not_found"
            return
        fi
        eecho "Found $numRes fitting generator targets for query"
        while read line; do
            getref "$line" $runmode
        done < <(grep "^$QUERY" "${SOURCE_DIR}/references.csv")
    fi

}

getref(){
    #echo "$#:$*"
    IFS=';' read -a args <<< $1
    local numargs=${#args[@]}
    [ $numargs -eq 0 ] && return
    local id=${args[0]}
    local type=${args[1]}
    case "$type" in
        ncbi* )
           ncbi_get "$1"
           ;;
        f*)
           file_get "$1"
           ;;
        *)
           file_get "$1"
           ;;
    esac
}

file_get(){
    IFS=';' read -a args <<< $1
    local numargs=${#args[@]}
    local id=${args[0]}
    local type=${args[1]}
    local file=${args[2]}
    if [ -f "$file" ]; then
        echo "$file"
    elif [ -f "${SOURCE_DIR}/${file}" ]; then
        realpath "${SOURCE_DIR}/${file}"
    elif [ -f "${SOURCE_DIR}/$( basename "${file}")" ]; then
        realpath "${SOURCE_DIR}/$( basename "${file}")"
    else
        eecho -e "FILE NOT FOUND: (Not as absolute path and not relative to REF DIR)\n${file}"
    fi
}


ncbi_get(){
    IFS=';' read -a args <<< $1
    local ncbi_scr="${SOURCE_DIR}/search_ncbi_nuc.py"
    local numargs=${#args[@]}
    local id=${args[0]}
    local type=${args[1]}
    local query=${args[2]}
    local outfile="${SOURCE_DIR}/${id}.fna"

    case $runmode in
        c*|cmd)
            echo "$ncbi_scr \"$query\" -o \"${outfile}.gz\""
            ;;
        g*|r*)
            [ -f "$outfile" ] || [ -f "${outfile}.gz" ]  || \
                "$ncbi_scr" "$query" -o "${outfile}.gz" 1>&2
            [ -f  "$outfile" ] || gunzip "${outfile}.gz"
            echo $outfile
            ;;
        i*)
            echo $outfile
            ;;
   esac

}


main $@
