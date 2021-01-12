#!/bin/bash

set -eEuo pipefail

eecho(){
    echo $@ 1>&2
}

main(){
    local OLD_IFS=$IFS
    local QUERY="${1-".*"}"
    local runmode="${2-"info"}"
    local SOURCE="$( realpath "${BASH_SOURCE[0]}" )"
    local SOURCE_DIR="$( dirname "$SOURCE" )"
    local REF_DIR="$( realpath "$SOURCE_DIR/../References/" )"

    local numRes="$(grep "^$QUERY" "${SOURCE_DIR}/contaminants.csv" | wc -l)"
    if [ $numRes -eq 0 ]; then
        echo "${QUERY}_not_found"
        return
    fi
    eecho "Found $numRes fitting generator targets for query"
    while read line; do
        generate_reads "$line" $runmode
    done < <(grep "^$QUERY" "${SOURCE_DIR}/contaminants.csv")

}

generate_reads(){
    #echo "$#:$*"
    IFS=';' read -a args <<< $1
    local numargs=${#args[@]}
    [ $numargs -eq 0 ] && return

    local runmode="${2-"info"}"
    local alt_runmode="info"
    [[ $runmode =~ g.* ]] && alt_runmode=$runmode
    local out=${args[0]}
    local ref=${args[1]}
    ref="$( "$REF_DIR/get_ref.sh" "$ref" "$alt_runmode" )"
    if [ ! -f $ref ]; then
        echo "${out}_error_input_reference_missing"
        return
    fi
    local tech=${args[2]}
    local opts=${args[3]}
    local seedopt=()
    local seed=${args[4]-"-1"}
    local outprefix="${SOURCE_DIR}/${out}"
    [ $seed -ne -1 ] && seedopt=( "-rs" "$seed" )
    local paired=0
    [[ "$opts" =~ "-m" ]] && paired=1

    opts=( $opts )
    case $2 in
        i*)
            [ $paired -eq 1 ] && printf "%s;%s\n" "${outprefix}"{1,2}.fq.gz \
                || echo "${outprefix}.fq.gz" 
            ;;
        r*|g*)
            if [ ! -f "${outprefix}1.fq.gz" ] && [ $paired -eq 1 ]; then
                art_illumina -ss "$tech" -i "$ref" ${opts[@]} ${seedopt[@]} -o "$outprefix" 1>&2
                gzip "${outprefix}1.fq" & gzip "${outprefix}2.fq"
            elif [ ! -f "${outprefix}.fq.gz" ] && [ $paired -eq 0 ]; then
                art_illumina -ss "$tech" -i "$ref" ${opts[@]} ${seedopt[@]} -o "$outprefix" 1>&2
                gzip "${outprefix}.fq"
            fi
            [ $paired -eq 1 ] && printf "%s;%s\n" "${outprefix}"{1,2}.fq.gz \
                || echo "${outprefix}.fq.gz" 
            ;;
        c*)
            echo art_illumina -ss "$tech" -i "$ref" ${opts[@]} ${seedopt[@]} -o "$outprefix"
            ;;



    esac
    #art_illumina -ss MSv3 -sam -i ../References/NC_045512.2.fasta -l 250 -f 50 -m 300 -s 30  -o "VIRAL_RAW_ALLVARS_50X"
}


main $@
