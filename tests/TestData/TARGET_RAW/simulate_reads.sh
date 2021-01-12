#!/bin/bash

set -eEuo pipefail


eecho(){
    echo $@ 1>&2
}

main(){
    local scrname="TARGET SCRIPT"
    local OLD_IFS=$IFS
    local QUERY="${1-".*"}"
    local runmode="${2-"info"}"
    local SOURCE="$( realpath "${BASH_SOURCE[0]}" )"
    local SOURCE_DIR="$( dirname "$SOURCE" )"
    local REF_DIR="$( realpath "$SOURCE_DIR/../References/" )"
    local MOD_REF="$( realpath "$SOURCE_DIR/../MODIFIED_REFERENCE/" )"
    local numRes="$(grep "^$QUERY" "${SOURCE_DIR}/targets.csv" | wc -l)"
    if [ $numRes -eq 0 ]; then
        eecho "${scrname}:${QUERY} not found in TARGET FOLDER"
        echo "${QUERY}_not_found"
        return
    fi
    eecho "${scrname}:Found $numRes fitting generator targets for query"
    eecho "${scrname}: Runmode $runmode"
    while read line; do
        generate_reads "$line" "$runmode"
    done < <(grep "^$QUERY" "${SOURCE_DIR}/targets.csv")
}



generate_reads(){
    #echo "$#:$*"
    IFS=';' read -a args <<< $1
    local numargs=${#args[@]}
    [ $numargs -eq 0 ] && return

    local out=${args[0]}
    local ref=${args[1]}
    local temp_ref="$( "$MOD_REF/get_modref.sh" "$ref" "$runmode" )"
    [ ! -f $temp_ref ] && \
        temp_ref="$( "$REF_DIR/get_ref.sh" "$ref" "$runmode" )"
    if [ ! -f $temp_ref ]; then
        echo "${out}_error_input_reference_missing"
        return
    fi
    ref="$temp_ref"
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
                gzip "${SOURCE_DIR}/${out}.fq"
            fi
            [ $paired -eq 1 ] && printf "%s;%s\n" "${outprefix}"{1,2}.fq.gz \
                || echo "${outprefix}.fq.gz" 
            ;;
        c*)
            echo art_illumina -ss "$tech" -i "$ref" ${opts[@]} ${seedopt[@]} -o "$outprefix"
            ;;



    esac

}



#[ -f "VIRAL_RAW_ALLVARS_50X_1.fq"   ] || art_illumina -ss MSv3 -sam -i ../../References/NC_045512.2.fasta -l 250 -f 50 -m 300 -s 30  -o "VIRAL_RAW_ALLVARS_50X_"
#[ -f "VIRAL_RAW_ALLVARS_50X_150l_MSv3_1.fq"   ] || art_illumina -ss MSv3 -sam -i ../../References/NC_045512.2.fasta -l 150 -f 50 -m 200 -s 30  -o "VIRAL_RAW_ALLVARS_50X_r150_MSv3_"
#art_illumina -ss MSv3 -sam -i ../References/hg19.fa 250 -f 50 -m 600 -s 20 -o "HUMAN_RAW_50X"
main $@

