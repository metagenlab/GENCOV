#!/bin/bash

set -eEuo pipefail


# @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>
# @SIM:1:FCX:1:15:6329:1045 1:N:0:2

eecho(){
    echo $@ 1>&2
}

main(){
    local scrname="GEN_SAMPLES SCRIPT"
    local QUERY="${1-".*"}"
    local runmode="${2-"info"}"
    local SOURCE="$( realpath "${BASH_SOURCE[0]}" )"
    local SOURCE_DIR="$( dirname "$SOURCE" )"

    local REG_FILE="${SOURCE_DIR}/mixture.csv"
    local CONTA_SCRIPT="$( realpath "$SOURCE_DIR/../CONTAMINENTS/simulate_reads.sh" )"
    local TARGET_SCRIPT="$( realpath "$SOURCE_DIR/../TARGET_RAW/simulate_reads.sh" )"

    local numRes="$( grep -i "^$QUERY" $REG_FILE | wc -l)"
    eecho "${scrname}:Found $numRes fitting generator targets for query"
    if [ $numRes -eq 0 ]; then
        eecho "${scrname}:${QUERY} not found in Mixture Folder"
        echo "${QUERY}_not_found"
        return
    fi
    while read line; do
        eecho -e "\n### Generation $(cut -d';' -f1 <<< "$line")  ###\n"
        generate_sample "$line" $runmode
    done < <(grep -i "^$QUERY" "$REG_FILE" )
}



generate_sample(){
    IFS=';' read -a args <<< $1
    local numargs=${#args[@]}
    local runmode="${2-"info"}"
    local alt_runmode="info"
    [[ $runmode =~ g.* ]] && alt_runmode=$runmode
    local id=${args[0]}
    local ratio=${args[3]}

    local src=${args[1]}
    local temp_src="$( $SOURCE "$src" "$alt_runmode" )"
    [[ "$temp_src" =~ 'not_found'  ]] \
        && temp_src="$( $TARGET_SCRIPT "$src" "$alt_runmode" )"
    [[ "$temp_src" =~ 'not_found'  ]] \
        && temp_src="$( $CONTA_SCRIPT "$src" "$alt_runmode" )"
    if [[ ! "$temp_src" =~ 'not_found'  ]]; then
        IFS=';' read -a src <<< $temp_src
    else
       eecho -e "Target Data: $src\n Was not found anywhere"
       echo "${id}_input_missing"
       return
    fi

    # Locate Contaminent Data
    local cont=${args[2]}
    local temp_cont="$( $SOURCE "$cont" "$alt_runmode" )"
    [[ "$temp_cont" =~ 'not_found'  ]] \
        && temp_cont="$( $CONTA_SCRIPT "$cont" "$alt_runmode" )"
    if [[ ! "$temp_cont" =~ 'not_found'  ]]; then
        IFS=';' read -a cont <<< $temp_cont
    else
       eecho "Contaminant Data: $cont\n Was not found anywhere"
       echo "${id}_input_missing"
       return
    fi
    eecho sample_mix_and_fix_paired "$id" ${src[@]} ${cont[@]} ${ratio}
    sample_mix_and_fix_paired "$id" ${src[@]} ${cont[@]} ${ratio}


}
sample_mix_and_fix_paired(){
    local s1="$1"
    local src1="$2"
    local src2="$3"
    local cont1="$4"
    local cont2="$5"
    local ratio="$6"
    local fastq_out1="${SOURCE_DIR}/${s1}_S1_L001_R1_001.fastq"
    local fastq_out2="${SOURCE_DIR}/${s1}_S1_L001_R2_001.fastq"
    local seed=1234
    local instrument="SIM_MSv3"
    local run_number="1"

    case $runmode in
        i*)
            printf "%s;%s\n" "${fastq_out1}.gz"\
                             "${fastq_out2}.gz"
            ;;
        r*|g* )
            if [ ! -f "${fastq_out1}.gz" ]; then
                local num_reads=( $(( $(zcat "${src1}" | wc -l  | cut -f1 -d' ') / 4 )) )
                local cont_main_rate="$ratio"
                local cont_main_reads=$( bc <<< "$num_reads / (1 - ${cont_main_rate})" )
                sample_fastq "${src1}" "${src2}" "${fastq_out1}" "${fastq_out2}" ${ratio}
                shuffle_fastq "${fastq_out1}" "${seed}"
                fix_paired_ordering "${fastq_out1}" "${fastq_out2}"
                local tot_reads=$(( $(wc -l "${fastq_out1}" | cut -f1 -d' ') / 4 ))
                illuminafy_readnames 1 "${fastq_out1}" "$tot_reads"
                illuminafy_readnames 2 "${fastq_out2}" "$tot_reads"
                gzip "${fastq_out1}" & gzip "${fastq_out2}"
            fi
            printf "%s;%s\n" "${fastq_out1}.gz"\
                             "${fastq_out2}.gz"
            ;;
        c*)
            eecho -e "    vvv  CMDLINE  vvv\n"
            printf "%s\n" \
                "local cont_main_rate=$ratio"\
                "local num_reads=( \$(( \$( zcat \"${src1}\" | wc -l  | cut -f1 -d' ') / 4 )) )"\
                "local cont_main_reads=\$( bc <<< \"\$num_reads / (1 - \${cont_main_rate})\" )"\
                "sample_fastq ${src1} ${src2} ${fastq_out1} ${fastq_out2} ${ratio}"\
                "shuffle_fastq ${fastq_out1} ${seed}"\
                "fix_paired_ordering \"${fastq_out1}\" \"${fastq_out2}\""\
                "local tot_reads=\$(( \$(wc -l "${fastq_out1}" | cut -f1 -d' ') / 4 ))"\
                "illuminafy_readnames 1 \"${fastq_out1}\" \$tot_reads"\
                "illuminafy_readnames 2 \"${fastq_out2}\" \$tot_reads"\
                "gzip \"${fastq_out1}\" & gzip \"${fastq_out2}\""

            ;;
    esac
}


sample_fastq(){

    [ -f "${fastq_out1}" ] || (
        reformat.sh \
            ow=t \
            upsample=t\
            in1="${cont1}" in2="${cont2}" \
            out1="${fastq_out1}" \
            out2="${fastq_out2}" \
            samplereadstarget="${cont_main_reads}"\
            sampleseed=$seed\
            -Xmx2G

        zcat "${src1}" >> "${fastq_out1}"
        zcat "${src2}" >> "${fastq_out2}"
        )


}


illuminafy_readnames(){
    printf "%s\n" \
       "Fixing Readnames to correspond to illumina specifcation" \
       " "\
       "Running For File:"\
       "    ${2}" 1>&2
    local read="$1"
    local fqin="$2"
    local fqout="$fqin"
    local tot_reads="$3"
    local fqtemp=${fqin%%.fastq*}_temp.fastq
    local rename_csv=${fqin%%.fastq*}_rename.tsv
    local lane=3
    local tiles="3,15,14,13,7,8"
    local flowcell="FCX"
    local run_number=1
    awk \
        -vinstrument="$instrument" \
        -vrun="$run_number" \
        -vflowcell="$flowcell" \
        -vlane="$lane"\
        -vtiles="$tiles"\
        -vtot_reads="$tot_reads"\
        -vread="$read"\
        -vrename_csv="$rename_csv"\
        -vfqout="$fqtemp"\
    'BEGIN{
        num_tiles=split(tiles, tile_arr, ",")
        tile_cap=sprintf("%.f",tot_reads/(num_tiles-1)+.5)
        limits=sprintf("%.f", sqrt(tile_cap)+.5)
        spacer=20
     }
     {
        if( (NR-1) % 4 == 0){
            readNum=((NR-1)/4)
            modNum=(readNum % tile_cap)
            newName=sprintf("@%s:%i:%s:%i:%i:%i:%i %i:N:0:2",
                            instrument, run, flowcell, lane,
                            tile_arr[int(readNum/tile_cap) +1],
                            ((modNum % limits) +1)*spacer,
                            (int(modNum / limits) +1)*spacer,
                            read)
            printf "%s\t%s\n", newName, $0 > rename_csv
            print newName > fqout
        } else {
            print $0 > fqout
        }
    }' "${fqin}"
    mv "$fqtemp" "$fqout"
}

shuffle_fastq(){
    local fastq="$1"
    local sort_seed="$2"
    fastq-sort -R --seed="${sort_seed}" "${fastq_out1}" > "${fastq_out1}_temp"
    mv "${fastq}_temp" "${fastq}"
}

fix_paired_ordering(){
    local fq1="$1"
    local fq2="$2"
    local fqtemp1=${fq1%%.fastq*}_temp.fastq
    local fqtemp2=${fq2%%.fastq*}_temp.fastq
    repair.sh \
        in="${fq1}" in2="${fq2}" \
        out="${fqtemp1}" out2="${fqtemp2}" \
        -Xmx2G
    mv "${fqtemp1}" "${fq1}"
    mv "${fqtemp2}" "${fq2}"

}


main $@

