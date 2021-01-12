#!/bin/bash

set -eEuo pipefail
RET_VAL=0
SOURCE="$( realpath "${BASH_SOURCE[0]}" )"
SOURCE_DIR="$( dirname "$SOURCE" )"
STATE=0
STATE_FILE=".release_in_progress"

eecho(){
    echo $@ 1>&2
}

stop_process(){
    local rv=${RET_VAL-0}
    if [ $STATE -ne 0 ]; then
       printf "%s\n" \
          "Something went wrong during the release process." \
          "You can pick up where you left by rerunning the script." \
          " "\
          "Note that if you made it past step 1 of the progress," \
          "the version file has already been updated to the newest release." \
          "Do not be tricked by your software already reporting the new" \
          "version back to you."
          " " 1>&2
    fi
    unset RET_VAL STATE SOURCE SOURCE_DIR STATE_FILE
    exit $rv
}

had_error(){
    set +x
    ERROR_TYPE=1
    local errortype=""
    local prov_msg="${2:-""}"
    local msg=""
    case $1 in
        1)
            eecho "$prov_msg"
            errortype="RUNTIME ERROR:"
            local call=( $(caller) )
            # print surrounding lines of code
            # Stack Overflow
            # https://unix.stackexchange.com/questions/39623/trap-err-and-echoing-the-error-line
            # Thanks to unpythonic
            # https://unix.stackexchange.com/users/8365/unpythonic
            msg=$(
                awk 'NR>L-4 && NR<L+4 {
                        printf "%-5d%3s%s\n",NR,(NR==L?">>>":""),$SOURCE
                     }' L="${call[0]}" "$SOURCE" )
            prov_msg="Line "${call[@]}"\n"
        ;;
        2)
            errortype="USER INTERUPT:"
        ;;
        3)
            errortype="STRANGE ERROR:"
        ;;
    esac
    eecho -e "$errortype" "$prov_msg" "$msg"
    stop_process
}

trap 'had_error 2 "User requested Interupt"' SIGINT
trap 'had_error 1' ERR

main(){
    local release_branch="master"
    local patch=0
    local minor=0
    local major=0
    local prior_version=( $(read_in_version_file) )
    get_state
    major=${prior_version[0]}
    minor=${prior_version[1]}
    patch=${prior_version[2]}
    local curr_branch="$(git rev-parse --abbrev-ref HEAD)"
    local last_tag="$(get_latest_tag)"
    local release_mode="patch"
    [ $STATE -eq 0 ] && increment_version ${1-$release_mode}
    eecho "Current Version File: v${prior_version[0]}.${prior_version[1]}.${prior_version[2]}"
    eecho "Last Tag: ${last_tag}"
    eecho "Trying to create new ${release_mode} release! (v${major}.${minor}.${patch})"
    eecho "   You are in branch: ${curr_branch}"
    if [[ ! "$curr_branch" =~ "development_unstable" ]] \
            && [[ ! "$curr_branch" =~ "master" ]]; then
        eecho "Error: Can only issue new releases from master or development_unstable branch"
    elif [[ $STATE -eq 0 ]]; then

        printf "%s\n" \
            "You are about to release Version ${major}.${minor}.${patch}"\
            "This is a ${release_mode} release."\
            " "\
            "Ensure that you have run all required tests"\
            "and cooridinated with the other members of the project"\
            " " 1>&2
        #local input="N"
        read -r -p"Do you still want to continue? [Y/n]" input
        case "$input" in
            [yY][eE][sS]|[yY])
                eecho "Triggered release"
                release
                ;;
            [nN][oO]|[nN])
                eecho "Aborted release"
                exit 0
                ;;
        esac
    else
        printf "%s\n" \
            "Continuing unfinisched Release: ${major}.${minor}.${patch}" \
            "You are currently in stage ${STATE} of 11" \
            " " 1>&2
        release
    fi
}

increment_version(){
    case $1 in
        ma*|major)
            release_mode="major"
            major=$(( major+1  ))
            minor=0
            patch=0
            ;;
        mi*|minor)
            release_mode="minor"
            minor=$(( minor+1 ))
            patch=0
            ;;
        p*|patch)
            release_mode="patch"
            patch=$(( patch+1 ))
            ;;
        *)
            release_mode="patch"
            patch=$(( patch+1 ))
            ;;
    esac
}

release(){
    checkpoint  1 write_version_file
	checkpoint  2 setup_old_conda 
	checkpoint  3 setup_new_conda 
    checkpoint  4 git add ${SOURCE_DIR}/__version__.py
    checkpoint  5 git add ${SOURCE_DIR}/.conda/*
    checkpoint  6 git commit --amend --no-edit
    checkpoint  7 git tag -a "v${major}.${minor}.${patch}"
    checkpoint  8 git pull origin ${release_branch}
    checkpoint  9 git checkout ${release_branch}
    checkpoint 10 git merge ${curr_branch} --ff-only
    checkpoint 11 git push
    checkpoint 12 git push origin "v${major}.${minor}.${patch}"
    checkpoint 13 git checkout ${curr_branch}
	checkpoint 14 git pull
    checkpoint 15 git push
    release_done
}

get_state(){
   if [ ! -f "${SOURCE_DIR}/${STATE_FILE}" ]; then
        STATE=0
   else
        local pre_state="$(grep -oP 'ReleaseState=\K\d+' "${SOURCE_DIR}/${STATE_FILE}")"
        if [ -z $pre_state ]; then
            eecho "Warning: Release Already in Progress but state file is corrupt"
            STATE=1
        else
            STATE="$pre_state"
        fi

   fi

}

checkpoint(){
    point=$1
    if [ $point -gt $STATE ]; then
        ${@:2}
        STATE=$point
        echo "ReleaseState=${STATE}" > "${SOURCE_DIR}/${STATE_FILE}"
    else
        eecho "SKIPPING: ${@:2}"
    fi
}

release_done(){
    eecho "Release Completed Successfully!!"
    rm "${SOURCE_DIR}/${STATE_FILE}"
    STATE=0
}

read_in_version_file(){
    awk -vFS='=' '/^MAJOR/{major=$2}
                  /^MINOR/{minor=$2}
                  /^PATCH/{patch=$2}
                  END{printf "%s %s %s", major, minor, patch}' \
        "${SOURCE_DIR}/__version__.py"
}

write_version_file(){
    sed -i 's/\(^MAJOR=\).*/\1'"${major}/g"'
            s/\(^MINOR=\).*/\1'"${minor}/g"'
            s/\(^PATCH=\).*/\1'"${patch}/g" "${SOURCE_DIR}/__version__.py"
}

setup_old_conda(){
	local OLD_RECIPE_DIR="${SOURCE_DIR}/.conda/${prior_version[0]}.${prior_version[1]}.${prior_version[2]}" 
	mkdir -p "${OLD_RECIPE_DIR}"
	cp -f .conda/meta.yaml .conda/build.sh "${OLD_RECIPE_DIR}/."
}

setup_new_conda(){
	first_line="{% set version = \"${major}.${minor}.${patch}\" %}"
    sed -i '1 s/^.*$/'"${first_line}/" .conda/meta.yaml	
}

get_latest_tag(){
    git describe --long  | cut -d'-' -f1
}


main $@
