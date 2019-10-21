#!/bin/bash

DEBUG=0
GITHUB_URL=https://github.com/OHBA-analysis

replace_repo () {
    if [ -d "$1" ]; then
        if (( $DEBUG )); then 
            echo "Cloning repository '${GITHUB_URL}/$1' into folder '$1'"
        else 
            rm -rf "$1"
            git clone "${GITHUB_URL}/$1"
        fi
    fi
} 

read -r -p 'This will delete all OSL folders and clone from git. Use with caution! Continue? [y/N] ' response
case "$response" in
    [yY][eE][sS]|[yY]) 

        read -r -p 'Is this computer associated with a GitHub account, which has access to the OHBA Analysis Group repositories? [y/N] ' response
        case "$response" in
            [yY][eE][sS]|[yY])
                GITHUB_URL=git@github.com:OHBA-analysis
                ;;
        esac

        cd ..
        replace_repo osl-core
        replace_repo ohba-external
        replace_repo GLEAN
        replace_repo MEG-ROI-nets
        replace_repo HMM-MAR
        cd osl-core
        ;;
    *)
        echo 'Cancelled'
        ;;
esac
