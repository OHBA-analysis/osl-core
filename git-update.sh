#!/bin/bash

HERE=$(pwd)
THIS=${BASH_SOURCE}
OSLCORE=$(dirname "$THIS")

repo_update() {
    echo "+ -------------------------------"
    echo "+ Folder $1..."
    [ ! -d "${OSLCORE}/../${1}/.git" ] && { echo "ERROR - Not a repository."; return 1; }
    cd "${OSLCORE}/../${1}"
    git branch -vv
    git pull
    echo "+ -------------------------------"
}

repo_update osl-core
repo_update ohba-external
repo_update GLEAN
repo_update MEG-ROI-nets 
repo_update HMM-MAR
cd "$HERE"
