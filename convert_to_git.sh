#!/bin/bash

replace_repo () {
	# If no release tag is provided, use master branch
	if [ -d "$1" ]; then
		rm -rf "$1"
		git clone git@github.com:OHBA-analysis/$1
	fi
} 

read -r -p "This will delete all OSL folders and clone from git. Use with caution! Continue? [y/N] " response
case "$response" in
    [yY][eE][sS]|[yY]) 
		cd ..
		replace_repo osl-core
		replace_repo ohba-external
		replace_repo GLEAN
		replace_repo MEG-ROI-nets
		replace_repo HMM-MAR
        ;;
    *)
        echo "Cancelled"
        ;;
esac


