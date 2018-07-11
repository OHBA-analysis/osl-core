#!/bin/bash

# use HTTPS by default to allow non-members to clone the repos
if [ $# -gt 0 ]; then 
	GITHUB_URL=git@github.com:OHBA-analysis
else 
	GITHUB_URL=https://github.com/OHBA-analysis
fi

replace_repo () {
	if [ -d "$1" ]; then
		rm -rf "$1"
		git clone "${GITHUB_URL}/$1"
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
		cd osl-core
        ;;
    *)
        echo "Cancelled"
        ;;
esac
