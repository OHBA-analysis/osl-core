#!/bin/bash
# Romesh Abeysuriya 2017

retrieve_repo () {
    # If no release tag is provided, use master branch
    REPONAME=$1
    curl -L -k https://github.com/OHBA-analysis/$REPONAME/tarball/master --output $REPONAME.tar.gz
    echo $REPONAME = `tar -tf $REPONAME.tar.gz | grep -o '^[^/]\+' | sort -u | rev | cut -d'-' -f1 | rev` >> version.txt
    rm -rf "$REPONAME"/*
    mkdir $REPONAME 
    tar -xvf $REPONAME.tar.gz --strip 1 -C $REPONAME 
    rm -rf $REPONAME.tar.gz
} 

echo "This script will upgrade your OSL code to the latest version on GitHub"
echo "Other files including example_data, std_masks, and spm12 will not be modified"
echo "As this is not a packaged release, it is possible (although unlikely) that this will result in internal incompatibilities"
echo "If you encounter any unexpected behaviour, you may need to obtain the latest *complete* OSL release"
echo
read -r -p "Perform in-place OSL upgrade? [y/N] " response
case "$response" in
    [yY][eE][sS]|[yY]) 
        wd=${PWD##*/}
        if [ "$wd" == "osl-core" ]; then

            if [ -d ".git" ]; then
                echo "ERROR - A .git folder is present. Updating should be done via git"
                exit 1
            fi

            cd .. # This is now osldir
            echo `head -n 1 version.txt | cut -f1 -d'+'`\+UPGRADE > version.txt
            retrieve_repo osl-core 
            retrieve_repo ohba-external 
            retrieve_repo GLEAN 
            retrieve_repo MEG-ROI-nets 
            retrieve_repo HMM-MAR 
            echo
            echo "Upgrade complete"
        else
            echo "ERROR - Update script must be run from within the osl-core folder"
        fi
        ;;
    *)
        echo "Cancelled"
        ;;
esac


