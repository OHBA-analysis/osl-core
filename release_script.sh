#!/bin/bash
# Romesh Abeysuriya 2017
# MW pre-2014

# Assemble the final zip file in this directory - relative to current folder
WORKINGDIR=$HOME/Desktop/osl-release

# Specify OSL version which also names this release
OSL2RELEASETAG=0.3 
RELEASENAME=osl2.${OSL2RELEASETAG}.1 

retrieve_repo () {
	REPONAME=$1
	GITRELEASETAG=$2
	curl -L https://github.com/OHBA-analysis/$REPONAME/archive/v$GITRELEASETAG.tar.gz --output osl/$REPONAME.tar.gz
	tar -xvf osl/$REPONAME.tar.gz
	rm -rf osl/$REPONAME.tar.gz
	mv $REPONAME-$GITRELEASETAG osl/$REPONAME
	echo "$REPONAME - $GITRELEASETAG" >> osl/version.txt
} 

# Set up
mkdir -p $WORKINGDIR/osl
cd $WORKINGDIR

# Copy supplementary directories - this will also clean any previous downloads
rsync -azP --delete hbapc33:/data/analysis_data/osl/current_supplements/ osl/

# Alternatively, may want to do it via a tar file and a version - something to consider
# SUPPDIRVERSION=0.2
# rsync -a hbapc33:/data/analysis_data/osl/supplements/supplement_v$SUPPDIRVERSION.tar.gz $WORKINDIR
# tar xvf supplement_v$SUPPDIRVERSION.tar.gz
# rm -rf supplement_v$SUPPDIRVERSION.tar.gz

# Clean supplements
rm osl/update.sh
find $WORKINGDIR -name .git | xargs rm -fr
find $WORKINGDIR -name .svn | xargs rm -fr
find $WORKINGDIR -name *DS_Store* | xargs rm -fr

# Create version file
echo "OSL ZIP FILE DISTRIBUTION $RELEASENAME" > $WORKINGDIR/osl/version.txt

# Download git repos with specified release names
retrieve_repo osl2 $OSL2RELEASETAG
retrieve_repo GLEAN 0.2
retrieve_repo MEG-ROI-nets 1.6.0
retrieve_repo HMM-MAR 0.9

if which pigz >/dev/null; then
    tar cvf - ./* | pigz -9 -p 8 > $WORKINGDIR/$RELEASENAME.tar.gz
else
	tar cvf - ./* | gzip -9 > $WORKINGDIR/$RELEASENAME.tar.gz
fi

echo "Build complete"

read -r -p "Upload to Woolrich Jalapeno [y/N] " response
case "$response" in
    [yY][eE][sS]|[yY]) 
		OSLUPLOADDIR=woolrich@jalapeno.fmrib.ox.ac.uk:/home/fs0/woolrich/www/osl2
		scp -r $WORKINGDIR/$RELEASENAME.tar.gz $OSLUPLOADDIR/
        echo Built release can be downloaded from $OSLUPLOADDIR/$RELEASENAME.tar.gz
        ;;
    *)
        echo "Cancelled"
        ;;
esac


