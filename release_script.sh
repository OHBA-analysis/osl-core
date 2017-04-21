#!/bin/bash
# Romesh Abeysuriya 2017
# MW pre-2014

# Assemble the final zip file in this directory - relative to current folder
WORKINGDIR=$HOME/Desktop/osl-release

# Specify OSL version which also names this release
OSLDISTVERSION=2.2.0 

retrieve_repo () {
	# If no release tag is provided, use master branch
	REPONAME=$1
	GITRELEASETAG=$2

	if [ -z "$GITRELEASETAG" ]
	then
		curl -L https://github.com/OHBA-analysis/$REPONAME/tarball/master --output osl/$REPONAME.tar.gz
		echo $REPONAME = `tar -tf osl/$REPONAME.tar.gz | grep -o '^[^/]\+' | sort -u | rev | cut -d'-' -f1 | rev` >> osl/version.txt
	else
		curl -L https://github.com/OHBA-analysis/$REPONAME/archive/v$GITRELEASETAG.tar.gz --output osl/$REPONAME.tar.gz
		echo "$REPONAME = $GITRELEASETAG" >> osl/version.txt
	fi

	rm -rf osl/$REPONAME
	mkdir osl/$REPONAME 
	tar -xvf osl/$REPONAME.tar.gz --strip 1 -C osl/$REPONAME 
	rm -rf osl/$REPONAME.tar.gz
} 

# Set up
mkdir -p $WORKINGDIR/osl
cd $WORKINGDIR

# Copy supplementary directories - this will also clean any previous downloads
read -r -p "Include example data? [y/N] " response
case "$response" in
    [yY][eE][sS]|[yY]) 
		rsync -azP --delete hbapc33:/data/analysis_data/osl/current_supplements/ osl/
        ;;
    *)
		rsync -azP --exclude="example_data" --delete hbapc33:/data/analysis_data/osl/current_supplements/ osl/
        ;;
esac

# Clean supplements
rm osl/update.sh
find $WORKINGDIR -name .git | xargs rm -fr
find $WORKINGDIR -name .svn | xargs rm -fr
find $WORKINGDIR -name *DS_Store* | xargs rm -fr

# Create version file
echo "osl = "$OSLDISTVERSION > $WORKINGDIR/osl/version.txt

# Download git repos with specified release names
# i.e. uncomment the version numbers below to use the specified tagged releases
retrieve_repo osl-core #0.3
retrieve_repo ohba-external #0.3
retrieve_repo GLEAN #0.2
retrieve_repo MEG-ROI-nets #1.6.0
retrieve_repo HMM-MAR #0.9

rm -f $WORKINGDIR/osl_$OSLDISTVERSION.tar.gz

# pigz does parallel gzip - can be installed with 'brew install pigz'
if which pigz >/dev/null; then
    tar --exclude=osl/example_data -cvf  - osl | pigz -9 -p 8 > ${WORKINGDIR}/osl_${OSLDISTVERSION}.tar.gz
else
	tar --exclude=osl/example_data -cvf  - osl | gzip -9 > ${WORKINGDIR}/osl_${OSLDISTVERSION}.tar.gz
fi

echo "Main build complete"

# Compressing the example data is time consuming and probably doesn't need to be done very often
read -r -p "Compress example data? [y/N] " response
case "$response" in
    [yY][eE][sS]|[yY]) 
		if which pigz >/dev/null; then
		    tar -C osl -cvf  - example_data | pigz -9 -p 8 > ${WORKINGDIR}/osl_${OSLDISTVERSION}_example_data.tar.gz
		else
			tar -C osl -cvf  - example_data | gzip -9 > ${WORKINGDIR}/osl_${OSLDISTVERSION}_example_data.tar.gz
		fi
        ;;
    *)
        echo "Cancelled"
        ;;
esac

read -r -p "Upload to Woolrich Jalapeno [y/N] " response
case "$response" in
    [yY][eE][sS]|[yY]) 
		OSLUPLOADDIR=woolrich@jalapeno.fmrib.ox.ac.uk:/home/fs0/woolrich/www/osl2
		scp -r $WORKINGDIR/osl_$OSLDISTVERSION.tar.gz $OSLUPLOADDIR/
        echo Built release can be downloaded from $OSLUPLOADDIR/osl_$OSLDISTVERSION.tar.gz
        ;;
    *)
        echo "Cancelled"
        ;;
esac


