######
# settings

WORKINGDIR=/Users/woolrich/Desktop

# specify these 5 things to build a release:
OSL2GITRELEASETAG=0.3  # e.g corresponds to release/tag (e.g. v0.3) created on github website
GLEANGITRELEASETAG=0.2 # e.g. corresponds to release/tag (e.g. v0.2) created on github website, note that GLEAN contains HMMMAR as a subproject
MEGROInetsGITRELEASETAG=1.6.0 # e.g. corresponds to release/tag (e.g. v1.6.0) created on github website 
SUPPDIRSNAME=osl2_supporting_dirs_v2 # name of supplementary directories to use, do not want .tar.gz file ext here
RELEASENAME=osl2.${OSL2GITRELEASETAG}.1 # overall name of built package

# location of repository for supplementary directories
SUPPDIRSDOWNLOADSITE=http://users.fmrib.ox.ac.uk/~woolrich/osl2

# where we will upload osl release to
OSLUPLOADDIR=woolrich@jalapeno.fmrib.ox.ac.uk:/home/fs0/woolrich/www/osl2

GITDOWNLOADSITE=https://github.com/OHBA-analysis/

######
# make dir to work in

mkdir $WORKINGDIR

######
# download current tar ball of supplementary dirs for osl2

cd $WORKINGDIR

rm -rf supporting_dirs.tar.gz
curl -L $SUPPDIRSDOWNLOADSITE/$SUPPDIRSNAME.tar.gz --output supporting_dirs.tar.gz

# untar it
rm -rf supporting_dirs

tar xvf supporting_dirs.tar.gz
rm -rf supporting_dirs.tar.gz

echo "function ret=osl2_supporting_dirs_version, ret=‘$SUPPDIRSDOWNLOADSITE/$SUPPDIRSNAME';" > osl2/osl2_supporting_dirs_version.m

######
# download current github osl release

REPONAME=osl2
GITRELEASETAG=$OSL2GITRELEASETAG

cd $WORKINGDIR

rm -rf osl2/$REPONAME.tar.gz
rm -rf osl2/$REPONAME

tmp=$GITDOWNLOADSITE$REPONAME/archive/v$GITRELEASETAG.tar.gz

echo Downloading $tmp

curl -L $tmp --output osl2/$REPONAME.tar.gz

cd osl2

# untar it
rm -rf $REPONAME-$GITRELEASETAG
tar xvf $REPONAME.tar.gz
rm -rf $REPONAME.tar.gz

# rename the internal osl dir osl2
mv $REPONAME-$GITRELEASETAG $REPONAME

# update $REPONAME_version.m
cd $WORKINGDIR
rm -f osl2/$REPONAME/${REPONAME}_version.m
echo "function ret=${REPONAME}_version, ret='$REPONAME.$GITRELEASETAG';" > osl2/$REPONAME/${REPONAME}_version.m

# also put a $REPONAME_version.m in the root osl2 dir
rm -f osl2/${REPONAME}_version.m
echo "function ret=${REPONAME}_version, ret='$tmp';" > osl2/$REPONAME/${REPONAME}_version.m

######
# download current GLEAN osl release
# note that GLEAN contains HMMMAR as a subproject

REPONAME=GLEAN
GITRELEASETAG=$GLEANGITRELEASETAG

cd $WORKINGDIR

rm -rf osl2/$REPONAME.tar.gz
rm -rf osl2/$REPONAME

tmp=$GITDOWNLOADSITE$REPONAME/archive/v$GITRELEASETAG.tar.gz

echo Downloading $tmp

curl -L $tmp --output osl2/$REPONAME.tar.gz

cd osl2

# untar it
rm -rf $REPONAME-$GITRELEASETAG
tar xvf $REPONAME.tar.gz
rm -rf $REPONAME.tar.gz

# rename the internal osl dir osl2
mv $REPONAME-$GITRELEASETAG $REPONAME

# update $REPONAME_version.m
cd $WORKINGDIR
rm -f osl2/$REPONAME/${REPONAME}_version.m
echo "function ret=${REPONAME}_version, ret='$tmp';" > osl2/$REPONAME/${REPONAME}_version.m

######
# download current MEG-ROI-nets osl release

REPONAME=MEG-ROI-nets
GITRELEASETAG=$MEGROInetsGITRELEASETAG

cd $WORKINGDIR

rm -rf osl2/$REPONAME.tar.gz
rm -rf osl2/$REPONAME

tmp=$GITDOWNLOADSITE$REPONAME/archive/v$GITRELEASETAG.tar.gz

echo Downloading $tmp

curl -L $tmp --output osl2/$REPONAME.tar.gz

cd osl2

# untar it
rm -rf $REPONAME-$GITRELEASETAG
tar xvf $REPONAME.tar.gz
rm -rf $REPONAME.tar.gz

# rename the internal osl dir osl2
mv $REPONAME-$GITRELEASETAG $REPONAME

# update $REPONAME_version.m
cd $WORKINGDIR
rm -f osl2/$REPONAME/${REPONAME}_version.m
echo "function ret=${REPONAME}_version, ret='$tmp';" > osl2/$REPONAME/${REPONAME}_version.m

######
# create tar ball of osl

echo "function ret=osl2_release_version, ret='$RELEASENAME';" > osl2/osl2_release_version.m
cd $WORKINGDIR
tar cvf $RELEASENAME.tar osl2
gzip -f $RELEASENAME.tar

######
# upload - this bit obviously requires write permission for $OSLUPLOADDIR

scp -r $RELEASENAME.tar.gz $OSLUPLOADDIR/

echo Built release can be downloaded from $OSLUPLOADDIR/$RELEASENAME.tar.gz

