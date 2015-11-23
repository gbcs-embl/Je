#!/bin/bash
#
# Makes sure that the custom and picard/htsjdk deps are push to the local mvn repo 
#
MVN=/usr/local/apache-maven/apache-maven-3.3.3/bin/mvn
MVN_ROOT_REPOS=~/.m2
# this need to be executed AFTER you ant build picard project sucessfully
PICARD_GIT_PATH="/Users/girardot/git/picard"
PICARD_JAR=$PICARD_GIT_PATH"/dist/picard.jar"
HTSJDK_JAR=$PICARD_GIT_PATH"/dist/htsjdk_lib_dir/htsjdk-1.140.jar"
JE_LIB_PATH="/Users/girardot/git/Je/lib/custom-picard"

cd $MVN_ROOT_REPOS
if [ -e "$PICARD_JAR" ]
then
	echo "copying $PICARD_JAR to $JE_LIB_PATH"
	cp $PICARD_JAR $JE_LIB_PATH
	echo "installing $PICARD_JAR in maven local repos"
	$MVN install:install-file -DgroupId=net.sf -DartifactId=picard -Dversion=1.140custom -Dfile=$PICARD_JAR -Dpackaging=jar -DgeneratePom=true 
fi

if [ -e "$HTSJDK_JAR" ]
then
	echo "copying $HTSJDK_JAR to $JE_LIB_PATH"
	cp $HTSJDK_JAR $JE_LIB_PATH
	echo "installing $HTSJDK_JAR in maven local repos"
	$MVN install:install-file -DgroupId=net.sf -DartifactId=htsjdk -Dversion=1.140custom -Dfile=$HTSJDK_JAR -Dpackaging=jar -DgeneratePom=true 
fi
