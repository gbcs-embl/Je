#jars found in this folder are artifact that are not found in maven central, you can then puch them in your local maven repo with the following commands:

#ADAPT fpath to YOUR Je/lib
LIBPATH="/Users/girardot/Work/eclipse_ws/Je/lib/"
cd ~/.m2 
mvn install:install-file -DgroupId=net.sf -DartifactId=htsjdk -Dversion=1.140custom -Dfile=$LIBPATH/custom-picard/htsjdk-1.140.jar -Dpackaging=jar -DgeneratePom=true 
mvn install:install-file -DgroupId=net.sf -DartifactId=picard -Dversion=1.140custom -Dfile=$LIBPATH/custom-picard/picard.jar -Dpackaging=jar -DgeneratePom=true

# Uncomment to ADD GBCS artifacts if needed (ie if you don t have access to these repos)
# IF you are at embl, you rather want to checkout the relevant projects and build them locally
# mvn install:install-file -DgroupId=org.embl.cg.utilitytools -DartifactId=ut_utils -Dversion=1.0 -Dfile=$LIBPATH/ut_utils-1.0.jar -Dpackaging=jar -DgeneratePom=true 
# mvn install:install-file -DgroupId=org.embl.gbcs.embase -DartifactId=jembase_api -Dversion=1.0 -Dfile=$LIBPATH/jembase_api-1.0.jar -Dpackaging=jar -DgeneratePom=true

