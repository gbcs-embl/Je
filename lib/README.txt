#jars found in this folder are artifact that are not found in maven central, you can then puch them in your local maven repo with the following commands:

#ADAPT fpath to Je/lib
LIBPATH="/Users/girardot/Work/eclipse_ws/Je/lib/custom-picard"
cd ~/.m2 
mvn install:install-file -DgroupId=net.sf -DartifactId=htsjdk -Dversion=1.140custom -Dfile=$LIBPATH/htsjdk-1.140.jar -Dpackaging=jar -DgeneratePom=true 
mvn install:install-file -DgroupId=net.sf -DartifactId=picard -Dversion=1.140custom -Dfile=$LIBPATH/picard.jar -Dpackaging=jar -DgeneratePom=true

