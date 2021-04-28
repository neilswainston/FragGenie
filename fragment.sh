BASEDIR=$(dirname "$0")

mvn -f $BASEDIR/pom_jar.xml install

java \
	-Xmx32000m \
	-cp $BASEDIR/target/liv-metfrag-0.1.0-SNAPSHOT.jar:$BASEDIR/target/lib/* \
	uk.ac.liverpool.metfrag.MetFragFragmenter \
	"$@"