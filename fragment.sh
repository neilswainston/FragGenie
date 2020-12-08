mvn -f pom_jar.xml install

java \
	-Xmx32000m \
	-cp target/liv-metfrag-0.1.0-SNAPSHOT.jar:target/lib/* \
	uk.ac.liverpool.metfrag.MetFragFragmenter \
	"$@"