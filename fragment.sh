mvn install

java \
	-Xmx32000m \
	-cp target/liv-metfrag-0.1.0-SNAPSHOT.war:target/lib/* \
	uk.ac.liverpool.metfrag.MetFragFragmenter \
	"$@"