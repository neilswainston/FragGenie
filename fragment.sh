java \
	-Xmx32000m \
	-cp target/liv-metfrag-0.1.0-SNAPSHOT/WEB-INF/classes/:target/liv-metfrag-0.1.0-SNAPSHOT/WEB-INF/lib/* \
	uk.ac.liverpool.metfrag.MetFragFragmenter \
	"$@"