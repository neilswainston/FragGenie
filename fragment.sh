java \
	-cp target/liv-metfrag-0.1.0-SNAPSHOT/WEB-INF/classes/:target/liv-metfrag-0.1.0-SNAPSHOT/WEB-INF/lib/* \
	uk.ac.liverpool.metfrag.MetFragFragmenter \
	$1 \
	out/spectra.csv \
	$2 \
	$3