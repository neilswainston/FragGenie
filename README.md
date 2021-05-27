# liv-metfrag

liv-metfrag provides a simple codebase and application for generating a list of
theoretical chemical fragments which can be used to generate theoretical mass
spectra of chemical compounds.

liv-metfrag can be accessed directly through the underlying Java code or run
from the command line.

An example script for running the application is provided (`test.sh`) along with
a simple example input file `test_input.csv`.

The `test.sh` example script runs the more general `fragment.sh` script with the
following parameters:

`./fragment.sh test_input.csv test_output.csv smiles 3 80.0 256 1000000000 true METFRAG_MZ`

where:
* `test_input.csv` is the input file name.
* `test_output.csv` is the output file name.
* `smiles` is the header name of the column in the input csv file containing
SMILES string representations of the chemicals of interest.
* `3` is the fragmentation recursion depth.
* `80.0` is the minimum fragment mass to return.
* `256` is the maximum length SMILES string to consider.
* `1000000000` is the maximum number of SMILES strings to analyse.
* `true` is a flag to indicate whether to run the script in silent mode.
* `METFRAG_MZ` specifies that masses (mz) values should be generated. Other
options are `METFRAG_FORMULAE` which generates fragment formulae and
`METFRAG_BROKEN_BONDS` which generates a list of bonds broken for each fragment.
Any combination of these three parameters may be specified.