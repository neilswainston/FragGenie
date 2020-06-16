package uk.ac.liverpool.metfraglib.candidate;

import org.openscience.cdk.interfaces.IAtomContainer;

import de.ipbhalle.metfraglib.additionals.MoleculeFunctions;
import uk.ac.liverpool.metfraglib.precursor.Precursor;

public class PrecursorCandidate {

	private final Precursor structure;

	public PrecursorCandidate(final String smiles) throws Exception {
		IAtomContainer molecule = null;

		int trials = 1;

		while (trials <= 10) {
			try {
				molecule = MoleculeFunctions.getAtomContainerFromSMILES(smiles);
			} catch (Exception e) {
				trials++;
				e.printStackTrace();
				continue;
			}
			break;
		}

		if (molecule == null) {
			throw new Exception("Could not read SMILES!");
		}

		MoleculeFunctions.prepareAtomContainer(molecule, true);
		MoleculeFunctions.convertExplicitToImplicitHydrogens(molecule);
		
		this.structure = new Precursor(molecule);
		this.structure.preprocessPrecursor();
	}

	/**
	 * 
	 * @return
	 */
	public Precursor getPrecursorMolecule() {
		return this.structure;
	}
}