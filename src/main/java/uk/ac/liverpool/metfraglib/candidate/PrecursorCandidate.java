package uk.ac.liverpool.metfraglib.candidate;

import org.openscience.cdk.interfaces.IAtomContainer;

import de.ipbhalle.metfraglib.additionals.MoleculeFunctions;
import de.ipbhalle.metfraglib.precursor.TopDownBitArrayPrecursor;

public class PrecursorCandidate {

	private final String smiles;
	private final TopDownBitArrayPrecursor structure;

	public PrecursorCandidate(final String smiles) throws Exception {
		this.smiles = smiles;
		this.structure = new TopDownBitArrayPrecursor(this.getAtomContainer());
		this.structure.preprocessPrecursor();
	}

	public TopDownBitArrayPrecursor getPrecursorMolecule() {
		return this.structure;
	}

	private IAtomContainer getAtomContainer() throws Exception {
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
		return molecule;
	}

	@Override
	public PrecursorCandidate clone() {
		try {
			return new PrecursorCandidate(this.smiles);
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
}
