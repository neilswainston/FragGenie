package uk.ac.liverpool.metfrag;

import java.util.ArrayList;
import java.util.Collection;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 * 
 * @author neilswainston
 */
public class Precursor {


	/**
	 * 
	 */
	private final IAtomContainer atomContainer;

	/**
	 * 
	 * @param smiles
	 * @throws CDKException
	 */
	Precursor(final String smiles) throws CDKException {
		this.atomContainer = getAtomContainer(smiles);
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((this.atomContainer == null) ? 0 : this.atomContainer.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
			
		if (obj == null) {
			return false;
		}
			
		if (getClass() != obj.getClass()) {
			return false;
		}
			
		final Precursor other = (Precursor) obj;
		
		if (this.atomContainer == null) {
			if (other.atomContainer != null) {
				return false;
			}
		}
		else if (!this.atomContainer.equals(other.atomContainer)) {
			return false;
		}
			
		return true;
	}

	/**
	 * 
	 * @return IAtomContainer
	 */
	IAtomContainer getAtomContainer() {
		return this.atomContainer;
	}
	
	/**
	 * 
	 * @param smiles
	 * @return IAtomContainer
	 * @throws CDKException 
	 * @throws Exception
	 */
	private static IAtomContainer getAtomContainer(final String smiles) throws CDKException {
		final SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
		final IAtomContainer molecule = parser.parseSmiles(smiles);
		final Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());

		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		aromaticity.apply(molecule);

		final CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(molecule.getBuilder());

		for (int i = 0; i < molecule.getAtomCount(); i++) {
			hAdder.addImplicitHydrogens(molecule, molecule.getAtom(i));
		}

		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);

		removeHydrogens(molecule);
		return molecule;
	}

	/**
	 * 
	 * @param molecule
	 */
	private static void removeHydrogens(final IAtomContainer molecule) {
		final Collection<IAtom> hydrogenAtoms = new ArrayList<>();

		for (IAtom atom : molecule.atoms()) {
			if (atom.getSymbol().equals("H")) { //$NON-NLS-1$
				hydrogenAtoms.add(atom);
			}

			int numberHydrogens = 0;

			for (IAtom neighbour : molecule.getConnectedAtomsList(atom)) {
				if (neighbour.getSymbol().equals("H")) { //$NON-NLS-1$
					numberHydrogens++;
				}
			}

			atom.setImplicitHydrogenCount(Integer.valueOf(numberHydrogens));
		}

		for (IAtom atom : hydrogenAtoms) {
			molecule.removeAtom(atom);
		}
	}
}