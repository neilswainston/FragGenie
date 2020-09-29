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
	 * @param smiles
	 * @return IAtomContainer
	 * @throws CDKException 
	 * @throws Exception
	 */
	public static IAtomContainer getAtomContainer(final String smiles) throws CDKException {
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