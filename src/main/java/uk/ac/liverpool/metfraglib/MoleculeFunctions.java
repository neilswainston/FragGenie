package uk.ac.liverpool.metfraglib;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

class MoleculeFunctions {

	static IAtomContainer getAtomContainer(final String smiles) throws Exception {
		final SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
		final IAtomContainer molecule = parser.parseSmiles(smiles);
		final Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		aromaticity.apply(molecule);
		
		final CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
        
        for(int i = 0; i < molecule.getAtomCount(); i++) {
        	hAdder.addImplicitHydrogens(molecule, molecule.getAtom(i));
        }
        
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
        
		removeHydrogens(molecule);
		return molecule;
	}
	
	@SuppressWarnings("boxing")
	private static void removeHydrogens(IAtomContainer molecule) {
		java.util.ArrayList<IAtom> hydrogenAtoms = new java.util.ArrayList<>();
		java.util.Iterator<IAtom> atoms = molecule.atoms().iterator();
		while(atoms.hasNext()) {
			IAtom currentAtom = atoms.next();
			if(currentAtom.getSymbol().equals("H")) hydrogenAtoms.add(currentAtom); //$NON-NLS-1$
			java.util.List<IAtom> neighbours = molecule.getConnectedAtomsList(currentAtom);
			int numberHydrogens = 0;
			for(int k = 0; k < neighbours.size(); k++) {
				if(neighbours.get(k).getSymbol().equals("H")) numberHydrogens++; //$NON-NLS-1$
			}
			currentAtom.setImplicitHydrogenCount(numberHydrogens);
		}
		for(IAtom atom : hydrogenAtoms) {
			molecule.removeAtom(atom);
		}
	}
}