package uk.ac.liverpool.metfraglib;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import uk.ac.liverpool.metfraglib.FastBitArray;
import uk.ac.liverpool.metfraglib.Constants;

/**
 * 
 * @author neilswainston
 */
public class Precursor {

	/**
	 * 
	 * @param smiles
	 * @return Precursor
	 * @throws Exception
	 */
	public static Precursor fromSmiles(final String smiles) throws Exception {
		final IAtomContainer molecule = getAtomContainer(smiles);
		return new Precursor(molecule);
	}
	
	private static IAtomContainer getAtomContainer(final String smiles) throws Exception {
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
	
	/**
	 * 
	 * @param molecule
	 */
	private static void removeHydrogens(final IAtomContainer molecule) {
		final Collection<IAtom> hydrogenAtoms = new ArrayList<>();
		
		for(IAtom atom : molecule.atoms()) {
			if(atom.getSymbol().equals("H")) { //$NON-NLS-1$
				hydrogenAtoms.add(atom);
			}
			
			int numberHydrogens = 0;
			
			for(IAtom neighbour : molecule.getConnectedAtomsList(atom)) {
				if(neighbour.getSymbol().equals("H")) { //$NON-NLS-1$
					numberHydrogens++; 
				}
			}
			
			atom.setImplicitHydrogenCount(Integer.valueOf(numberHydrogens));
		}
		
		for(IAtom atom : hydrogenAtoms) {
			molecule.removeAtom(atom);
		}
	}

	/**
	 * 
	 */
	private final FastBitArray aromaticBonds;

	/**
	 * 
	 */
	private final short[] atomAdjacencyList;

	/**
	 * 
	 */
	private final List<short[]> atomIndexToConnectedAtomIndeces = new ArrayList<>();

	/**
	 * 
	 */
	private final short[][] bondIndexToConnectedAtomIndeces;
	
	/**
	 * 
	 */
	private final double[] massesOfAtoms;

	/**
	 * 
	 */
	private final Integer[] numberHydrogensConnectedToAtom;

	/**
	 * 
	 */
	private final IAtomContainer precursorMolecule;

	/**
	 * 
	 */
	private final FastBitArray[] ringBondToBelongingRingBondIndeces;

	/**
	 * 
	 * @param precMolecule
	 * @throws Exception
	 */
	private Precursor(final IAtomContainer precMolecule) throws Exception {
		this.precursorMolecule = precMolecule;
		this.bondIndexToConnectedAtomIndeces = new short[this.getNonHydrogenBondCount()][2];
		this.ringBondToBelongingRingBondIndeces = new FastBitArray[this.precursorMolecule.getBondCount() + 1];
		this.aromaticBonds = new FastBitArray(this.getNonHydrogenBondCount());
		this.atomAdjacencyList = new short[getIndex(this.getNonHydrogenAtomCount() - 2,
				this.getNonHydrogenAtomCount() - 1) + 1];
		this.numberHydrogensConnectedToAtom = new Integer[this.getNonHydrogenAtomCount()];
		this.massesOfAtoms = new double[this.getNonHydrogenAtomCount()];
		
		this.initiliseAtomIndexToConnectedAtomIndeces();
		this.initialiseNumberHydrogens();
		this.initiliseBondIndexToConnectedAtomIndeces();
		this.initialiseRingBondsFastBitArray();
		this.initialiseAtomAdjacencyList();
		this.initialiseAtomMasses();
	}

	/**
	 * Returns bond index + 1
	 * 
	 * @param x
	 * @param y
	 * @return short
	 */
	public short getBondIndexFromAtomAdjacencyList(final short x, final short y) {
		return this.atomAdjacencyList[this.getIndex(x, y)];
	}

	/**
	 * returns atom indeces that are connected by bond with bondIndex
	 * 
	 * @param atomIndex
	 * @return short[]
	 */
	public short[] getConnectedAtomIndecesOfAtomIndex(final short atomIndex) {
		return this.atomIndexToConnectedAtomIndeces.get(atomIndex);
	}

	/**
	 * 
	 * @param bondIndex
	 * @return short[]
	 */
	public short[] getConnectedAtomIndecesOfBondIndex(final short bondIndex) {
		return this.bondIndexToConnectedAtomIndeces[bondIndex];
	}

	/**
	 * 
	 * @param index
	 * @return double
	 */
	public double getMassOfAtom(final int index) {
		return this.massesOfAtoms[index] + this.getNumberHydrogensConnectedToAtomIndex(index)
				* Constants.MONOISOTOPIC_MASSES.get("H").doubleValue(); //$NON-NLS-1$
	}
	
	/**
	 * 
	 * @param index
	 * @return String
	 */
	public String getAtom(final int index) {
		return this.precursorMolecule.getAtom(index).getSymbol();
	}

	/**
	 * 
	 * @return int
	 */
	public int getNonHydrogenAtomCount() {
		return this.precursorMolecule.getAtomCount();
	}

	/**
	 * 
	 * @return int
	 */
	public int getNonHydrogenBondCount() {
		return this.precursorMolecule.getBondCount();
	}

	/**
	 * 
	 * @param atomIndex
	 * @return
	 */
	public int getNumberHydrogensConnectedToAtomIndex(final int atomIndex) {
		return this.numberHydrogensConnectedToAtom[atomIndex].intValue();
	}

	/**
	 * 
	 * @return IAtomContainer
	 */
	public IAtomContainer getStructureAsIAtomContainer() {
		return this.precursorMolecule;
	}
	
	/**
	 * 
	 * @param idx
	 * @return List<Object>
	 */
	public List<Object> getBond(final short idx) {
		final List<Object> bondDefinition = new ArrayList<>();
		final List<String> atoms = new ArrayList<>();
		final IBond bond = this.precursorMolecule.getBond(idx);
		
		for(IAtom atom : bond.atoms()) {
			atoms.add(atom.getSymbol());
		}
		
		Collections.sort(atoms);
		
		bondDefinition.add(atoms);
		bondDefinition.add(bond.getOrder());
		bondDefinition.add(Boolean.valueOf(bond.isAromatic()));
		return bondDefinition;
	}

	/**
	 * Convert 2D matrix coordinates to 1D adjacency list coordinate.
	 * 
	 * @param a
	 * @param b
	 * @return int
	 */
	private int getIndex(final int a, final int b) {
		int row = a;
		int col = b;
		
		if(a > b) {
			row = b;
			col = a;
		}
		
		return row * this.getNonHydrogenAtomCount() + col - ((row + 1) * (row + 2)) / 2;
	}

	/**
	 * Initialise 1D atom adjacency list.
	 */
	private void initialiseAtomAdjacencyList() {
		for(int i = 0; i < this.getNonHydrogenBondCount(); i++) {
			final Iterator<IAtom> atoms = this.precursorMolecule.getBond(i).atoms().iterator();
			this.atomAdjacencyList[getIndex(this.precursorMolecule.indexOf(atoms.next()),
					this.precursorMolecule.indexOf(atoms.next()))] = (short) (i + 1);
		}
	}

	/**
	 * 
	 */
	private void initialiseAtomMasses() {
		for(int i = 0; i < this.getNonHydrogenAtomCount(); i++) {
			this.massesOfAtoms[i] = Constants.MONOISOTOPIC_MASSES
					.get(this.precursorMolecule.getAtom(i).getSymbol()).doubleValue();
		}
	}

	/**
	 * 
	 */
	private void initialiseNumberHydrogens() {
		for(int i = 0; i < this.getNonHydrogenAtomCount(); i++) {
			this.numberHydrogensConnectedToAtom[i] = this.precursorMolecule.getAtom(i).getImplicitHydrogenCount();
		}
	}

	/**
	 * Initialise indeces belonging to a ring in the precursor molecule.
	 */
	private void initialiseRingBondsFastBitArray() throws Exception {
		final AllRingsFinder allRingsFinder = new AllRingsFinder();
		final IRingSet ringSet = allRingsFinder.findAllRings(this.precursorMolecule);

		this.initialiseRingBondToBelongingRingBondIndecesFastBitArrays(ringSet);

		if(ringSet.getAtomContainerCount() != 0) {
			final Aromaticity arom = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
			final Set<IBond> aromBonds = arom.findBonds(this.precursorMolecule);
			final Iterator<IBond> it = aromBonds.iterator();
			
			while(it.hasNext()) {
				this.aromaticBonds.set(this.precursorMolecule.indexOf(it.next()));
			}
		}
	}

	/**
	 * initialises ringBondToBelongingRingBondIndeces FastBitArray array fast and
	 * easy way to retrieve all bond indices belonging to a ring including the bond
	 * at specified index of that array.
	 * 
	 * @param ringSet
	 */
	private void initialiseRingBondToBelongingRingBondIndecesFastBitArrays(final IRingSet ringSet) {
		for(int i = 0; i < this.ringBondToBelongingRingBondIndeces.length; i++) {
			this.ringBondToBelongingRingBondIndeces[i] = new FastBitArray(this.precursorMolecule.getBondCount() + 1);
		}
			
		for(int i = 0; i < ringSet.getAtomContainerCount(); i++) {
			final int[] indexes = new int[ringSet.getAtomContainer(i).getBondCount()];
			
			for(int j = 0; j < ringSet.getAtomContainer(i).getBondCount(); j++) {
				indexes[j] = this.precursorMolecule.indexOf(ringSet.getAtomContainer(i).getBond(j));
			}
			
			for(int j = 0; j < indexes.length; j++) {
				this.ringBondToBelongingRingBondIndeces[indexes[j]].setBits(indexes);
			}
		}
	}

	/**
	 * 
	 */
	private void initiliseAtomIndexToConnectedAtomIndeces() {
		for (int i = 0; i < this.getNonHydrogenAtomCount(); i++) {
			final List<IAtom> connectedAtoms = this.precursorMolecule
					.getConnectedAtomsList(this.precursorMolecule.getAtom(i));
			final short[] connectedAtomIndeces = new short[connectedAtoms.size()];
			
			for(int k = 0; k < connectedAtoms.size(); k++) {
				connectedAtomIndeces[k] = (short) this.precursorMolecule.indexOf(connectedAtoms.get(k));
			}
				
			this.atomIndexToConnectedAtomIndeces.add(i, connectedAtomIndeces);
		}
	}

	/**
	 * 
	 */
	private void initiliseBondIndexToConnectedAtomIndeces() {
		for (int i = 0; i < this.getNonHydrogenBondCount(); i++) {
			this.bondIndexToConnectedAtomIndeces[i][0] = (short) this.precursorMolecule
					.indexOf(this.precursorMolecule.getBond(i).getAtom(0));
			this.bondIndexToConnectedAtomIndeces[i][1] = (short) this.precursorMolecule
					.indexOf(this.precursorMolecule.getBond(i).getAtom(1));
		}
	}
}