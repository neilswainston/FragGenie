package uk.ac.liverpool.metfraglib;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
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

/**
 * 
 * @author neilswainston
 */
@SuppressWarnings("boxing")
public class Precursor {
	
	/**
	 * 
	 */
	private static final Map<String, Double> MONOISOTOPIC_MASSES = new HashMap<>();

	/**
	 * 
	 */
	private final FastBitArray aromaticBonds;

	/**
	 * 
	 */
	private final int[] atomAdjacencyList;

	/**
	 * 
	 */
	private final List<int[]> atomIndexToConnectedAtomIndeces = new ArrayList<>();

	/**
	 * 
	 */
	private final int[][] bondIndexToConnectedAtomIndeces;

	/**
	 * 
	 */
	private final IAtomContainer atomContainer;

	/**
	 * 
	 */
	private final FastBitArray[] ringBondToBelongingRingBondIndeces;

	/**
	 * 
	 * @param precMolecule
	 * @throws Exception
	 */
	Precursor(final String smiles) throws Exception {
		this.atomContainer = getAtomContainer(smiles);
		this.bondIndexToConnectedAtomIndeces = new int[this.getNonHydrogenBondCount()][2];
		this.ringBondToBelongingRingBondIndeces = new FastBitArray[this.atomContainer.getBondCount() + 1];
		this.aromaticBonds = new FastBitArray(this.getNonHydrogenBondCount(), false);
		this.atomAdjacencyList = new int[getIndex(this.getNonHydrogenAtomCount() - 2, this.getNonHydrogenAtomCount() - 1) + 1];
		
		this.initiliseAtomIndexToConnectedAtomIndeces();
		this.initiliseBondIndexToConnectedAtomIndeces();
		this.initialiseRingBondsFastBitArray();
		this.initialiseAtomAdjacencyList();
	}

	/**
	 * Returns bond index + 1
	 * 
	 * @param x
	 * @param y
	 * @return int
	 */
	int getBondIndexFromAtomAdjacencyList(final int x, final int y) {
		return this.atomAdjacencyList[this.getIndex(x, y)];
	}

	/**
	 * returns atom indeces that are connected by bond with bondIndex
	 * 
	 * @param atomIndex
	 * @return int[]
	 */
	int[] getConnectedAtomIndecesOfAtomIndex(final int atomIndex) {
		return this.atomIndexToConnectedAtomIndeces.get(atomIndex);
	}

	/**
	 * 
	 * @param bondIndex
	 * @return int[]
	 */
	int[] getConnectedAtomIndecesOfBondIndex(final int bondIndex) {
		return this.bondIndexToConnectedAtomIndeces[bondIndex];
	}

	/**
	 * 
	 * @param idx
	 * @return double
	 */
	double getMassOfAtom(final int idx) {
		return MONOISOTOPIC_MASSES.get(this.getAtom(idx)).doubleValue()
				+ this.getNumberHydrogensConnectedToAtomIndex(idx)
				* MONOISOTOPIC_MASSES.get("H").doubleValue(); //$NON-NLS-1$
	}
	
	/**
	 * 
	 * @param index
	 * @return String
	 */
	String getAtom(final int index) {
		return this.atomContainer.getAtom(index).getSymbol();
	}

	/**
	 * 
	 * @return int
	 */
	int getNonHydrogenAtomCount() {
		return this.atomContainer.getAtomCount();
	}

	/**
	 * 
	 * @return int
	 */
	int getNonHydrogenBondCount() {
		return this.atomContainer.getBondCount();
	}

	/**
	 * 
	 * @param atomIndex
	 * @return int
	 */
	int getNumberHydrogensConnectedToAtomIndex(final int idx) {
		return this.atomContainer.getAtom(idx).getImplicitHydrogenCount();
	}

	/**
	 * 
	 * @return IAtomContainer
	 */
	IAtomContainer getStructureAsIAtomContainer() {
		return this.atomContainer;
	}
	
	/**
	 * 
	 * @param i
	 * @return List<Object>
	 */
	List<Object> getBond(final int i) {
		final List<Object> bondDefinition = new ArrayList<>();
		final List<String> atoms = new ArrayList<>();
		final IBond bond = this.atomContainer.getBond(i);
		
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
			final Iterator<IAtom> atoms = this.atomContainer.getBond(i).atoms().iterator();
			this.atomAdjacencyList[getIndex(this.atomContainer.indexOf(atoms.next()),
					this.atomContainer.indexOf(atoms.next()))] = (short) (i + 1);
		}
	}

	/**
	 * Initialise indeces belonging to a ring in the precursor molecule.
	 */
	private void initialiseRingBondsFastBitArray() throws Exception {
		final AllRingsFinder allRingsFinder = new AllRingsFinder();
		final IRingSet ringSet = allRingsFinder.findAllRings(this.atomContainer);

		this.initialiseRingBondToBelongingRingBondIndecesFastBitArrays(ringSet);

		if(ringSet.getAtomContainerCount() != 0) {
			final Aromaticity arom = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
			final Set<IBond> aromBonds = arom.findBonds(this.atomContainer);
			final Iterator<IBond> it = aromBonds.iterator();
			
			while(it.hasNext()) {
				this.aromaticBonds.set(this.atomContainer.indexOf(it.next()));
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
			this.ringBondToBelongingRingBondIndeces[i] = new FastBitArray(this.atomContainer.getBondCount() + 1, false);
		}
			
		for(int i = 0; i < ringSet.getAtomContainerCount(); i++) {
			final int[] indexes = new int[ringSet.getAtomContainer(i).getBondCount()];
			
			for(int j = 0; j < ringSet.getAtomContainer(i).getBondCount(); j++) {
				indexes[j] = this.atomContainer.indexOf(ringSet.getAtomContainer(i).getBond(j));
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
			final List<IAtom> connectedAtoms = this.atomContainer
					.getConnectedAtomsList(this.atomContainer.getAtom(i));
			final int[] connectedAtomIndeces = new int[connectedAtoms.size()];
			
			for(int k = 0; k < connectedAtoms.size(); k++) {
				connectedAtomIndeces[k] = this.atomContainer.indexOf(connectedAtoms.get(k));
			}
				
			this.atomIndexToConnectedAtomIndeces.add(i, connectedAtomIndeces);
		}
	}

	/**
	 * 
	 */
	private void initiliseBondIndexToConnectedAtomIndeces() {
		for (int i = 0; i < this.getNonHydrogenBondCount(); i++) {
			final IBond bond = this.atomContainer.getBond(i);
			this.bondIndexToConnectedAtomIndeces[i][0] = this.atomContainer.indexOf(bond.getAtom(0));
			this.bondIndexToConnectedAtomIndeces[i][1] = this.atomContainer.indexOf(bond.getAtom(1));
		}
	}
	
	/**
	 * 
	 * @param smiles
	 * @return IAtomContainer
	 * @throws Exception
	 */
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
	
	static {
		MONOISOTOPIC_MASSES.put("[13C]", 13.00335); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("C", 12.00000); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Al", 26.98154); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Am", 243.06140); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ar", 39.96238); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("As", 74.92160); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("At", 209.98710); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Au", 196.96660); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("B", 11.00931); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ba", 137.90520); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Bi", 208.98040); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Br", 78.91834); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ca", 39.96259); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Cd", 113.90340); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ce", 139.90540); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Cl", 34.96885); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Co", 58.93320); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Cr", 51.94051);  //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Cu", 62.92960); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("D", 2.01410); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("F", 18.99840); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Fe", 55.93494); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ga", 68.92558); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Gd", 157.92410); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ge", 73.92118); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("H", 1.00783); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("He", 4.00260); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Hg", 201.97060); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("I", 126.90450); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("In", 114.90390); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("K", 38.96371); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Li", 7.01600); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Mg", 23.98504); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Mn", 54.93805); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Mo", 97.90541); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("[15N]", 15.00011);  //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("N", 14.00307); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Na", 22.98977); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ne", 19.99244); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ni", 57.93535); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("[18O]", 17.99916); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("O", 15.99491); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("P", 30.97376); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Pb", 207.97660); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Po", 208.98240); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Pt", 194.96480); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ra", 226.02540); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Rb", 84.91179); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Re", 186.95580); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ru", 101.90430); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("S", 31.97207); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Sb", 120.90380); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Se", 79.91652); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Si", 27.97693); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Sn", 119.90220); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Tc", 97.90722); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Te", 129.90620); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ti", 47.94795); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Tl", 204.97440); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("U", 238.05080); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("V", 50.94396); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("W", 183.95090); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Y", 88.90585); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Zn", 63.92915); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Zr", 89.90470);  //$NON-NLS-1$
	}
}