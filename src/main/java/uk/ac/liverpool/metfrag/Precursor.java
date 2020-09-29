package uk.ac.liverpool.metfrag;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
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
	private final IAtomContainer atomContainer;
	
	/**
	 * 
	 */
	private final FastBitArray aromaticBonds;

	/**
	 * 
	 */
	private final int[] atomAdjacencies;

	/**
	 * 
	 * @param smiles
	 * @throws CDKException
	 */
	Precursor(final String smiles) throws CDKException {
		this.atomContainer = getAtomContainer(smiles);
		
		this.aromaticBonds = new FastBitArray(this.getBondCount(), false);
		this.initialiseAromaticBonds();
		
		// Initialise atomAdjacencies:
		this.atomAdjacencies = new int[getIndex(this.getAtomCount() - 2, this.getAtomCount() - 1) + 1];
		
		for (int i = 0; i < this.getBondCount(); i++) {
			final Iterator<IAtom> atoms = this.atomContainer.getBond(i).atoms().iterator();
			this.atomAdjacencies[getIndex(this.atomContainer.indexOf(atoms.next()), this.atomContainer.indexOf(atoms.next()))] = i + 1;
		}
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
	 * Returns bond index + 1.
	 * 
	 * @param atomIdx1
	 * @param atomIdx2
	 * @return int
	 */
	int getBondIndex(final int atomIdx1, final int atomIdx2) {
		return this.atomAdjacencies[this.getIndex(atomIdx1, atomIdx2)];
	}

	/**
	 * Returns atom indices that are connected to atom with atomIndex.
	 * 
	 * @param atomIndex
	 * @return int[]
	 */
	int[] getConnectedAtomIndicesFromAtomIdx(final int atomIdx) {
		final List<IAtom> connected = this.atomContainer.getConnectedAtomsList(this.atomContainer.getAtom(atomIdx));
		final int[] connectedIdxs = new int[connected.size()];

		for (int k = 0; k < connected.size(); k++) {
			connectedIdxs[k] = this.atomContainer.indexOf(connected.get(k));
		}

		return connectedIdxs;
	}

	/**
	 * Returns atom indices that are connected to bond with bondIdx.
	 * 
	 * @param bondIdx
	 * @return int[]
	 */
	int[] getConnectedAtomIndecesFromBondIndex(final int bondIdx) {
		final IBond bond = this.atomContainer.getBond(bondIdx);
		return new int[] { this.atomContainer.indexOf(bond.getAtom(0)), this.atomContainer.indexOf(bond.getAtom(1)) };
	}

	/**
	 * 
	 * @param idx
	 * @return double
	 */
	double getMassOfAtom(final int idx) {
		return MONOISOTOPIC_MASSES.get(this.getAtomSymbol(idx)).doubleValue()
				+ this.getImplicitHydrogenCount(idx) * MONOISOTOPIC_MASSES.get("H").doubleValue(); //$NON-NLS-1$
	}

	/**
	 * 
	 * @param atomIdx
	 * @return String
	 */
	String getAtomSymbol(final int atomIdx) {
		return this.atomContainer.getAtom(atomIdx).getSymbol();
	}

	/**
	 * 
	 * @return int
	 */
	int getAtomCount() {
		return this.atomContainer.getAtomCount();
	}

	/**
	 * 
	 * @return int
	 */
	int getBondCount() {
		return this.atomContainer.getBondCount();
	}

	/**
	 * 
	 * @param atomIdx
	 * @return int
	 */
	int getImplicitHydrogenCount(final int atomIdx) {
		return this.atomContainer.getAtom(atomIdx).getImplicitHydrogenCount();
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
	 * @param bondIdx
	 * @return Bond
	 */
	Bond getBond(final int bondIdx) {
		final IBond iBond = this.atomContainer.getBond(bondIdx);
		final Iterator<IAtom> atoms = iBond.atoms().iterator();
		return new Bond(atoms.next(), atoms.next(), iBond.getOrder(), iBond.isAromatic());
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

		if (a > b) {
			row = b;
			col = a;
		}

		return row * this.getAtomCount() + col - ((row + 1) * (row + 2)) / 2;
	}

	/**
	 * Initialise indices belonging to a ring in the precursor molecule.
	 * 
	 * @throws CDKException 
	 */
	private void initialiseAromaticBonds() throws CDKException {
		final IRingSet ringSet = new AllRingsFinder().findAllRings(this.atomContainer);

		if (ringSet.getAtomContainerCount() != 0) {
			final Aromaticity arom = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
			final Set<IBond> aromBonds = arom.findBonds(this.atomContainer);
			final Iterator<IBond> it = aromBonds.iterator();

			while (it.hasNext()) {
				this.aromaticBonds.set(this.atomContainer.indexOf(it.next()));
			}
		}
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
		MONOISOTOPIC_MASSES.put("Cr", 51.94051); //$NON-NLS-1$
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
		MONOISOTOPIC_MASSES.put("[15N]", 15.00011); //$NON-NLS-1$
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
		MONOISOTOPIC_MASSES.put("Zr", 89.90470); //$NON-NLS-1$
	}
}