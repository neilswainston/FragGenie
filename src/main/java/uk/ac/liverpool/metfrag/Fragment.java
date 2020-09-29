package uk.ac.liverpool.metfrag;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.TreeMap;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 * 
 * @author neilswainston
 */
@SuppressWarnings("boxing")
public class Fragment implements Comparable<Fragment> {

	/**
	 * 
	 */
	private static final Map<String, Double> MONOISOTOPIC_MASSES = new HashMap<>();
	
	/**
	 * 
	 */
	private final IAtomContainer prec;

	/**
	 * 
	 */
	private boolean[] atomsArray;

	/**
	 * 
	 */
	private boolean[] bondsArray;

	/**
	 * 
	 */
	private boolean[] brokenBondsArray;

	/**
	 * 
	 */
	private int addedToQueueCounts = 0;

	/**
	 * 
	 */
	private int lastSkippedBond = -1;

	/**
	 * 
	 */
	private int treeDepth = 0;

	/**
	 * Constructor.
	 * 
	 * @param precursor
	 */
	Fragment(final IAtomContainer precursor) {
		this(precursor, new boolean[precursor.getAtomCount()],
				new boolean[precursor.getBondCount()],
				new boolean[precursor.getBondCount()]);
		
		Arrays.fill(this.atomsArray, true);
		Arrays.fill(this.bondsArray, true);
	}

	/**
	 * Constructor.
	 * 
	 * @param precursor
	 * @param atoms
	 * @param bonds
	 * @param brokenBonds
	 */
	private Fragment(final IAtomContainer precursor, final boolean[] atoms, final boolean[] bonds,
			final boolean[] brokenBonds) {
		this.prec = precursor;
		this.atomsArray = atoms;
		this.bondsArray = bonds;
		this.brokenBondsArray = brokenBonds;
	}
	
	@Override
	public String toString() {
		return getFormula();
	}
	
	/**
	 * 
	 * @return float
	 */
	public float getMonoisotopicMass() {
		float mass = 0.0f;

		for (int i = 0; i < this.atomsArray.length; i++) {
			if (this.atomsArray[i]) {
				mass += this.getAtomMass(i);
			}
		}
		return mass;
	}

	/**
	 * 
	 * @return String
	 */
	public String getFormula() {
		final String HYDROGEN = "H"; //$NON-NLS-1$
		final Map<String, Integer> elementCount = new TreeMap<>();

		for (int i = 0; i < this.atomsArray.length; i++) {
			if (this.atomsArray[i]) {
				final IAtom atom = this.prec.getAtom(i);
				final String element = atom.getSymbol();

				if (elementCount.get(element) == null) {
					elementCount.put(element, Integer.valueOf(1));
				} else {
					elementCount.put(element, Integer.valueOf(elementCount.get(element).intValue() + 1));
				}

				final int hCount = atom.getImplicitHydrogenCount();

				if (elementCount.get(HYDROGEN) == null) {
					elementCount.put(HYDROGEN, Integer.valueOf(hCount));
				} else {
					elementCount.put(HYDROGEN, Integer.valueOf(elementCount.get(HYDROGEN).intValue() + hCount));
				}
			}
		}

		final StringBuilder builder = new StringBuilder();

		for (Map.Entry<String, Integer> entry : elementCount.entrySet()) {
			builder.append(entry.getKey());
			final int count = entry.getValue().intValue();

			if (count > 1) {
				builder.append(count);
			}
		}

		return builder.toString();
	}
	
	/**
	 * 
	 * @return Collection<Bond>
	 */
	public Collection<Bond> getBrokenBonds() {
		final Collection<Bond> bonds = new ArrayList<>();
		
		for (int i = 0; i < this.brokenBondsArray.length; i++) {
			if(this.brokenBondsArray[i]) {
				final IBond iBond = this.prec.getBond(i);
				final Iterator<IAtom> atoms = iBond.atoms().iterator();
				bonds.add(new Bond(atoms.next(), atoms.next(), iBond.getOrder(), iBond.isAromatic()));
			}
		}
		
		return bonds;
	}

	/**
	 * 
	 * @return int
	 */
	int getAddedToQueueCounts() {
		return this.addedToQueueCounts;
	}

	/**
	 * 
	 * @return boolean[]
	 */
	boolean[] getAtomsArray() {
		return this.atomsArray;
	}

	/**
	 * 
	 * @return boolean[]
	 */
	boolean[] getBondsArray() {
		return this.bondsArray;
	}

	/**
	 * 
	 * @return boolean[]
	 */
	boolean[] getBrokenBondsArray() {
		return this.brokenBondsArray;
	}

	/**
	 * 
	 * @return int
	 */
	int getLastSkippedBond() {
		return this.lastSkippedBond;
	}

	/**
	 * 
	 * @param i
	 */
	void setAddedToQueueCounts(final int i) {
		this.addedToQueueCounts = i;
	}

	/**
	 * 
	 * @param i
	 */
	void setLastSkippedBond(final int i) {
		this.lastSkippedBond = i;
	}

	/**
	 * 
	 * @param precursorMolecule
	 * @param removeBond
	 * @param bondConnectedAtoms
	 * @return Fragment[]
	 * @throws Exception
	 */
	Fragment[] fragment(final IAtomContainer precursorMolecule, final int removeBond, final int[] bondConnectedAtoms) throws Exception {

		// Generate first fragment:
		// Traverse to first direction from atomIndex connected by broken bond:
		final Object[] result1 = this.traverse(precursorMolecule, bondConnectedAtoms[0], bondConnectedAtoms[1],
				removeBond, this.brokenBondsArray.clone());

		final Fragment fragment1 = (Fragment)result1[1];

		// Only one fragment is generated when a ring bond was broken:
		if (((Boolean)result1[0]).booleanValue()) {
			fragment1.treeDepth = this.treeDepth;
			fragment1.addedToQueueCounts = this.addedToQueueCounts + 1;
			return new Fragment[] { fragment1 };
		}
		
		fragment1.treeDepth = this.treeDepth + 1;
		
		// Generate second fragment:
		// Traverse the second direction from atomIndex connected by broken bond:
		final Object[] result2 = this.traverse(precursorMolecule, bondConnectedAtoms[1], bondConnectedAtoms[0],
				removeBond, this.brokenBondsArray.clone());

		final Fragment fragment2 = (Fragment)result2[1];
		fragment2.treeDepth = this.treeDepth + 1;

		return new Fragment[] { fragment1, fragment2 };
	}
	
	/**
	 * 
	 * @param atomIdx
	 * @return double
	 */
	private double getAtomMass(final int atomIdx) {
		final IAtom atom = this.prec.getAtom(atomIdx);
		final String symbol = atom.getSymbol();
		final int hCount = atom.getImplicitHydrogenCount().intValue();
		
		return MONOISOTOPIC_MASSES.get(symbol).doubleValue()
				+ hCount * MONOISOTOPIC_MASSES.get("H").doubleValue(); //$NON-NLS-1$
	}

	/**
	 * Traverse the fragment to one direction starting from startAtom.
	 * 
	 * @param precursorMolecule
	 * @param startAtom
	 * @param endAtom
	 * @param removeBond
	 * @param brokenBondArrayOfNewFragment
	 * @return Object[]
	 */
	private Object[] traverse(final IAtomContainer precursorMolecule, final int startAtom, final int endAtom,
			final int removeBond, final boolean[] brokenBondArrayOfNewFragment) {
		final boolean[] newAtomArray = new boolean[precursorMolecule.getAtomCount()];
		final boolean[] newBondArray = new boolean[precursorMolecule.getBondCount()];
		final boolean[] currentBondArray = this.bondsArray;
		
		// When traversing the fragment graph, we want to know if we already visited
		// a node (atom) to check for ringed structures.
		// If traversed an already visited atom, then no new fragment was generated.
		final boolean[] visited = new boolean[precursorMolecule.getAtomCount()];

		// Traverse molecule in the first direction:
		final Stack<int[]> toProcessConnectedAtoms = new Stack<>();
		final Stack<Integer> toProcessAtom = new Stack<>();
		
		toProcessConnectedAtoms.push(this.getConnectedAtomIndicesFromAtomIdx(startAtom));
		toProcessAtom.push(Integer.valueOf(startAtom));
		visited[startAtom] = true;
		
		boolean singleFragment = false;
		
		// Set the first atom of possible new fragment atom as the one direction of cut bond.
		newAtomArray[startAtom] = true;
		
		while (!toProcessConnectedAtoms.isEmpty()) {
			final int midAtom = toProcessAtom.pop().intValue();
			
			for (int nextAtom : toProcessConnectedAtoms.pop()) {
				// Did we visit the current atom already?
				final int currentBond = this.getBondIndex(nextAtom, midAtom);

				if (!currentBondArray[currentBond] || currentBond == removeBond) {
					continue;
				}
				
				// If we visited the current atom already, then we do not have to check it again:
				if (visited[nextAtom]) {
					newBondArray[currentBond] = true;
					continue;
				}
				
				// If we reach the second atom of the cleaved bond then still one fragment is present:
				if (nextAtom == endAtom) {
					singleFragment = true;
				}

				visited[nextAtom] = true;
				newAtomArray[nextAtom] = true;

				newBondArray[this.getBondIndex(midAtom, nextAtom)] = true;
				toProcessConnectedAtoms.push(this.getConnectedAtomIndicesFromAtomIdx(nextAtom));
				toProcessAtom.push(Integer.valueOf(nextAtom));
			}
		}

		brokenBondArrayOfNewFragment[removeBond] = true;
		newBondArray[removeBond] = false;
		
		final Fragment newFragment = new Fragment(precursorMolecule, newAtomArray, newBondArray, brokenBondArrayOfNewFragment);

		return new Object[] {Boolean.valueOf(singleFragment), newFragment};
	}
	
	/**
	 * Returns bond index + 1.
	 * 
	 * @param atomIdx1
	 * @param atomIdx2
	 * @return int
	 */
	private int getBondIndex(final int atomIdx1, final int atomIdx2) {
		final int[] queryAtomIdxs = new int[] {atomIdx1, atomIdx2};
		Arrays.sort(queryAtomIdxs);
		
		for (int i = 0; i < this.prec.getBondCount(); i++) {
			final Iterator<IAtom> atoms = this.prec.getBond(i).atoms().iterator();
			final int[] bondAtomIdxs = new int[] {this.prec.indexOf(atoms.next()), this.prec.indexOf(atoms.next())};
			Arrays.sort(bondAtomIdxs);
			
			if(Arrays.equals(queryAtomIdxs, bondAtomIdxs)) {
				return i;
			}
		}
		
		return -1;
	}
	
	/**
	 * Returns atom indices that are connected to atom with atomIndex.
	 * 
	 * @param atomIndex
	 * @return int[]
	 */
	private int[] getConnectedAtomIndicesFromAtomIdx(final int atomIdx) {
		final List<IAtom> connected = this.prec.getConnectedAtomsList(this.prec.getAtom(atomIdx));
		final int[] connectedIdxs = new int[connected.size()];

		for (int k = 0; k < connected.size(); k++) {
			connectedIdxs[k] = this.prec.indexOf(connected.get(k));
		}

		return connectedIdxs;
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((this.atomsArray == null) ? 0 : this.atomsArray.hashCode());
		result = prime * result + ((this.bondsArray == null) ? 0 : this.bondsArray.hashCode());
		result = prime * result + ((this.brokenBondsArray == null) ? 0 : this.brokenBondsArray.hashCode());
		result = prime * result + ((this.prec == null) ? 0 : this.prec.hashCode());
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
			
		
		final Fragment other = (Fragment) obj;
		
		if (this.atomsArray == null) {
			if (other.atomsArray != null) {
				return false;
			}	
		}
		else if (!this.atomsArray.equals(other.atomsArray)) {
			return false;
		}
			
		
		if (this.bondsArray == null) {
			if (other.bondsArray != null) {
				return false;
			}	
		}
		else if (!this.bondsArray.equals(other.bondsArray)) {
			return false;
		}
			
		if (this.brokenBondsArray == null) {
			if (other.brokenBondsArray != null) {
				return false;
			}
		}
		else if (!this.brokenBondsArray.equals(other.brokenBondsArray)) {
			return false;
		}
			
		if (this.prec == null) {
			if (other.prec != null) {
				return false;
			}
				
		}
		else if (!this.prec.equals(other.prec)) {
			return false;
		}
			
		return true;
	}

	@Override
	public int compareTo(final Fragment obj) {
		return (int)((this.getMonoisotopicMass() - obj.getMonoisotopicMass()) * 1000000);
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