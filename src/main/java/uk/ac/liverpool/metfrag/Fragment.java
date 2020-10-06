package uk.ac.liverpool.metfrag;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 * 
 * @author neilswainston
 */
public class Fragment implements Comparable<Fragment> {

	/**
	 * 
	 */
	private static final Map<String, Float> MONOISOTOPIC_MASSES = new HashMap<>();
	
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
	private final int[][] bondAtomIdxs;

	/**
	 * 
	 */
	private float monoisotopicMass = Float.NaN;

	/**
	 * Constructor.
	 * 
	 * @param precursor
	 */
	Fragment(final IAtomContainer precursor) {
		this(precursor,
				new boolean[precursor.getAtomCount()],
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
	private Fragment(final IAtomContainer precursor, final boolean[] atoms, final boolean[] bonds, final boolean[] brokenBonds) {
		this.prec = precursor;
		this.atomsArray = atoms;
		this.bondsArray = bonds;
		this.brokenBondsArray = brokenBonds;
		this.bondAtomIdxs = new int[atoms.length][atoms.length];
		
		for (int i = 0; i < this.prec.getBondCount(); i++) {
			final Iterator<IAtom> atomsIterator = this.prec.getBond(i).atoms().iterator();
			final int idx1 = this.prec.indexOf(atomsIterator.next());
			final int idx2 = this.prec.indexOf(atomsIterator.next());
			this.bondAtomIdxs[idx1][idx2] = i;
			this.bondAtomIdxs[idx2][idx1] = i;
		}
	}
	
	/**
	 * 
	 * @return Collection<List<Object>>
	 */
	public Collection<Fragment> fragment(final int maxBrokenBonds) {
		final Collection<Fragment> fragments = new TreeSet<>();
		fragment(this, fragments, maxBrokenBonds);
		return fragments;
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
		if(Float.isNaN(this.monoisotopicMass)) {
			float mass = 0.0f;

			for (int i = 0; i < this.atomsArray.length; i++) {
				if (this.atomsArray[i]) {
					mass += this.getAtomMass(i);
				}
			}
			
			this.monoisotopicMass = mass;
		}
		
		return this.monoisotopicMass;
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

				final int hCount = atom.getImplicitHydrogenCount().intValue();

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
		final List<Bond> bonds = new ArrayList<>();
		
		for (int i = 0; i < this.brokenBondsArray.length; i++) {
			if(this.brokenBondsArray[i]) {
				final IBond iBond = this.prec.getBond(i);
				final Iterator<IAtom> atoms = iBond.atoms().iterator();
				bonds.add(new Bond(atoms.next(), atoms.next(), iBond.getOrder(), iBond.isAromatic()));
			}
		}
		
		Collections.sort(bonds);
		
		return bonds;
	}
	
	/**
	 * 
	 * @return int
	 */
	int getNumBrokenBonds() {
		int numBrokenBonds = 0;
		
		for (int i = 0; i < this.brokenBondsArray.length; i++) {
			if(this.brokenBondsArray[i]) {
				numBrokenBonds++;
			}
		}
		
		return numBrokenBonds;
	}
	
	/**
	 * 
	 * @param fragment
	 * @param fragments
	 * @param maxBrokenBonds
	 */
	private static void fragment(final Fragment fragment, final Collection<Fragment> fragments, final int maxBrokenBonds) {
		if(fragment.getNumBrokenBonds() <= maxBrokenBonds) {
			fragments.add(fragment);
			
			if(fragment.getNumBrokenBonds() < maxBrokenBonds) {
				for(int bondIdx = 0; bondIdx < fragment.bondsArray.length; bondIdx++) {
					if(fragment.bondsArray[bondIdx] && !fragment.brokenBondsArray[bondIdx]) {
						for(final Fragment childFragment : fragment.fragmentBond(bondIdx)) {
							fragment(childFragment, fragments, maxBrokenBonds);
						}
					}
				}
			}
		}
	}

	/**
	 * 
	 * @param bondIdx
	 * @return Fragment[]
	 */
	Fragment[] fragmentBond(final int bondIdx) {
		final IBond bond = this.prec.getBond(bondIdx);
		final int[] bondConnectedAtoms =  new int[] { this.prec.indexOf(bond.getAtom(0)), this.prec.indexOf(bond.getAtom(1)) };
		
		// Generate first fragment:
		// Traverse to first direction from atomIndex connected by broken bond:
		final Object[] result1 = this.traverse(bondConnectedAtoms[0], bondConnectedAtoms[1],
				bondIdx, this.brokenBondsArray.clone());

		final Fragment fragment1 = (Fragment)result1[1];

		// Only one fragment is generated when a ring bond was broken:
		if (((Boolean)result1[0]).booleanValue()) {
			return new Fragment[] { fragment1 };
		}
		
		// Generate second fragment:
		// Traverse the second direction from atomIndex connected by broken bond:
		final Object[] result2 = this.traverse(bondConnectedAtoms[1], bondConnectedAtoms[0],
				bondIdx, this.brokenBondsArray.clone());

		final Fragment fragment2 = (Fragment)result2[1];

		return new Fragment[] { fragment1, fragment2 };
	}
	
	/**
	 * 
	 * @param atomIdx
	 * @return float
	 */
	private float getAtomMass(final int atomIdx) {
		final IAtom atom = this.prec.getAtom(atomIdx);
		final String symbol = atom.getSymbol();
		final int hCount = atom.getImplicitHydrogenCount().intValue();
		
		return MONOISOTOPIC_MASSES.get(symbol).floatValue()
				+ hCount * MONOISOTOPIC_MASSES.get("H").floatValue(); //$NON-NLS-1$
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
	private Object[] traverse(final int startAtom, final int endAtom,
			final int removeBond, final boolean[] brokenBondArrayOfNewFragment) {
		final boolean[] newAtomArray = new boolean[this.prec.getAtomCount()];
		final boolean[] newBondArray = new boolean[this.prec.getBondCount()];
		final boolean[] currentBondArray = this.bondsArray;
		
		// When traversing the fragment graph, we want to know if we already visited
		// a node (atom) to check for ringed structures.
		// If traversed an already visited atom, then no new fragment was generated.
		final boolean[] visited = new boolean[this.prec.getAtomCount()];

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
				final int currentBond = this.bondAtomIdxs[midAtom][nextAtom];

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

				newBondArray[this.bondAtomIdxs[midAtom][nextAtom]] = true;
				toProcessConnectedAtoms.push(this.getConnectedAtomIndicesFromAtomIdx(nextAtom));
				toProcessAtom.push(Integer.valueOf(nextAtom));
			}
		}

		brokenBondArrayOfNewFragment[removeBond] = true;
		newBondArray[removeBond] = false;
		
		final Fragment newFragment = new Fragment(this.prec, newAtomArray, newBondArray, brokenBondArrayOfNewFragment);

		return new Object[] {Boolean.valueOf(singleFragment), newFragment};
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
		MONOISOTOPIC_MASSES.put("[13C]", Float.valueOf(13.00335f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("C", Float.valueOf(12.00000f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Al", Float.valueOf(26.98154f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Am", Float.valueOf(243.06140f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ar", Float.valueOf(39.96238f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("As", Float.valueOf(74.92160f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("At", Float.valueOf(209.98710f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Au", Float.valueOf(196.96660f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("B", Float.valueOf(11.00931f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ba", Float.valueOf(137.90520f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Bi", Float.valueOf(208.98040f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Br", Float.valueOf(78.91834f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ca", Float.valueOf(39.96259f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Cd", Float.valueOf(113.90340f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ce", Float.valueOf(139.90540f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Cl", Float.valueOf(34.96885f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Co", Float.valueOf(58.93320f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Cr", Float.valueOf(51.94051f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Cu", Float.valueOf(62.92960f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("D", Float.valueOf(2.01410f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("F", Float.valueOf(18.99840f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Fe", Float.valueOf(55.93494f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ga", Float.valueOf(68.92558f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Gd", Float.valueOf(157.92410f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ge", Float.valueOf(73.92118f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("H", Float.valueOf(1.00783f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("He", Float.valueOf(4.00260f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Hg", Float.valueOf(201.97060f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("I", Float.valueOf(126.90450f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("In", Float.valueOf(114.90390f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("K", Float.valueOf(38.96371f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Li", Float.valueOf(7.01600f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Mg", Float.valueOf(23.98504f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Mn", Float.valueOf(54.93805f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Mo", Float.valueOf(97.90541f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("[15N]", Float.valueOf(15.00011f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("N", Float.valueOf(14.00307f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Na", Float.valueOf(22.98977f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ne", Float.valueOf(19.99244f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ni", Float.valueOf(57.93535f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("[18O]", Float.valueOf(17.99916f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("O", Float.valueOf(15.99491f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("P", Float.valueOf(30.97376f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Pb", Float.valueOf(207.97660f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Po", Float.valueOf(208.98240f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Pt", Float.valueOf(194.96480f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ra", Float.valueOf(226.02540f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Rb", Float.valueOf(84.91179f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Re", Float.valueOf(186.95580f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ru", Float.valueOf(101.90430f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("S", Float.valueOf(31.97207f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Sb", Float.valueOf(120.90380f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Se", Float.valueOf(79.91652f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Si", Float.valueOf(27.97693f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Sn", Float.valueOf(119.90220f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Tc", Float.valueOf(97.90722f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Te", Float.valueOf(129.90620f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Ti", Float.valueOf(47.94795f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Tl", Float.valueOf(204.97440f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("U", Float.valueOf(238.05080f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("V", Float.valueOf(50.94396f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("W", Float.valueOf(183.95090f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Y", Float.valueOf(88.90585f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Zn", Float.valueOf(63.92915f)); //$NON-NLS-1$
		MONOISOTOPIC_MASSES.put("Zr", Float.valueOf(89.90470f)); //$NON-NLS-1$
	}
}