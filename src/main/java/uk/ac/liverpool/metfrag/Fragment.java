package uk.ac.liverpool.metfrag;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.Stack;
import java.util.TreeMap;

/**
 * 
 * @author neilswainston
 */
public class Fragment implements Comparable<Fragment> {

	/**
	 * 
	 */
	private final Precursor prec;

	/**
	 * 
	 */
	private FastBitArray atomsArray;

	/**
	 * 
	 */
	private FastBitArray bondsArray;

	/**
	 * 
	 */
	private FastBitArray brokenBondsArray;

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
	Fragment(final Precursor precursor) {
		this(precursor, new FastBitArray(precursor.getAtomCount(), true),
				new FastBitArray(precursor.getBondCount(), true),
				new FastBitArray(precursor.getBondCount(), false));
	}

	/**
	 * Constructor.
	 * 
	 * @param precursor
	 * @param atoms
	 * @param bonds
	 * @param brokenBonds
	 */
	private Fragment(final Precursor precursor, final FastBitArray atoms, final FastBitArray bonds,
			final FastBitArray brokenBonds) {
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

		for (int i = 0; i < this.atomsArray.getSize(); i++) {
			if (this.atomsArray.get(i)) {
				mass += this.prec.getMassOfAtom(i);
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

		for (int i = 0; i < this.atomsArray.getSize(); i++) {
			if (this.atomsArray.get(i)) {
				final String element = this.prec.getAtomSymbol(i);

				if (elementCount.get(element) == null) {
					elementCount.put(element, Integer.valueOf(1));
				} else {
					elementCount.put(element, Integer.valueOf(elementCount.get(element).intValue() + 1));
				}

				final int hCount = this.prec.getImplicitHydrogenCount(i);

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
		
		for (int i = 0; i < this.brokenBondsArray.getSize(); i++) {
			if(this.brokenBondsArray.get(i)) {
				bonds.add(this.prec.getBond(i));
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
	 * @return FastBitArray
	 */
	FastBitArray getAtomsArray() {
		return this.atomsArray;
	}

	/**
	 * 
	 * @return FastBitArray
	 */
	FastBitArray getBondsArray() {
		return this.bondsArray;
	}

	/**
	 * 
	 * @return FastBitArray
	 */
	FastBitArray getBrokenBondsArray() {
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
	Fragment[] fragment(final Precursor precursorMolecule, final int removeBond, final int[] bondConnectedAtoms) throws Exception {

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
	 * Traverse the fragment to one direction starting from startAtom.
	 * 
	 * @param precursorMolecule
	 * @param startAtom
	 * @param endAtom
	 * @param removeBond
	 * @param brokenBondArrayOfNewFragment
	 * @return Object[]
	 */
	private Object[] traverse(final Precursor precursorMolecule, final int startAtom, final int endAtom,
			final int removeBond, final FastBitArray brokenBondArrayOfNewFragment) {
		final FastBitArray newAtomArray = new FastBitArray(precursorMolecule.getAtomCount(), false);
		final FastBitArray newBondArray = new FastBitArray(precursorMolecule.getBondCount(), false);
		final FastBitArray currentBondArray = this.bondsArray;
		
		// When traversing the fragment graph, we want to know if we already visited
		// a node (atom) to check for ringed structures.
		// If traversed an already visited atom, then no new fragment was generated.
		final FastBitArray visited = new FastBitArray(precursorMolecule.getAtomCount(), false);

		// Traverse molecule in the first direction:
		final Stack<int[]> toProcessConnectedAtoms = new Stack<>();
		final Stack<Integer> toProcessAtom = new Stack<>();
		
		toProcessConnectedAtoms.push(precursorMolecule.getConnectedAtomIndicesFromAtomIdx(startAtom));
		toProcessAtom.push(Integer.valueOf(startAtom));
		visited.set(startAtom);
		
		boolean singleFragment = false;
		
		// Set the first atom of possible new fragment atom as the one direction of cut bond.
		newAtomArray.set(startAtom);
		
		while (!toProcessConnectedAtoms.isEmpty()) {
			final int midAtom = toProcessAtom.pop().intValue();
			
			for (int nextAtom : toProcessConnectedAtoms.pop()) {
				// Did we visit the current atom already?
				final int currentBond = precursorMolecule.getBondIndex(nextAtom, midAtom) - 1;

				if (!currentBondArray.get(currentBond) || currentBond == removeBond) {
					continue;
				}
				
				// If we visited the current atom already, then we do not have to check it again:
				if (visited.get(nextAtom)) {
					newBondArray.set(currentBond);
					continue;
				}
				
				// If we reach the second atom of the cleaved bond then still one fragment is present:
				if (nextAtom == endAtom) {
					singleFragment = true;
				}

				visited.set(nextAtom);
				newAtomArray.set(nextAtom);

				newBondArray.set(precursorMolecule.getBondIndex(midAtom, nextAtom) - 1);
				toProcessConnectedAtoms.push(precursorMolecule.getConnectedAtomIndicesFromAtomIdx(nextAtom));
				toProcessAtom.push(Integer.valueOf(nextAtom));
			}
		}

		brokenBondArrayOfNewFragment.set(removeBond);
		newBondArray.set(removeBond, false);
		
		final Fragment newFragment = new Fragment(precursorMolecule, newAtomArray, newBondArray, brokenBondArrayOfNewFragment);

		return new Object[] {Boolean.valueOf(singleFragment), newFragment};
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
}