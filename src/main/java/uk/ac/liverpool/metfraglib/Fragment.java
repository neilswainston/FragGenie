package uk.ac.liverpool.metfraglib;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.TreeMap;

/**
 * 
 * @author neilswainston
 */
public class Fragment {

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
		this(precursor, new FastBitArray(precursor.getNonHydrogenAtomCount(), true),
				new FastBitArray(precursor.getNonHydrogenBondCount(), true),
				new FastBitArray(precursor.getNonHydrogenBondCount(), false));
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
				final String element = this.prec.getAtom(i);

				if (elementCount.get(element) == null) {
					elementCount.put(element, Integer.valueOf(1));
				} else {
					elementCount.put(element, Integer.valueOf(elementCount.get(element).intValue() + 1));
				}

				final int hCount = this.prec.getNumberHydrogensConnectedToAtomIndex(i);

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
	 * @return Collection<List<Object>>
	 */
	public Collection<List<Object>> getBrokenBonds() {
		final Collection<List<Object>> bonds = new ArrayList<>();
		
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
		final FastBitArray newAtomArray = new FastBitArray(precursorMolecule.getNonHydrogenAtomCount(), false);
		final FastBitArray newBondArray = new FastBitArray(precursorMolecule.getNonHydrogenBondCount(), false);
		final FastBitArray currentBondArray = this.bondsArray;
		
		// When traversing the fragment graph, we want to know if we already visited
		// a node (atom) to check for ringed structures.
		// If traversed an already visited atom, then no new fragment was generated.
		final FastBitArray visited = new FastBitArray(precursorMolecule.getNonHydrogenAtomCount(), false);

		// Traverse molecule in the first direction:
		final Stack<int[]> toProcessConnectedAtoms = new Stack<>();
		final Stack<Integer> toProcessAtom = new Stack<>();
		
		toProcessConnectedAtoms.push(precursorMolecule.getConnectedAtomIndecesOfAtomIndex(startAtom));
		toProcessAtom.push(Integer.valueOf(startAtom));
		visited.set(startAtom);
		
		boolean singleFragment = false;
		
		// Set the first atom of possible new fragment atom as the one direction of cut bond.
		newAtomArray.set(startAtom);
		
		while (!toProcessConnectedAtoms.isEmpty()) {
			final int midAtom = toProcessAtom.pop().intValue();
			
			for (int nextAtom : toProcessConnectedAtoms.pop()) {
				// Did we visit the current atom already?
				final int currentBond = precursorMolecule.getBondIndexFromAtomAdjacencyList(nextAtom, midAtom) - 1;

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

				newBondArray.set(precursorMolecule.getBondIndexFromAtomAdjacencyList(midAtom, nextAtom) - 1);
				toProcessConnectedAtoms.push(precursorMolecule.getConnectedAtomIndecesOfAtomIndex(nextAtom));
				toProcessAtom.push(Integer.valueOf(nextAtom));
			}
		}

		brokenBondArrayOfNewFragment.set(removeBond);
		newBondArray.set(removeBond, false);
		
		final Fragment newFragment = new Fragment(precursorMolecule, newAtomArray, newBondArray, brokenBondArrayOfNewFragment);

		return new Object[] {Boolean.valueOf(singleFragment), newFragment};
	}
}