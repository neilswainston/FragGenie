package uk.ac.liverpool.metfraglib;

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
	 * @param atomsArray
	 * @param bondsArray
	 * @param brokenBondsArray
	 */
	private Fragment(final Precursor precursor, final FastBitArray atomsArray, final FastBitArray bondsArray,
			final FastBitArray brokenBondsArray) {
		this.prec = precursor;
		this.atomsArray = atomsArray;
		this.bondsArray = bondsArray;
		this.brokenBondsArray = brokenBondsArray;
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
	 * @return float
	 */
	float getMonoisotopicMass() {
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
	String getFormula() {
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

				if (elementCount.get("H") == null) { //$NON-NLS-1$
					elementCount.put("H", Integer.valueOf(hCount)); //$NON-NLS-1$
				} else {
					elementCount.put("H", Integer.valueOf(elementCount.get("H").intValue() + hCount)); //$NON-NLS-1$ //$NON-NLS-2$
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
	 * main function of fragment generation traverse the given fragment and return
	 * two new fragments by removing bond with bondIndexToRemove
	 * 
	 * @param fragment
	 * @param bondNumber
	 * @param bondAtoms
	 * @return Fragment[]
	 * @throws Exception
	 */
	Fragment[] traverseMolecule(final Precursor precursorMolecule, final int bondIndexToRemove, final int[] bondConnectedAtoms) throws Exception {

		// Generate first fragment:
		final int[] numberHydrogensOfNewFragment = new int[1];

		/*
		 * traverse to first direction from atomIndex connected by broken bond
		 */
		final Object[] result1 = this.traverseSingleDirection(precursorMolecule, bondConnectedAtoms[0], bondConnectedAtoms[1],
				bondIndexToRemove, this.brokenBondsArray.clone(), numberHydrogensOfNewFragment);

		final Fragment fragment1 = (Fragment)result1[1];

		// Only one fragment is generated when a ring bond was broken:
		if (((Boolean)result1[0]).booleanValue()) {
			fragment1.treeDepth = this.treeDepth;
			fragment1.addedToQueueCounts = this.getAddedToQueueCounts() + 1;
			return new Fragment[] { fragment1 };
		}
		
		fragment1.treeDepth = this.treeDepth + 1;
		
		// Generate second fragment:
		numberHydrogensOfNewFragment[0] = 0;

		// Traverse the second direction from atomIndex connected by broken bond:
		final Object[] result2 = this.traverseSingleDirection(precursorMolecule, bondConnectedAtoms[1], bondConnectedAtoms[0],
				bondIndexToRemove, this.brokenBondsArray.clone(),
				numberHydrogensOfNewFragment);

		final Fragment fragment2 = (Fragment)result2[1];
		fragment2.treeDepth = this.treeDepth + 1;

		return new Fragment[] { fragment1, fragment2 };

	}

	/**
	 * traverse the fragment to one direction starting from startAtomIndex to set
	 * FastBitArrays of new fragment
	 * 
	 * @param startAtomIndex
	 * @param fragment
	 * @param bondIndexToRemove
	 * @param atomArrayOfNewFragment
	 * @param bondArrayOfNewFragment
	 * @param brokenBondArrayOfNewFragment
	 * @param numberHydrogensOfNewFragment
	 * @return
	 */
	private Object[] traverseSingleDirection(Precursor precursorMolecule, int startAtomIndex, int endAtomIndex,
			int bondIndexToRemove, FastBitArray brokenBondArrayOfNewFragment, int[] numberHydrogens) {
		final FastBitArray atomArrayOfNewFragment = new FastBitArray(precursorMolecule.getNonHydrogenAtomCount(), false);
		final FastBitArray bondArrayOfNewFragment = new FastBitArray(precursorMolecule.getNonHydrogenBondCount(), false);
		final FastBitArray bondFastBitArrayOfCurrentFragment = this.bondsArray;
		/*
		 * when traversing the fragment graph then we want to know if we already visited
		 * a node (atom) need to be done for checking of ringed structures if traversed
		 * an already visited atom, then no new fragment was generated
		 */
		FastBitArray visited = new FastBitArray(precursorMolecule.getNonHydrogenAtomCount(), false);
		numberHydrogens[0] = 0;

		/*
		 * traverse molecule in the first direction
		 */
		final Stack<int[]> toProcessConnectedAtoms = new Stack<>();
		final Stack<Integer> toProcessAtom = new Stack<>();
		
		toProcessConnectedAtoms.push(precursorMolecule.getConnectedAtomIndecesOfAtomIndex(startAtomIndex));
		toProcessAtom.push(Integer.valueOf(startAtomIndex));
		visited.set(startAtomIndex);
		boolean stillOneFragment = false;
		/*
		 * set the first atom of possible new fragment atom is of the one direction of
		 * cutted bond
		 */
		atomArrayOfNewFragment.set(startAtomIndex);
		numberHydrogens[0] += precursorMolecule.getNumberHydrogensConnectedToAtomIndex(startAtomIndex);
		while (!toProcessConnectedAtoms.isEmpty()) {
			int[] nextAtoms = toProcessConnectedAtoms.pop();
			int midAtom = toProcessAtom.pop().shortValue();
			for (int i = 0; i < nextAtoms.length; i++) {
				/*
				 * did we visit the current atom already?
				 */
				int currentBondNumber = precursorMolecule.getBondIndexFromAtomAdjacencyList(nextAtoms[i], midAtom) - 1;

				if (!bondFastBitArrayOfCurrentFragment.get(currentBondNumber)
						|| currentBondNumber == bondIndexToRemove) {
					continue;
				}
				/*
				 * if we visited the current atom already then we do not have to check it again
				 */
				if (visited.get(nextAtoms[i])) {
					bondArrayOfNewFragment.set(currentBondNumber);
					continue;
				}
				/*
				 * if we reach the second atom of the cleaved bond then still one fragment is
				 * present
				 */
				if (nextAtoms[i] == endAtomIndex) {
					stillOneFragment = true;
				}

				visited.set(nextAtoms[i]);
				atomArrayOfNewFragment.set(nextAtoms[i]);
				/*
				 * add number of hydrogens of current atom
				 */
				numberHydrogens[0] += precursorMolecule
						.getNumberHydrogensConnectedToAtomIndex(nextAtoms[i]);
				bondArrayOfNewFragment
						.set(precursorMolecule.getBondIndexFromAtomAdjacencyList(midAtom, nextAtoms[i]) - 1);
				toProcessConnectedAtoms.push(precursorMolecule.getConnectedAtomIndecesOfAtomIndex(nextAtoms[i]));
				toProcessAtom.push(Integer.valueOf(nextAtoms[i]));
			}
		}

		brokenBondArrayOfNewFragment.set(bondIndexToRemove);
		bondArrayOfNewFragment.set(bondIndexToRemove, false);
		
		final Fragment newFragment = new Fragment(precursorMolecule, atomArrayOfNewFragment, bondArrayOfNewFragment, brokenBondArrayOfNewFragment);

		return new Object[] {Boolean.valueOf(stillOneFragment), newFragment};
	}
}