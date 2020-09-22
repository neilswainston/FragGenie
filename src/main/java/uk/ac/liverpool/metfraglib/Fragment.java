package uk.ac.liverpool.metfraglib;

import java.util.Map;
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
	private short lastSkippedBond = -1;

	/**
	 * 
	 */
	private byte treeDepth = 0;

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
	 * @return short
	 */
	short getLastSkippedBond() {
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
	 * @param addedToQueueCounts
	 */
	void setAddedToQueueCounts(final byte addedToQueueCounts) {
		this.addedToQueueCounts = addedToQueueCounts;
	}

	/**
	 * 
	 * @param lastSkippedBond
	 */
	void setLastSkippedBond(final short lastSkippedBond) {
		this.lastSkippedBond = lastSkippedBond;
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
	Fragment[] traverseMolecule(final Precursor precursorMolecule, final int bondIndexToRemove,
			final int[] indecesOfBondConnectedAtoms) throws Exception {

		/*
		 * generate first fragment
		 */
		FastBitArray atomArrayOfNewFragment_1 = new FastBitArray(precursorMolecule.getNonHydrogenAtomCount(), false);
		FastBitArray bondArrayOfNewFragment_1 = new FastBitArray(precursorMolecule.getNonHydrogenBondCount(), false);
		FastBitArray brokenBondArrayOfNewFragment_1 = this.brokenBondsArray.clone();
		int[] numberHydrogensOfNewFragment = new int[1];

		/*
		 * traverse to first direction from atomIndex connected by broken bond
		 */
		boolean stillOneFragment = this.traverseSingleDirection(precursorMolecule, indecesOfBondConnectedAtoms[0],
				indecesOfBondConnectedAtoms[1], bondIndexToRemove, atomArrayOfNewFragment_1, bondArrayOfNewFragment_1,
				brokenBondArrayOfNewFragment_1, numberHydrogensOfNewFragment);

		Fragment firstNewGeneratedFragment = new Fragment(this.prec, atomArrayOfNewFragment_1, bondArrayOfNewFragment_1,
				brokenBondArrayOfNewFragment_1);

		/*
		 * only one fragment is generated when a ring bond was broken
		 */
		if (stillOneFragment) {
			firstNewGeneratedFragment.treeDepth = this.treeDepth;
			firstNewGeneratedFragment.setAddedToQueueCounts((byte) (this.getAddedToQueueCounts() + 1));
			Fragment[] newFrags = { firstNewGeneratedFragment };
			return newFrags;
		}
		/*
		 * generate second fragment
		 */
		FastBitArray atomArrayOfNewFragment_2 = new FastBitArray(precursorMolecule.getNonHydrogenAtomCount(), false);
		FastBitArray bondArrayOfNewFragment_2 = new FastBitArray(precursorMolecule.getNonHydrogenBondCount(), false);
		FastBitArray brokenBondArrayOfNewFragment_2 = this.brokenBondsArray.clone();
		numberHydrogensOfNewFragment[0] = 0;

		/*
		 * traverse the second direction from atomIndex connected by broken bond
		 */
		this.traverseSingleDirection(precursorMolecule, indecesOfBondConnectedAtoms[1], indecesOfBondConnectedAtoms[0],
				bondIndexToRemove, atomArrayOfNewFragment_2, bondArrayOfNewFragment_2, brokenBondArrayOfNewFragment_2,
				numberHydrogensOfNewFragment);

		Fragment secondNewGeneratedFragment = new Fragment(this.prec, atomArrayOfNewFragment_2,
				bondArrayOfNewFragment_2, brokenBondArrayOfNewFragment_2);

		firstNewGeneratedFragment.treeDepth = (byte) (this.treeDepth + 1);

		secondNewGeneratedFragment.treeDepth = (byte) (this.treeDepth + 1);

		Fragment[] newFrags = { firstNewGeneratedFragment, secondNewGeneratedFragment };

		return newFrags;

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
	private boolean traverseSingleDirection(Precursor precursorMolecule, int startAtomIndex, int endAtomIndex,
			int bondIndexToRemove, FastBitArray atomArrayOfNewFragment, FastBitArray bondArrayOfNewFragment,
			FastBitArray brokenBondArrayOfNewFragment, int[] numberHydrogensOfNewFragment) {
		final FastBitArray bondFastBitArrayOfCurrentFragment = this.bondsArray;
		/*
		 * when traversing the fragment graph then we want to know if we already visited
		 * a node (atom) need to be done for checking of ringed structures if traversed
		 * an already visited atom, then no new fragment was generated
		 */
		FastBitArray visited = new FastBitArray(precursorMolecule.getNonHydrogenAtomCount(), false);
		numberHydrogensOfNewFragment[0] = 0;

		/*
		 * traverse molecule in the first direction
		 */
		java.util.Stack<int[]> toProcessConnectedAtoms = new java.util.Stack<>();
		java.util.Stack<Integer> toProcessAtom = new java.util.Stack<>();
		toProcessConnectedAtoms.push(precursorMolecule.getConnectedAtomIndecesOfAtomIndex(startAtomIndex));
		toProcessAtom.push(Integer.valueOf(startAtomIndex));
		visited.set(startAtomIndex);
		boolean stillOneFragment = false;
		/*
		 * set the first atom of possible new fragment atom is of the one direction of
		 * cutted bond
		 */
		atomArrayOfNewFragment.set(startAtomIndex);
		numberHydrogensOfNewFragment[0] += precursorMolecule.getNumberHydrogensConnectedToAtomIndex(startAtomIndex);
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
				numberHydrogensOfNewFragment[0] += precursorMolecule
						.getNumberHydrogensConnectedToAtomIndex(nextAtoms[i]);
				bondArrayOfNewFragment
						.set(precursorMolecule.getBondIndexFromAtomAdjacencyList(midAtom, nextAtoms[i]) - 1);
				toProcessConnectedAtoms.push(precursorMolecule.getConnectedAtomIndecesOfAtomIndex(nextAtoms[i]));
				toProcessAtom.push(Integer.valueOf(nextAtoms[i]));
			}
		}

		brokenBondArrayOfNewFragment.set(bondIndexToRemove);
		bondArrayOfNewFragment.set(bondIndexToRemove, false);

		return stillOneFragment;
	}
}