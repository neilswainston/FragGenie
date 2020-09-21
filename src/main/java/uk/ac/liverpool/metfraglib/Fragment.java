package uk.ac.liverpool.metfraglib;


import java.util.Map;
import java.util.TreeMap;

import de.ipbhalle.metfraglib.FastBitArray;
import uk.ac.liverpool.metfraglib.Precursor;

public class Fragment {

	private int addedToQueueCounts;

	private FastBitArray atomsFastBitArray;
	private FastBitArray bondsFastBitArray;
	private FastBitArray brokenBondsFastBitArray;
	private short lastSkippedBond = -1;
	/**
	 * 
	 */
	private final Precursor prec;
	private byte treeDepth = 0;

	/**
	 * constructor setting all bits of atomsFastBitArray and bondsFastBitArray to
	 * true entire structure is represented
	 * 
	 * @param precursor
	 * @throws AtomTypeNotKnownFromInputListException
	 */
	public Fragment(final Precursor precursor) {
		this.prec = precursor;
		this.addedToQueueCounts = 0;
		this.atomsFastBitArray = new FastBitArray(precursor.getNonHydrogenAtomCount(), true);
		this.bondsFastBitArray = new FastBitArray(precursor.getNonHydrogenBondCount(), true);
		this.brokenBondsFastBitArray = new FastBitArray(precursor.getNonHydrogenBondCount());
	}

	/**
	 * constructor setting bits of atomsFastBitArray and bondsFastBitArray by given
	 * ones
	 * 
	 * @param precursor
	 * @param atomsFastBitArray
	 * @param bondsFastBitArray
	 * @param brokenBondsFastBitArray
	 * @param numberHydrogens
	 * @throws Exception
	 */
	private Fragment(final Precursor precursor, final FastBitArray atomsFastBitArray,
			final FastBitArray bondsFastBitArray, final FastBitArray brokenBondsFastBitArray) throws Exception {
		this.prec = precursor;
		this.atomsFastBitArray = atomsFastBitArray;
		this.bondsFastBitArray = bondsFastBitArray;
		this.brokenBondsFastBitArray = brokenBondsFastBitArray;
	}

	/**
	 * 
	 * @return int
	 */
	public int getAddedToQueueCounts() {
		return this.addedToQueueCounts;
	}

	/**
	 * 
	 * @return FastBitArray
	 */
	public FastBitArray getAtomsFastBitArray() {
		return this.atomsFastBitArray;
	}

	/**
	 * 
	 * @return FastBitArray
	 */
	public FastBitArray getBondsFastBitArray() {
		return this.bondsFastBitArray;
	}

	/**
	 * 
	 * @return int[]
	 */
	public int[] getBrokenBondIndeces() {
		return this.brokenBondsFastBitArray.getSetIndeces();
	}

	/**
	 * 
	 * @return FastBitArray
	 */
	public FastBitArray getBrokenBondsFastBitArray() {
		return this.brokenBondsFastBitArray;
	}

	/**
	 * 
	 * @return short
	 */
	public short getLastSkippedBond() {
		return this.lastSkippedBond;
	}

	/**
	 * 
	 * @return float
	 */
	public float getMonoisotopicMass() {
		float mass = 0.0f;

		for (int i = 0; i < this.atomsFastBitArray.getSize(); i++) {
			if (this.atomsFastBitArray.get(i)) {
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
		final Map<String,Integer> elementCount = new TreeMap<>();

		for (int i = 0; i < this.atomsFastBitArray.getSize(); i++) {
			if (this.atomsFastBitArray.get(i)) {
				final String element = this.prec.getAtom(i);
				
				if(elementCount.get(element) == null) {
					elementCount.put(element, Integer.valueOf(1));
				}
				else {
					elementCount.put(element, Integer.valueOf(elementCount.get(element).intValue() + 1));
				}
				
				final int hCount = this.prec.getNumberHydrogensConnectedToAtomIndex(i);
				
				if(elementCount.get("H") == null) { //$NON-NLS-1$
					elementCount.put("H", Integer.valueOf(hCount)); //$NON-NLS-1$
				}
				else {
					elementCount.put("H", Integer.valueOf(elementCount.get("H").intValue() + hCount)); //$NON-NLS-1$ //$NON-NLS-2$
				}
			}
		}
		
		final StringBuilder builder = new StringBuilder();
		
		for(Map.Entry<String,Integer> entry : elementCount.entrySet()) {
			builder.append(entry.getKey());
			final int count = entry.getValue().intValue();
			
			if(count > 1) {
				builder.append(count);
			}
		}
		
		return builder.toString();
	}

	/**
	 * 
	 * @return int
	 */
	public int getMaximalIndexOfRemovedBond() {
		return this.brokenBondsFastBitArray.getLastSetBit();
	}

	/**
	 * 
	 * @param addedToQueueCounts
	 */
	public void setAddedToQueueCounts(final byte addedToQueueCounts) {
		this.addedToQueueCounts = addedToQueueCounts;
	}

	/**
	 * 
	 * @param lastSkippedBond
	 */
	public void setLastSkippedBond(final short lastSkippedBond) {
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
	public Fragment[] traverseMolecule(final Precursor precursorMolecule, final short bondIndexToRemove,
			final short[] indecesOfBondConnectedAtoms) throws Exception {

		/*
		 * generate first fragment
		 */
		FastBitArray atomArrayOfNewFragment_1 = new FastBitArray(precursorMolecule.getNonHydrogenAtomCount());
		FastBitArray bondArrayOfNewFragment_1 = new FastBitArray(precursorMolecule.getNonHydrogenBondCount());
		FastBitArray brokenBondArrayOfNewFragment_1 = this.getBrokenBondsFastBitArray().clone();
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
		FastBitArray atomArrayOfNewFragment_2 = new FastBitArray(precursorMolecule.getNonHydrogenAtomCount());
		FastBitArray bondArrayOfNewFragment_2 = new FastBitArray(precursorMolecule.getNonHydrogenBondCount());
		FastBitArray brokenBondArrayOfNewFragment_2 = this.getBrokenBondsFastBitArray().clone();
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
	private boolean traverseSingleDirection(Precursor precursorMolecule, short startAtomIndex, int endAtomIndex,
			int bondIndexToRemove, FastBitArray atomArrayOfNewFragment, FastBitArray bondArrayOfNewFragment,
			FastBitArray brokenBondArrayOfNewFragment, int[] numberHydrogensOfNewFragment) {
		FastBitArray bondFastBitArrayOfCurrentFragment = this.getBondsFastBitArray();
		/*
		 * when traversing the fragment graph then we want to know if we already visited
		 * a node (atom) need to be done for checking of ringed structures if traversed
		 * an already visited atom, then no new fragment was generated
		 */
		FastBitArray visited = new FastBitArray(precursorMolecule.getNonHydrogenAtomCount());
		numberHydrogensOfNewFragment[0] = 0;

		/*
		 * traverse molecule in the first direction
		 */
		java.util.Stack<short[]> toProcessConnectedAtoms = new java.util.Stack<>();
		java.util.Stack<Short> toProcessAtom = new java.util.Stack<>();
		toProcessConnectedAtoms.push(precursorMolecule.getConnectedAtomIndecesOfAtomIndex(startAtomIndex));
		toProcessAtom.push(Short.valueOf(startAtomIndex));
		visited.set(startAtomIndex);
		boolean stillOneFragment = false;
		/*
		 * set the first atom of possible new fragment atom is of the one direction of
		 * cutted bond
		 */
		atomArrayOfNewFragment.set(startAtomIndex);
		numberHydrogensOfNewFragment[0] += precursorMolecule.getNumberHydrogensConnectedToAtomIndex(startAtomIndex);
		while (!toProcessConnectedAtoms.isEmpty()) {
			short[] nextAtoms = toProcessConnectedAtoms.pop();
			short midAtom = toProcessAtom.pop().shortValue();
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
				toProcessAtom.push(Short.valueOf(nextAtoms[i]));
			}
		}

		brokenBondArrayOfNewFragment.set(bondIndexToRemove);
		bondArrayOfNewFragment.set(bondIndexToRemove, false);

		return stillOneFragment;
	}
}