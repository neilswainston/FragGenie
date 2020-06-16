package uk.ac.liverpool.metfraglib.fragment;

import de.ipbhalle.metfraglib.FastBitArray;
import de.ipbhalle.metfraglib.interfaces.IMolecularStructure;
import de.ipbhalle.metfraglib.precursor.AbstractTopDownBitArrayPrecursor;

public abstract class AbstractTopDownBitArrayFragment extends DefaultBitArrayFragment {

	/*
	 * value needed during fragment generation stores whether the fragment was
	 * generated by a ring bond cleavage without losing an additional atom which
	 * won't result in a new valid fragment this variable can then be used to set
	 * the proper precursor fragment
	 */
	protected boolean wasRingCleavedFragment;
	protected byte addedToQueueCounts;
	protected short lastSkippedBond;
	protected boolean hasMatchedChild;

	public AbstractTopDownBitArrayFragment(AbstractTopDownBitArrayPrecursor precursor) {
		super(precursor);
		this.wasRingCleavedFragment = false;
		this.addedToQueueCounts = 0;
		this.lastSkippedBond = -1;
		this.hasMatchedChild = false;
	}

	public AbstractTopDownBitArrayFragment(FastBitArray atomsFastBitArray, FastBitArray bondsFastBitArray,
			FastBitArray brokenBondsFastBitArray) {
		super(atomsFastBitArray, bondsFastBitArray, brokenBondsFastBitArray);
		this.lastSkippedBond = -1;
	}

	public abstract AbstractTopDownBitArrayFragment[] traverseMolecule(IMolecularStructure precursorMolecule,
			short bondIndexToRemove, short[] indecesOfBondConnectedAtoms) throws Exception;

	public int getMaximalIndexOfRemovedBond() {
		return this.brokenBondsFastBitArray.getLastSetBit();
	}

	public boolean hasMatchedChild() {
		return this.hasMatchedChild;
	}

	public void setHasMatchedChild(boolean hasMatchedChild) {
		this.hasMatchedChild = hasMatchedChild;
	}

	public short getLastSkippedBond() {
		return this.lastSkippedBond;
	}

	public void setLastSkippedBond(short lastSkippedBond) {
		this.lastSkippedBond = lastSkippedBond;
	}

	public byte getAddedToQueueCounts() {
		return this.addedToQueueCounts;
	}

	public void setAddedToQueueCounts(byte addedToQueueCounts) {
		this.addedToQueueCounts = addedToQueueCounts;
	}

	public boolean isWasRingCleavedFragment() {
		return this.wasRingCleavedFragment;
	}

	public void setWasRingCleavedFragment(boolean wasRingCleavedFragment) {
		this.wasRingCleavedFragment = wasRingCleavedFragment;
	}
}