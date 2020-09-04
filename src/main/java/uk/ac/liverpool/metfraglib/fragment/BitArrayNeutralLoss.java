package uk.ac.liverpool.metfraglib.fragment;

import de.ipbhalle.metfraglib.FastBitArray;

/**
 * 
 * @author neilswainston
 */
public class BitArrayNeutralLoss {
	/**
	 * 
	 */
	private final FastBitArray[] neutralLossAtoms;

	/**
	 * 
	 * @param numberNeutralLosses
	 */
	public BitArrayNeutralLoss(final int numberNeutralLosses) {
		this.neutralLossAtoms = new FastBitArray[numberNeutralLosses];
	}

	/**
	 * 
	 * @param index
	 * @return
	 */
	public FastBitArray getNeutralLossAtomFastBitArray(final int index) {
		return this.neutralLossAtoms[index];
	}

	/**
	 * 
	 * @return int
	 */
	public int getNumberNeutralLosses() {
		return this.neutralLossAtoms.length;
	}

	/**
	 * 
	 * @param index
	 * @param neutralLossAtoms
	 */
	public void setNeutralLoss(final int index, final FastBitArray neutralLossAtoms) {
		this.neutralLossAtoms[index] = neutralLossAtoms;
	}
}