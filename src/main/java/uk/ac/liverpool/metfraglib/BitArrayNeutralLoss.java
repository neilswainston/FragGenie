package uk.ac.liverpool.metfraglib;

import uk.ac.liverpool.metfraglib.FastBitArray;

/**
 * 
 * @author neilswainston
 */
class BitArrayNeutralLoss {
	/**
	 * 
	 */
	private final FastBitArray[] neutralLossAtoms;

	/**
	 * 
	 * @param numberNeutralLosses
	 */
	BitArrayNeutralLoss(final int numberNeutralLosses) {
		this.neutralLossAtoms = new FastBitArray[numberNeutralLosses];
	}

	/**
	 * 
	 * @param index
	 * @return
	 */
	FastBitArray getNeutralLossAtomFastBitArray(final int index) {
		return this.neutralLossAtoms[index];
	}

	/**
	 * 
	 * @return int
	 */
	int getNumberNeutralLosses() {
		return this.neutralLossAtoms.length;
	}

	/**
	 * 
	 * @param index
	 * @param neutralLossAtoms
	 */
	void setNeutralLoss(final int index, final FastBitArray neutralLossAtoms) {
		this.neutralLossAtoms[index] = neutralLossAtoms;
	}
}