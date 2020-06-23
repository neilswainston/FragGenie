package uk.ac.liverpool.metfraglib.fragment;

import de.ipbhalle.metfraglib.FastBitArray;
import de.ipbhalle.metfraglib.additionals.NeutralLosses;

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
	 */
	private final byte neutralLossIndex;

	/**
	 * 
	 * @param numberNeutralLosses
	 * @param neutralLossIndex
	 */
	public BitArrayNeutralLoss(final int numberNeutralLosses, final byte neutralLossIndex) {
		this.neutralLossAtoms = new FastBitArray[numberNeutralLosses];
		this.neutralLossIndex = neutralLossIndex;
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

	/**
	 * 
	 * @return byte
	 */
	public byte getHydrogenDifference() {
		return new NeutralLosses().getHydrogenDifference(this.neutralLossIndex);
	}
}