package uk.ac.liverpool.metfraglib.fragment;

import de.ipbhalle.metfraglib.FastBitArray;
import de.ipbhalle.metfraglib.additionals.NeutralLosses;

public class BitArrayNeutralLoss {

	private FastBitArray[] neutralLossAtoms;
	/*
	 * relates to de.ipbhalle.metfrag.additionals.NeutralLosses
	 */
	private byte neutralLossIndex;
	
	public BitArrayNeutralLoss(int numberNeutralLosses, byte neutralLossIndex) {
		this.neutralLossAtoms = new FastBitArray[numberNeutralLosses];
		this.neutralLossIndex = neutralLossIndex;
	}
	
	public FastBitArray getNeutralLossAtomFastBitArray(int index) {
		return this.neutralLossAtoms[index];
	}
	
	public int getNumberNeutralLosses() {
		return this.neutralLossAtoms.length;
	}
	
	public void setNeutralLoss(int index, FastBitArray neutralLossAtoms) {
		this.neutralLossAtoms[index] = neutralLossAtoms;
	}
	
	public byte getHydrogenDifference() {
		return new NeutralLosses().getHydrogenDifference(this.neutralLossIndex);
	}
	
}
