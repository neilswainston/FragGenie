package uk.ac.liverpool.metfraglib;

import javolution.util.FastBitSet;

/**
 * 
 * @author neilswainston
 */
class FastBitArray extends FastBitSet {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * 
	 */
	private final int size;

	/**
	 * Initialises BitArray with specified number of bits.
	 * 
	 * @param size
	 * @param set
	 */
	FastBitArray(final int size, final boolean set) {
		this.size = size;
		
		if(set) {
			this.flip(0, size);
		}
	}
	
	/**
	 * Returns number of bits.
	 * 
	 * @return int
	 */
	int getSize() {
		return this.size;
	}
	
	/**
	 * Get integer array of indices that are set to true.
	 * 
	 * @return int[]
	 */
	int[] getSetIndices() {
		final int[] setIndices = new int[this.cardinality()];
		
		int index = 0;
		
		for(int i = 0; i < this.size; i++) {
			if(this.get(i)) {
				setIndices[index++] = i;
			}
		}
		return setIndices;
	}
	
	/**
	 * returns last position of BitArray set to true
	 * returns -1 if there is no bit set to true
	 */
	int getLastSetBit() {
		for(int i = this.getSize() - 1; i >= 0; i--) {
			if(this.get(i)) {
				return i;
			}
		}
		
		return -1;
	}
	
	/**
	 * sets indeces of BitArray in the given integer array to true 
	 * 
	 * @param bitIndexes
	 */
	void setBits(final int[] bitIndexes) {
		for(int i = 0; i < bitIndexes.length; i++) {
			this.set(bitIndexes[i]);
		}
	}
	
	@Override
	public String toString() {
		final char[] set = new char[this.size];
		
		for(int i = 0; i < this.size; i++) {
			set[i] = this.get(i) ? '1' : '0';
		}
		
		return String.valueOf(set);
	}
	
	@Override
	public FastBitArray clone() {
		final FastBitArray clone = new FastBitArray(this.size, false);
		
		for(int i = 0; i < this.getSize(); i++) {
			if(this.get(i)) {
				clone.set(i);
			}
		}
		
		return clone;
	}
}