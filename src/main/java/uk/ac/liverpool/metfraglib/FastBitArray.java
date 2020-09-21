package uk.ac.liverpool.metfraglib;

import javolution.util.FastBitSet;

/**
 * 
 * @author neilswainston
 */
class FastBitArray {

	/**
	 * array to store bits
	 */
	private final FastBitSet bitArray;

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
		this.bitArray = new FastBitSet();
		this.size = size;
		
		if(set) {
			this.bitArray.flip(0, size);
		}
	}
	
	/**
	 * Set bit at position idx to true.
	 * 
	 * Warning: function does not check for IndexOutOfBounds
	 * 
	 * @param idx
	 */
	void set(final int idx) {
		this.bitArray.set(idx);
	}

	/**
	 * set bit at position n to given value val
	 * Warning: function does not check for IndexOutOfBounds
	 * 
	 * @param index
	 */
	void set(int n, boolean val) {
		this.bitArray.set(n, val);
	}
	
	/**
	 * returns bit at index n
	 * returns false if n < 0 and n > size of BitArray
	 * 
	 * @param n
	 * @return
	 */
	boolean get(int n) {
		return this.bitArray.get(n);
	}
	
	/**
	 * returns BitArray number of bits
	 * 
	 * @return int
	 */
	int getSize() {
		return this.size;
	}
	
	/**
	 * get integer array of indices that are set to true
	 * 
	 * @return int[]
	 */
	int[] getSetIndices() {
		int[] setIndeces = new int[this.bitArray.cardinality()];
		int index = 0;
		for(int i = 0; i < this.getSize(); i++) {
			if(this.bitArray.get(i)) {
				setIndeces[index] = i;
				index++;
			}
		}
		return setIndeces;
	}
	
	/**
	 * returns string with with true positions '1' and false positions '0' 
	 */
	@Override
	public String toString() {
		char[] set = new char[this.getSize()];
		for(int i = 0; i < this.getSize(); i++) {
			set[i] = this.bitArray.get(i) ? '1' : '0';
		}
		return String.valueOf(set);
	}
	
	/**
	 * returns last position of BitArray set to true
	 * returns -1 if there is no bit set to true
	 * 
	 */
	int getLastSetBit() {
		for(int i = this.getSize() - 1; i >= 0; i--)
			if(this.bitArray.get(i)) return i;
		return -1;
	}
	
	/**
	 * sets indeces of BitArray in the given integer array to true 
	 * 
	 * @param bitIndexes
	 */
	void setBits(int[] bitIndexes) {
		for(int i = 0; i < bitIndexes.length; i++) {
			this.bitArray.set(bitIndexes[i]);
		}
	}
	
	@Override
	public FastBitArray clone() {
		FastBitArray clone = new FastBitArray(this.getSize(), false);
		for(int i = 0; i < this.getSize(); i++)
			if(this.bitArray.get(i)) clone.set(i);
		return clone;
	}
}