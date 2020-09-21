package uk.ac.liverpool.metfraglib;

import javolution.util.FastBitSet;

/**
 * 
 * @author neilswainston
 */
public class FastBitArray {

	/**
	 * array to store bits
	 */
	private final FastBitSet bitArray;

	/**
	 * 
	 */
	private final short size;
	
	/**
	 * Initialises BitArray with specified number of bits.
	 * 
	 * @param size
	 */
	public FastBitArray(final int size) {
		this.bitArray = new FastBitSet();
		this.size = (short)size;
	}

	/**
	 * Initialises BitArray with specified number of bits and value at all positions.
	 * 
	 * @param size
	 * @param value
	 */
	public FastBitArray(final int size, final boolean value) {
		this.bitArray = new FastBitSet();
		if(value) this.bitArray.flip(0, size);
		this.size = (short)size;
	}
	
	/**
	 * Set bit at position idx to true.
	 * 
	 * Warning: function does not check for IndexOutOfBounds
	 * 
	 * @param idx
	 */
	public void set(final int idx) {
		this.bitArray.set(idx);
	}

	/**
	 * set bit at position n to given value val
	 * Warning: function does not check for IndexOutOfBounds
	 * 
	 * @param index
	 */
	public void set(int n, boolean val) {
		this.bitArray.set(n, val);
	}
	
	/**
	 * returns bit at index n
	 * returns false if n < 0 and n > size of BitArray
	 * 
	 * @param n
	 * @return
	 */
	public boolean get(int n) {
		return this.bitArray.get(n);
	}
	
	/**
	 * returns BitArray number of bits
	 * 
	 * @return
	 */
	public short getSize() {
		return this.size;
	}
	
	/**
	 * get integer array of indices that are set to true
	 * 
	 * @return int[]
	 */
	public int[] getSetIndices() {
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
	
	@Override
	public boolean equals(Object object) {
		FastBitArray bitArray = (FastBitArray)object;
		if(this.getSize() != bitArray.getSize()) return false;
		for(int i = 0; i < this.getSize(); i++) 
			if(this.bitArray.get(i) != bitArray.get(i)) return false;
		return true;
	}
	
	@Override
	public int hashCode() {
		java.util.BitSet bitset = new java.util.BitSet(this.size);
		for(int i = 0; i < this.getSize(); i++) 
			if(this.bitArray.get(i)) bitset.set(i);
		return bitset.hashCode();
	}
	
	/**
	 * returns true if the current and the given BitArray are of same length and have exactly the same
	 * bits set to true 
	 * 
	 * @param bitArray
	 * @return
	 */
	public boolean equals(String bitString) {
		if(this.getSize() != bitString.length()) return false;
		for(int i = 0; i < this.getSize(); i++) 
			if((this.bitArray.get(i) && bitString.charAt(i) == '0') || (!this.bitArray.get(i) && bitString.charAt(i) == '1')) return false;
		return true;
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
	public int getLastSetBit() {
		for(int i = this.getSize() - 1; i >= 0; i--)
			if(this.bitArray.get(i)) return i;
		return -1;
	}
	
	/**
	 * sets indeces of BitArray in the given integer array to true 
	 * 
	 * @param bitIndexes
	 */
	public void setBits(int[] bitIndexes) {
		for(int i = 0; i < bitIndexes.length; i++) {
			if(bitIndexes[i] < this.getSize() && bitIndexes[i] >= 0)
				this.bitArray.set(bitIndexes[i]);
			else
				System.err.println("Warning: Could not set bit at position " + bitIndexes[i] + " to true. Out of range!"); //$NON-NLS-1$ //$NON-NLS-2$
		}
	}
	
	/**
	 * 
	 */
	@Override
	public FastBitArray clone() {
		FastBitArray clone = new FastBitArray(this.getSize());
		for(int i = 0; i < this.getSize(); i++)
			if(this.bitArray.get(i)) clone.set(i);
		return clone;
	}
}