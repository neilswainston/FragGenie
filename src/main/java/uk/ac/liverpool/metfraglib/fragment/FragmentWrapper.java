package uk.ac.liverpool.metfraglib.fragment;

import de.ipbhalle.metfraglib.interfaces.IFragment;

/**
 * 
 * @author neilswainston
 */
public class FragmentWrapper {

	/**
	 * 
	 */
	private int peakIndex;
	
	/**
	 * 
	 */
	private IFragment fragment;

	/**
	 * 
	 * @param wrappedFragment
	 */
	public FragmentWrapper(final IFragment fragment) {
		this.fragment = fragment;
	}

	/**
	 * 
	 * @param root
	 * @param currentPeakIndexPointer
	 */
	public FragmentWrapper(final IFragment fragment, final int peakIndex) {
		this.fragment = fragment;
		this.peakIndex = peakIndex;
	}

	/**
	 * 
	 * @return
	 */
	public int getCurrentPeakIndexPointer() {
		return this.peakIndex;
	}

	/**
	 * 
	 * @param peakIndex
	 */
	public void setCurrentPeakIndexPointer(final int peakIndex) {
		this.peakIndex = peakIndex;
	}

	/**
	 * 
	 * @return
	 */
	public IFragment getWrappedFragment() {
		return this.fragment;
	}

	/**
	 * 
	 * @param fragment
	 */
	public void setWrappedFragment(final AbstractTopDownBitArrayFragment fragment) {
		this.fragment = fragment;
	}
}