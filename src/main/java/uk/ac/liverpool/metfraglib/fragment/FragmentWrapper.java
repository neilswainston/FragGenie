package uk.ac.liverpool.metfraglib.fragment;

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
	private Fragment fragment;

	/**
	 * 
	 * @param wrappedFragment
	 */
	public FragmentWrapper(final Fragment fragment) {
		this.fragment = fragment;
	}

	/**
	 * 
	 * @param root
	 * @param currentPeakIndexPointer
	 */
	public FragmentWrapper(final Fragment fragment, final int peakIndex) {
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
	public Fragment getWrappedFragment() {
		return this.fragment;
	}

	/**
	 * 
	 * @param fragment
	 */
	public void setWrappedFragment(final Fragment fragment) {
		this.fragment = fragment;
	}
}