package uk.ac.liverpool.metfraglib.fragment;

/**
 * 
 * @author neilswainston
 */
public class FragmentWrapper {

	/**
	 * 
	 */
	private final Fragment fragment;

	/**
	 * 
	 */
	private final int peakIndex;

	/**
	 * 
	 * @param fragment
	 * @param peakIndex
	 */
	public FragmentWrapper(final Fragment fragment, final int peakIndex) {
		this.fragment = fragment;
		this.peakIndex = peakIndex;
	}

	/**
	 * 
	 * @return
	 */
	public Fragment getFragment() {
		return this.fragment;
	}

	/**
	 * 
	 * @return
	 */
	public int getPeakIndex() {
		return this.peakIndex;
	}
}