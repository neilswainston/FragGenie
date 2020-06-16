package uk.ac.liverpool.metfraglib.fragment;

import de.ipbhalle.metfraglib.interfaces.IFragment;

public class AbstractTopDownBitArrayFragmentWrapper {

	protected Integer currentPeakIndexPointer;
	protected IFragment wrappedFragment;

	public AbstractTopDownBitArrayFragmentWrapper(IFragment wrappedFragment) {
		this.wrappedFragment = wrappedFragment;
	}

	public AbstractTopDownBitArrayFragmentWrapper(IFragment root, Integer currentPeakIndexPointer) {
		this.wrappedFragment = root;
		this.currentPeakIndexPointer = currentPeakIndexPointer;
	}

	public Integer getCurrentPeakIndexPointer() {
		return this.currentPeakIndexPointer;
	}

	public void setCurrentPeakIndexPointer(Integer currentPeakIndexPointer) {
		this.currentPeakIndexPointer = currentPeakIndexPointer;
	}

	public IFragment getWrappedFragment() {
		return this.wrappedFragment;
	}

	public void setWrappedFragment(AbstractTopDownBitArrayFragment wrappedFragment) {
		this.wrappedFragment = wrappedFragment;
	}

	public void shallowNullify() {
		this.currentPeakIndexPointer = null;
	}
}
