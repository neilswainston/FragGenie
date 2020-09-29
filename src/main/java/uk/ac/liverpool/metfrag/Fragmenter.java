package uk.ac.liverpool.metfrag;

import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * 
 * @author neilswainston
 */
public class Fragmenter {
	
	/**
	 * 
	 */
	private final Fragment precursorFragment;
	
	/**
	 * 
	 * @param precursor
	 * @throws IOException
	 * @throws CDKException
	 */
	public Fragmenter(final IAtomContainer precursor) throws IOException, CDKException {
		this.precursorFragment = new Fragment(precursor);
	}

	/**
	 * 
	 * @return Collection<List<Object>>
	 */
	public Collection<Fragment> getFragments(final int maxBrokenBonds) {
		final Collection<Fragment> fragments = new TreeSet<>();
		fragment(this.precursorFragment, fragments, maxBrokenBonds);
		return fragments;
	}
	
	/**
	 * 
	 * @param fragment
	 * @param fragments
	 * @param maxBrokenBonds
	 */
	private void fragment(final Fragment fragment, final Collection<Fragment> fragments, final int maxBrokenBonds) {
		if(fragment.getNumBrokenBonds() <= maxBrokenBonds) {
			fragments.add(fragment);
			
			final boolean[] bondsArray = fragment.getBondsArray();
			final boolean[] brokenBondsArray = fragment.getBrokenBondsArray();
			
			if(fragment.getNumBrokenBonds() < maxBrokenBonds) {
				for(int bondIdx = 0; bondIdx < bondsArray.length; bondIdx++) {
					if(bondsArray[bondIdx] && !brokenBondsArray[bondIdx]) {
						for(final Fragment childFragment : fragment.fragment(bondIdx)) {
							fragment(childFragment, fragments, maxBrokenBonds);
						}
					}
				}
			}
		}
	}
}