package uk.ac.liverpool.metfraglib.fragmenterassignerscorer;

import java.util.LinkedList;
import java.util.Queue;
import java.util.Set;
import java.util.TreeSet;

import uk.ac.liverpool.metfraglib.fragment.Fragment;
import uk.ac.liverpool.metfraglib.fragmenter.Fragmenter;
import uk.ac.liverpool.metfraglib.precursor.Precursor;

public class TopDownFragmenterAssignerScorer {

	/**
	 * 
	 */
	private final static double PROTON_MASS = 1.00728;
	
	/**
	 * 
	 */
	private final int maximumTreeDepth;

	/**
	 * 
	 */
	private final Precursor precursor;

	/**
	 * 
	 */
	private final Fragmenter fragmenter;

	/**
	 * 
	 * @param prec
	 * @param maxTreeDepth
	 * @throws Exception
	 */
	public TopDownFragmenterAssignerScorer(final Precursor prec, int maxTreeDepth) throws Exception {
		this.precursor = prec;
		this.maximumTreeDepth = maxTreeDepth;
		this.fragmenter = new Fragmenter(this.precursor);
	}

	/**
	 * 
	 * @return
	 * @throws Exception
	 */
	public double[] getMasses() throws Exception {
		final Set<Double> masses = new TreeSet<>();

		// Generate root fragment to start fragmentation:
		Queue<Fragment> fragments = new LinkedList<>();
		fragments.add(new Fragment(this.precursor));

		for(int k = 1; k <= this.maximumTreeDepth; k++) {
			Queue<Fragment> newFragments = new LinkedList<>();
			
			// Use each fragment that is marked as to be processed:
			while(!fragments.isEmpty()) {
				// Generate fragments of new tree depth:
				final Fragment fragment = fragments.poll();
				
				// Loop over all child fragments from precursor fragment:
				for(final Fragment childFragment : this.fragmenter.getFragmentsOfNextTreeDepth(fragment)) {
					masses.add(Double.valueOf(childFragment.getMonoisotopicMass(this.precursor) + PROTON_MASS));

					// Mark current fragment for further fragmentation:
					if(this.maximumTreeDepth > 0) {
						newFragments.add(childFragment);
					}	
				}
			}
			
			fragments = newFragments;
		}

		final double[] massesArray = new double[masses.size()];
		int i = 0;

		for (Double mass : masses) {
			massesArray[i++] = mass.doubleValue();
		}

		return massesArray;
	}
}