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
	private final int maxTreeDepth;

	/**
	 * 
	 */
	private final Precursor prec;

	/**
	 * 
	 */
	private final Fragmenter fragmenter;

	/**
	 * 
	 * @param precursor
	 * @param maximumTreeDepth
	 * @throws Exception
	 */
	public TopDownFragmenterAssignerScorer(final Precursor precursor, int maximumTreeDepth) throws Exception {
		this.prec = precursor;
		this.maxTreeDepth = maximumTreeDepth;
		this.fragmenter = new Fragmenter(this.prec);
	}

	/**
	 * 
	 * @return
	 * @throws Exception
	 */
	public double[] getMasses() throws Exception {
		final Set<Double> masses = new TreeSet<>();

		Queue<Fragment> fragments = new LinkedList<>();
		fragments.add(new Fragment(this.prec));

		for (int k = 1; k <= this.maxTreeDepth; k++) {
			Queue<Fragment> newFragments = new LinkedList<>();

			while(!fragments.isEmpty()) {
				final Fragment fragment = fragments.poll();

				for(final Fragment childFragment : this.fragmenter.getFragmentsOfNextTreeDepth(fragment)) {
					masses.addAll(childFragment.getMasses(this.prec, 2, true));

					if(this.maxTreeDepth > 0) {
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