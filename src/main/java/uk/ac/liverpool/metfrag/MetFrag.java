/**
 * 
 */
package uk.ac.liverpool.metfrag;

import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;
import java.util.TreeMap;
import java.util.Map.Entry;

import uk.ac.liverpool.metfraglib.Fragment;
import uk.ac.liverpool.metfraglib.Fragmenter;
import uk.ac.liverpool.metfraglib.Precursor;

/**
 * 
 * @author neilswainston
 */
public class MetFrag {

	/**
	 * 
	 * @param smiles
	 * @param maximumTreeDepth
	 * @return double[]
	 * @throws Exception
	 */
	public static float[] getFragmentMasses(final String smiles, final int maximumTreeDepth) throws Exception {
		final Map<String, Float> formulaToMasses = getFormulaToMasses(smiles, maximumTreeDepth);
		final float[] ionMassCorrections = new float[] { 1.00728f };
		final float[] correctedMasses = new float[formulaToMasses.size() * ionMassCorrections.length];
		int i = 0;
		
		for(final Entry<String, Float> entry : formulaToMasses.entrySet()) {
			for (float ionMassCorrection : ionMassCorrections) {
				correctedMasses[i++] = entry.getValue().floatValue() + ionMassCorrection;
			}
		}
		
		return correctedMasses;
	}
	
	/**
	 * 
	 * @return Map<String, Float>
	 * @throws Exception
	 */
	private static Map<String, Float> getFormulaToMasses(final String smiles, final int maxTreeDepth) throws Exception {
		final Map<String, Float> formulaToMasses = new TreeMap<>();
		final Precursor precursor = Precursor.fromSmiles(smiles);
		final Fragmenter fragmenter = new Fragmenter(precursor);
		Queue<Fragment> fragments = new LinkedList<>();
		
		final Fragment precursorFragment = new Fragment(precursor);
		fragments.add(precursorFragment);
		
		formulaToMasses.put(precursorFragment.getFormula(), Float.valueOf(precursorFragment.getMonoisotopicMass()));
		
		for (int k = 1; k <= maxTreeDepth; k++) {
			Queue<Fragment> newFragments = new LinkedList<>();

			while (!fragments.isEmpty()) {
				final Fragment fragment = fragments.poll();

				for (final Fragment childFragment : fragmenter.getFragmentsOfNextTreeDepth(fragment)) {
					formulaToMasses.put(childFragment.getFormula(), Float.valueOf(childFragment.getMonoisotopicMass()));

					if (maxTreeDepth > 0) {
						newFragments.add(childFragment);
					}
				}
			}

			fragments = newFragments;
		}
		
		return formulaToMasses;
	}
}