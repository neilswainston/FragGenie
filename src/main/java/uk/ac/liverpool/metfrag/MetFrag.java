/**
 * 
 */
package uk.ac.liverpool.metfrag;

import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;
import java.util.TreeMap;
import java.util.Map.Entry;

import org.apache.log4j.Level;

import de.ipbhalle.metfraglib.interfaces.ICandidate;
import de.ipbhalle.metfraglib.list.CandidateList;
import de.ipbhalle.metfraglib.parameter.VariableNames;
import de.ipbhalle.metfraglib.process.CombinedMetFragProcess;
import de.ipbhalle.metfraglib.settings.MetFragGlobalSettings;
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
		final Map<String, Float> formulaToMasses = getFormulaToMasses(Precursor.fromSmiles(smiles), maximumTreeDepth);
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
	private static Map<String, Float> getFormulaToMasses(final Precursor precursor, final int maxTreeDepth) throws Exception {
		final Map<String, Float> formulaToMasses = new TreeMap<>();
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

	/**
	 * 
	 * @return MetFragGlobalSettings
	 */
	private static MetFragGlobalSettings getSettings() {
		final MetFragGlobalSettings settings = new MetFragGlobalSettings();

		// Set logging:
		settings.set(VariableNames.LOG_LEVEL_NAME, Level.ALL);

		// Use SMILES:
		settings.set(VariableNames.USE_SMILES_NAME, Boolean.TRUE);

		// Remove pre-filter:
		settings.set(VariableNames.METFRAG_PRE_PROCESSING_CANDIDATE_FILTER_NAME, new String[0]);

		return settings;
	}

	/**
	 * 
	 * @param smiles
	 * @param mz
	 * @param inten
	 * @return Collection<Map<String, Object>>
	 * @throws Exception
	 */
	public static Collection<Map<String, Object>> match(final String[] smiles, final double[] mz, final int[] inten)
			throws Exception {
		final String candidateList = writeCandidateList(smiles);
		final String peakList = writePeakList(mz, inten);

		// Get settings:
		final MetFragGlobalSettings settings = getSettings();

		// Set peak list and candidate list:
		settings.set(VariableNames.PEAK_LIST_PATH_NAME, peakList);
		settings.set(VariableNames.LOCAL_DATABASE_PATH_NAME, candidateList);
		settings.set(VariableNames.METFRAG_DATABASE_TYPE_NAME, "LocalCSV"); //$NON-NLS-1$

		// Set other parameters:
		settings.set(VariableNames.PRECURSOR_NEUTRAL_MASS_NAME, Double.valueOf(mz[mz.length - 1]));

		final CombinedMetFragProcess metfragProcess = new CombinedMetFragProcess(settings);
		metfragProcess.retrieveCompounds();
		metfragProcess.run();

		final CandidateList scoredCandidateList = metfragProcess.getCandidateList();

		final Collection<Map<String, Object>> results = new ArrayList<>();

		for (int i = 0; i < scoredCandidateList.getNumberElements(); i++) {
			final ICandidate candidate = scoredCandidateList.getElement(i);
			final Map<String, Object> properties = candidate.getProperties();

			// Get index:
			final String identifier = (String) properties.remove("Identifier"); //$NON-NLS-1$
			properties.put("index", Integer.valueOf(Integer.parseInt(identifier.split("|")[0]))); //$NON-NLS-1$ //$NON-NLS-2$

			// Remove unnecessary fields:
			properties.remove("InChI"); //$NON-NLS-1$
			properties.remove("empty"); //$NON-NLS-1$

			results.add(properties);
		}

		return results;
	}

	/**
	 * 
	 * @param smiles
	 * @return String
	 * @throws IOException
	 */
	private static String writeCandidateList(final String[] smiles) throws IOException {
		try (final StringWriter writer = new StringWriter()) {
			// Write header:
			writer.write("empty," + VariableNames.IDENTIFIER_NAME + "," + VariableNames.SMILES_NAME + "," //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
					+ VariableNames.INCHI_NAME);
			writer.write(System.lineSeparator());

			// Write data:
			for (int i = 0; i < smiles.length; i++) {
				writer.write("," + i + "," + smiles[i] + ",\"INVALID_INCHI\""); //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
				writer.write(System.lineSeparator());
			}

			return writer.toString();
		}
	}

	/**
	 * 
	 * @param mz
	 * @param inten
	 * @return String
	 * @throws IOException
	 */
	private static String writePeakList(final double[] mz, final int[] inten) throws IOException {
		// Ensure length of mz and inten are equal.
		assert mz.length == inten.length;

		try (final StringWriter writer = new StringWriter()) {
			for (int i = 0; i < mz.length; i++) {
				writer.write(mz[i] + " " + inten[i]); //$NON-NLS-1$
				writer.write(System.lineSeparator());
			}

			return writer.toString();
		}
	}
}