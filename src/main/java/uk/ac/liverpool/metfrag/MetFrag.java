/**
 * 
 */
package uk.ac.liverpool.metfrag;

import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;

import org.apache.log4j.Level;

import de.ipbhalle.metfraglib.interfaces.ICandidate;
import de.ipbhalle.metfraglib.list.ScoredCandidateList;
import de.ipbhalle.metfraglib.parameter.VariableNames;
import de.ipbhalle.metfraglib.process.CombinedMetFragProcess;
import de.ipbhalle.metfraglib.settings.MetFragGlobalSettings;

/**
 * 
 * @author neilswainston
 */
public class MetFrag {

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

		final MetFragGlobalSettings settings = new MetFragGlobalSettings();

		// Set logging:
		settings.set(VariableNames.LOG_LEVEL_NAME, Level.ALL);
	
		// Use SMILES:
		settings.set(VariableNames.USE_SMILES_NAME, Boolean.TRUE);
		
		// Remove pre-filter:
		settings.set(VariableNames.METFRAG_PRE_PROCESSING_CANDIDATE_FILTER_NAME, new String[0]);
		
		// Set peak list and candidate list:
		settings.set(VariableNames.PEAK_LIST_PATH_NAME, peakList);
		settings.set(VariableNames.LOCAL_DATABASE_PATH_NAME, candidateList);
		settings.set(VariableNames.METFRAG_DATABASE_TYPE_NAME, "LocalCSV");

		// Set other parameters:
		settings.set(VariableNames.RELATIVE_MASS_DEVIATION_NAME, 5.0);
		settings.set(VariableNames.PRECURSOR_NEUTRAL_MASS_NAME, 253.966126);

		final CombinedMetFragProcess metfragProcess = new CombinedMetFragProcess(settings);

		metfragProcess.retrieveCompounds();
		metfragProcess.run();

		final ScoredCandidateList scoredCandidateList = (ScoredCandidateList) metfragProcess.getCandidateList();

		final Collection<Map<String, Object>> results = new ArrayList<>();

		for (int i = 0; i < scoredCandidateList.getNumberElements(); i++) {
			final ICandidate candidate = scoredCandidateList.getElement(i);
			final Map<String, Object> properties = candidate.getProperties();
			
			// Remove unnecessary field:
			properties.remove("empty");
			
			results.add(properties);
		}

		return results;
	}

	/**
	 * 
	 * @param smiles
	 * @param inchis
	 * @return String
	 * @throws IOException
	 */
	private static String writeCandidateList(final String[] smiles) throws IOException {
		try (final StringWriter writer = new StringWriter()) {
			// Write header:
			writer.write("empty," + VariableNames.IDENTIFIER_NAME + "," + VariableNames.SMILES_NAME + "," + VariableNames.INCHI_NAME);
			writer.write(System.lineSeparator());
			
			// Write data:
			for (int i = 0; i < smiles.length; i++) {
				writer.write("," + i + "," + smiles[i] + ",\"INVALID_INCHI\"");
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
				writer.write(mz[i] + " " + inten[i]);
				writer.write(System.lineSeparator());
			}
			
			return writer.toString();
		}
	}
}