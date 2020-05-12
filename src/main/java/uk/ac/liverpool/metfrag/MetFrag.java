/**
 * 
 */
package uk.ac.liverpool.metfrag;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
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
	 * @param inchis
	 * @param mz
	 * @param inten
	 * @return Collection<Map<String, Object>>
	 * @throws Exception
	 */
	public static Collection<Map<String, Object>> match(final String[] smiles, final String[] inchis, final float[] mz, final int[] inten)
			throws Exception {
		final File candidateListFile = writeCandidateList(smiles, inchis);
		final File peakListFile = writePeakList(mz, inten);

		final MetFragGlobalSettings settings = new MetFragGlobalSettings();

		// Set logging:
		settings.set(VariableNames.LOG_LEVEL_NAME, Level.ALL);
		
		// Set peak list and candidate list:
		settings.set(VariableNames.PEAK_LIST_PATH_NAME, peakListFile.getAbsolutePath());
		settings.set(VariableNames.LOCAL_DATABASE_PATH_NAME, candidateListFile.getAbsolutePath());
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
			
			// Map Identifier to index:
			final String identifier = (String)properties.remove(VariableNames.IDENTIFIER_NAME);
			properties.put("index", Integer.parseInt(identifier));
			
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
	 * @return File
	 * @throws IOException
	 */
	private static File writeCandidateList(final String[] smiles, final String[] inchis) throws IOException {
		// Ensure length of smiles and inchis are equal.
		assert smiles.length == inchis.length;
		
		final File temp = File.createTempFile("candidateList", ".tmp");

		try (final PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(temp)))) {
			// Write header:
			writer.println("empty," + VariableNames.IDENTIFIER_NAME + "," + VariableNames.SMILES_NAME + "," + VariableNames.INCHI_NAME);

			// Write data:
			for (int i = 0; i < smiles.length; i++) {
				writer.println("," + i + "," + smiles[i] + ",\"" + inchis[i] + "\"");
			}
		}

		return temp;
	}

	/**
	 * 
	 * @param mz
	 * @param inten
	 * @return File
	 * @throws IOException
	 */
	private static File writePeakList(final float[] mz, final int[] inten) throws IOException {
		// Ensure length of mz and inten are equal.
		assert mz.length == inten.length;

		final File temp = File.createTempFile("peakList", ".tmp");

		try (final PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(temp)))) {
			for (int i = 0; i < mz.length; i++) {
				writer.println(mz[i] + " " + inten[i]);
			}
		}

		return temp;
	}
}