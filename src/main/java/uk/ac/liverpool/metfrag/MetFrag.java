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
	 * @return Collection<Map<String,Object>>
	 * @throws Exception
	 */
	public static Collection<Map<String, Object>> match(final String[] inchis, final float[] mz, final int[] inten)
			throws Exception {
		final File candidateListFile = writeCandidateList(inchis);
		final File peakListFile = writePeakList(mz, inten);

		final MetFragGlobalSettings settings = new MetFragGlobalSettings();

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
			properties.remove("empty");
			results.add(properties);
		}

		return results;
	}

	/**
	 * 
	 * @param inchis
	 * @return File
	 * @throws IOException
	 */
	private static File writeCandidateList(final String[] inchis) throws IOException {
		final File temp = File.createTempFile("candidateList", ".tmp");

		try (final PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(temp)))) {
			// Write header:
			writer.println("empty,Identifier,InChI");

			// Write data:
			for (int i = 0; i < inchis.length; i++) {
				writer.println("," + i + ",\"" + inchis[i] + "\"");
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