/**
 * 
 */
package uk.ac.liverpool.metfrag;

import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Level;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import de.ipbhalle.metfraglib.interfaces.ICandidate;
import de.ipbhalle.metfraglib.list.CandidateList;
import de.ipbhalle.metfraglib.parameter.VariableNames;
import de.ipbhalle.metfraglib.process.CombinedMetFragProcess;
import de.ipbhalle.metfraglib.settings.MetFragGlobalSettings;
import uk.ac.liverpool.metfraglib.candidate.PrecursorCandidate;
import uk.ac.liverpool.metfraglib.fragmenterassignerscorer.TopDownFragmenterAssignerScorer;

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
	 * @param maximumTreeDepth
	 * @return double[]
	 * @throws Exception
	 */
	public static double[] getFragments(final String smiles, final int maximumTreeDepth) throws Exception {
		// final double PROTON_MASS = 1.00727647;
		// final SmilesParser parser = new
		// SmilesParser(SilentChemObjectBuilder.getInstance());
		// final IAtomContainer molecule = parser.parseSmiles(smiles);

		// Get settings:
		final IAtomContainer[] fragments = generateAllFragments(smiles, maximumTreeDepth);
		final Collection<Double> massesSet = new TreeSet<>();

		for (int i = 0; i < fragments.length; i++) {
			IAtomContainer fragment = fragments[i];
			final IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(fragment);
			massesSet.add(Double.valueOf(MolecularFormulaManipulator.getMajorIsotopeMass(formula)));
		}

		final double[] masses = new double[massesSet.size()];
		int i = 0;

		for (final Double mass : massesSet) {
			masses[i++] = mass.doubleValue();
		}

		return masses;
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

	/**
	 * 
	 * @param smiles
	 * @param maximumTreeDepth
	 * @return IAtomContainer[]
	 * @throws Exception
	 */
	private static IAtomContainer[] generateAllFragments(final String smiles, final int maximumTreeDepth)
			throws Exception {
		final PrecursorCandidate candidate = new PrecursorCandidate(smiles);

		final TopDownFragmenterAssignerScorer scorer = new TopDownFragmenterAssignerScorer(candidate);
		scorer.calculate();

		final IAtomContainer[] fragments = new IAtomContainer[0];

		return fragments;
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
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		final double[] result = MetFrag.getFragments("C(C(=O)O)OC1=NC(=C(C(=C1Cl)N)Cl)F", 2); //$NON-NLS-1$
		System.out.println(result);
	}
}