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

import de.ipbhalle.metfraglib.candidate.TopDownPrecursorCandidate;
import de.ipbhalle.metfraglib.fragmenter.TopDownNeutralLossFragmenter;
import de.ipbhalle.metfraglib.interfaces.ICandidate;
import de.ipbhalle.metfraglib.list.FragmentList;
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

		// Get settings:
		final MetFragGlobalSettings settings = getSettings();
		
		// Set peak list and candidate list:
		settings.set(VariableNames.PEAK_LIST_PATH_NAME, peakList);
		settings.set(VariableNames.LOCAL_DATABASE_PATH_NAME, candidateList);
		settings.set(VariableNames.METFRAG_DATABASE_TYPE_NAME, "LocalCSV");

		// Set other parameters:
		settings.set(VariableNames.RELATIVE_MASS_DEVIATION_NAME, 5.0);
		settings.set(VariableNames.PRECURSOR_NEUTRAL_MASS_NAME, 253.966126);

		final CombinedMetFragProcess metfragProcess = new CombinedMetFragProcess(settings, Level.ALL);

		metfragProcess.retrieveCompounds();
		metfragProcess.run();

		final ScoredCandidateList scoredCandidateList = (ScoredCandidateList) metfragProcess.getCandidateList();

		final Collection<Map<String, Object>> results = new ArrayList<>();

		for (int i = 0; i < scoredCandidateList.getNumberElements(); i++) {
			final ICandidate candidate = scoredCandidateList.getElement(i);
			final Map<String, Object> properties = candidate.getProperties();
			
			// Get index:
			final String identifier = (String)properties.remove("Identifier");
			properties.put("index", Integer.parseInt(identifier.split("|")[0]));
			
			// Remove unnecessary fields:
			properties.remove("InChI");
			properties.remove("empty");
			
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
		// final SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
		// final IAtomContainer molecule = parser.parseSmiles(smiles);
		
		// Get settings:
		final IAtomContainer[] fragments = generateAllFragments(smiles, maximumTreeDepth);
		final Collection<Double> massesSet = new TreeSet<>();
		
		for(int i = 0; i < fragments.length; i++)  {
			IAtomContainer fragment = fragments[i];
			final IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(fragment);
			massesSet.add(MolecularFormulaManipulator.getMajorIsotopeMass(formula));
		}
		
		final double[] masses = new double[massesSet.size()];
		int i = 0;
		
		for(final Double mass : massesSet) {
			masses[i++] = mass.doubleValue();
		}
		
		return masses;
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
	
	/**
	 * 
	 * @param molecule
	 * @param maximumTreeDepth
	 * @return IAtomContainer
	 * @throws Exception 
	 */
	private static IAtomContainer[] generateAllFragments(String smiles, int maximumTreeDepth) throws Exception {
		final MetFragGlobalSettings settings = getSettings();
		final ICandidate candidate = new TopDownPrecursorCandidate(null, "IDENTIFIER", smiles);
		candidate.setUseSmiles(true);
		candidate.initialisePrecursorCandidate();

		settings.set(VariableNames.CANDIDATE_NAME, candidate);
		settings.set(VariableNames.MAXIMUM_TREE_DEPTH_NAME, (byte)2);
		settings.set(VariableNames.MINIMUM_FRAGMENT_MASS_LIMIT_NAME, 0.0);
		settings.set(VariableNames.MAXIMUM_NUMBER_OF_TOPDOWN_FRAGMENT_ADDED_TO_QUEUE, (byte)maximumTreeDepth);
		
		final TopDownNeutralLossFragmenter fragmenter = new TopDownNeutralLossFragmenter(settings);
		final FragmentList fragmentList = fragmenter.generateFragments();
		final IAtomContainer[] fragments = new IAtomContainer[fragmentList.getNumberElements()];
		
		for(int i = 0; i < fragmentList.getNumberElements(); i++) {
			fragments[i] = fragmentList.getElement(i).getStructureAsIAtomContainer(candidate.getPrecursorMolecule());
		}
		
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
		final double[] result = MetFrag.getFragments("C(C(=O)O)OC1=NC(=C(C(=C1Cl)N)Cl)F", 2);
		System.out.println(result);
	}
}