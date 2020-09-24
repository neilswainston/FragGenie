package uk.ac.liverpool.metfrag;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVPrinter;
import org.apache.commons.csv.CSVRecord;

import uk.ac.liverpool.metfraglib.Fragment;
import uk.ac.liverpool.metfraglib.Fragmenter;

/**
 * 
 * @author neilswainston
 */
public class MetFragFragmenter {

	/**
	 * 
	 */
	private final static String METFRAG_HEADER = "MetFrag m/z"; //$NON-NLS-1$

	/**
	 * 
	 * @param smiles
	 * @param maximumTreeDepth
	 * @return double[]
	 * @throws Exception
	 */
	public static float[] getFragmentData(final String smiles, final int maximumTreeDepth) throws Exception {
		final Fragmenter fragmenter = new Fragmenter(smiles);
		final Collection<Fragment> fragments = fragmenter.getFragments(maximumTreeDepth);
		final float[] ionMassCorrections = new float[] { 1.00728f };
		final float[] correctedMasses = new float[fragments.size() * ionMassCorrections.length];
		int i = 0;

		for (final Fragment fragment : fragments) {
			for (float ionMassCorrection : ionMassCorrections) {
				correctedMasses[i++] = fragment.getMonoisotopicMass() + ionMassCorrection;
			}
		}

		return correctedMasses;
	}

	/**
	 * 
	 * @param inFile
	 * @param outFile
	 * @throws Exception
	 */
	private static void fragment(final File inFile, final File outFile, final String smilesHeader,
			final int maxLenSmiles, final int maxRecords) throws Exception {

		outFile.getParentFile().mkdirs();
		outFile.createNewFile();

		try (final InputStreamReader input = new InputStreamReader(new FileInputStream(inFile));
				final CSVParser csvParser = CSVFormat.DEFAULT.withFirstRecordAsHeader().parse(input);
				final CSVPrinter csvPrinter = new CSVPrinter(new OutputStreamWriter(new FileOutputStream(outFile)),
						CSVFormat.DEFAULT)) {

			final List<String> headerNames = new ArrayList<>(csvParser.getHeaderNames());
			headerNames.add(METFRAG_HEADER);

			csvPrinter.printRecord(headerNames);

			int count = 0;

			for (CSVRecord record : csvParser) {
				final String smiles = record.get(smilesHeader);

				if (smiles.length() < maxLenSmiles) {
					try {
						final float[] fragments = getFragmentData(smiles, 2);
						final Map<String, String> recordMap = record.toMap();
						recordMap.put(METFRAG_HEADER, Arrays.toString(fragments));

						final List<String> values = new ArrayList<>();

						for (String headerName : headerNames) {
							values.add(recordMap.get(headerName));
						}

						csvPrinter.printRecord(values);
					} catch (Exception e) {
						e.printStackTrace();
					}
				}

				if (count % 100 == 0) {
					System.out.println("Records fragmented: " + Integer.toString(count)); //$NON-NLS-1$
				}

				if (count++ == maxRecords) {
					break;
				}
			}
		}
	}

	/**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		final int maxLenSmiles = args.length > 3 ? Integer.parseInt(args[3]) : Integer.MAX_VALUE;
		final int maxRecords = args.length > 4 ? Integer.parseInt(args[4]) : Integer.MAX_VALUE;

		fragment(new File(new File(args[0]).getAbsolutePath()), new File(new File(args[1]).getAbsolutePath()), args[2],
				maxLenSmiles, maxRecords);
	}
}
