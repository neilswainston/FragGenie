package uk.ac.liverpool.metfrag;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
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

	enum Headers {
		METFRAG_MZ,
		METFRAG_FORMULAE,
		METFRAG_BROKEN_BONDS
	}
	
	/**
	 * 
	 */
	private final static float[] ION_MASS_CORRECTIONS = new float[] { 1.00728f };

	/**
	 * 
	 * @param smiles
	 * @param maximumTreeDepth
	 * @return double[]
	 * @throws Exception
	 */
	public static Object[] getFragmentData(final String smiles, final int maximumTreeDepth, final List<String> fields) throws Exception {
		final Fragmenter fragmenter = new Fragmenter(smiles);
		final Collection<Fragment> fragments = fragmenter.getFragments(maximumTreeDepth);
		final boolean getMasses = fields.indexOf(Headers.METFRAG_MZ.name()) != -1;
		final boolean getFormulae = fields.indexOf(Headers.METFRAG_FORMULAE.name()) != -1;
		final boolean getBonds = fields.indexOf(Headers.METFRAG_BROKEN_BONDS.name()) != -1;
		
		final List<Float> masses = getMasses ? new ArrayList<>() : null;
		final List<String> formulae = getFormulae ? new ArrayList<>() : null;
		
		// @SuppressWarnings({ "unchecked", "rawtypes" })
		final List<Collection<List<Object>>> brokenBonds = getBonds ? new ArrayList<>() : null;
		
		for (final Fragment fragment : fragments) {
			for (float ionMassCorrection : ION_MASS_CORRECTIONS) {
				if(masses != null) {
					masses.add(Float.valueOf(fragment.getMonoisotopicMass() + ionMassCorrection));
				}
			}
			
			if(formulae != null) {
				formulae.add(fragment.getFormula());
			}
			
			if(brokenBonds != null) {
				brokenBonds.add(fragment.getBrokenBonds());
			}
		}
		
		final Object[] data = new Object[3];
		data[Headers.METFRAG_MZ.ordinal()] = masses;
		data[Headers.METFRAG_FORMULAE.ordinal()] = formulae;
		data[Headers.METFRAG_BROKEN_BONDS.ordinal()] = brokenBonds;
		return data;
	}

	/**
	 * 
	 * @param inFile
	 * @param outFile
	 * @throws Exception
	 */
	private static void fragment(final File inFile, final File outFile, final String smilesHeader,
			final List<String> fields, final int maxLenSmiles, final int maxRecords) throws Exception {

		outFile.getParentFile().mkdirs();
		outFile.createNewFile();

		try (final InputStreamReader input = new InputStreamReader(new FileInputStream(inFile));
				final CSVParser csvParser = CSVFormat.DEFAULT.withFirstRecordAsHeader().parse(input);
				final CSVPrinter csvPrinter = new CSVPrinter(new OutputStreamWriter(new FileOutputStream(outFile)),
						CSVFormat.DEFAULT)) {

			final List<String> headerNames = new ArrayList<>(csvParser.getHeaderNames());
			headerNames.addAll(fields);

			csvPrinter.printRecord(headerNames);

			int count = 0;

			for (CSVRecord record : csvParser) {
				final String smiles = record.get(smilesHeader);

				if (smiles.length() < maxLenSmiles) {
					try {
						final Object[] fragmentData = getFragmentData(smiles, 2, fields);
						final Map<String, String> recordMap = record.toMap();
						
						for(final String field : fields) {
							final int idx = Headers.valueOf(field).ordinal();
							recordMap.put(field, fragmentData[idx].toString());
						}
						
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
		
		final List<String> headers = new ArrayList<>();
		
		for(final Headers header : Headers.values()) {
			headers.add(header.name());
		}
		
		fragment(new File(new File(args[0]).getAbsolutePath()),
				new File(new File(args[1]).getAbsolutePath()),
				args[2],
				headers,
				maxLenSmiles,
				maxRecords);
	}
}
