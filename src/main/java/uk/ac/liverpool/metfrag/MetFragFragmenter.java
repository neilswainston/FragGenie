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
	public static Object[] getFragmentData(final String smiles, final int maximumTreeDepth, final List<String> fields, final List<List<Object>> brokenBondsFilter) throws Exception {
		final Fragmenter fragmenter = new Fragmenter(smiles);
		final Collection<Fragment> fragments = fragmenter.getFragments(maximumTreeDepth);
		final boolean getMasses = fields.indexOf(Headers.METFRAG_MZ.name()) != -1;
		final boolean getFormulae = fields.indexOf(Headers.METFRAG_FORMULAE.name()) != -1;
		final boolean getBonds = fields.indexOf(Headers.METFRAG_BROKEN_BONDS.name()) != -1;
		
		final List<Float> masses = getMasses ? new ArrayList<>() : null;
		final List<String> formulae = getFormulae ? new ArrayList<>() : null;
		final List<Collection<List<Object>>> brokenBonds = getBonds ? new ArrayList<>() : null;
		
		for (final Fragment fragment : fragments) {
			final Collection<List<Object>> fragBrokenBonds = fragment.getBrokenBonds();
			
			if(filter(fragBrokenBonds, brokenBondsFilter)) {
				for (float ionMassCorrection : ION_MASS_CORRECTIONS) {
					if(masses != null) {
						masses.add(Float.valueOf(fragment.getMonoisotopicMass() + ionMassCorrection));
					}
				}
				
				if(formulae != null) {
					formulae.add(fragment.getFormula());
				}
				
				if(brokenBonds != null) {
					brokenBonds.add(fragBrokenBonds);
				}
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
	public static void fragment(final File inFile, final File outFile, final String smilesHeader,
			final List<String> fields, final List<List<Object>> brokenBondsFilter,
			final int maxLenSmiles, final int maxRecords) throws Exception {

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
						final Object[] fragmentData = getFragmentData(smiles, 2, fields, brokenBondsFilter);
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
	 * @param fragBrokenBonds
	 * @param brokenBondsFilter
	 * @return boolean
	 */
	private static boolean filter(final Collection<List<Object>> fragBrokenBonds, final List<List<Object>> brokenBondsFilter) {
		for(final List<Object> fragBrokenBond : fragBrokenBonds) {
			for(int i = 0; i < fragBrokenBond.size(); i++) {
				final Object value = fragBrokenBond.get(i);
				final List<Object> filter = brokenBondsFilter.get(i);
				
				if(filter != null && !filter.contains(value)) {
					return false;
				}
			}
		}
		return true;
	}

	/**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		final int maxLenSmiles = Integer.parseInt(args[3]);
		final int maxRecords = Integer.parseInt(args[4]);
		
		final List<String> fields = Arrays.asList(Arrays.copyOfRange(args, 5, args.length));
		
		final List<List<Object>> brokenBondsFilter = new ArrayList<>();
		brokenBondsFilter.add(null);
		brokenBondsFilter.add(Arrays.asList(new Object[] {"SINGLE"})); //$NON-NLS-1$
		brokenBondsFilter.add(Arrays.asList(new Object[] {Boolean.FALSE}));
		
		fragment(new File(new File(args[0]).getAbsolutePath()),
				new File(new File(args[1]).getAbsolutePath()),
				args[2],
				fields,
				brokenBondsFilter,
				maxLenSmiles,
				maxRecords);
	}
}
