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
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

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
	private final static float[] ION_MASS_CORRECTIONS = new float[] { -0.00055f, 1.00728f }; // { [M]+, [M]+H+

	/**
	 * 
	 * @param smiles
	 * @param maximumTreeDepth
	 * @param minMass
	 * @param fields
	 * @return
	 * @throws Exception
	 */
	public static Object[] getFragmentData(final String smiles, final int maximumTreeDepth, final float minMass, final List<String> fields) throws Exception {
		final List<List<Object>> brokenBondsFilter = null;
		return getFragmentData(smiles, maximumTreeDepth, minMass, fields, brokenBondsFilter);
	}
		
	/**
	 * 
	 * @param smiles
	 * @param maximumTreeDepth
	 * @param minMass
	 * @param fields
	 * @param brokenBondsFilter
	 * @return Object[]
	 * @throws Exception
	 */
	public static Object[] getFragmentData(final String smiles, final int maximumTreeDepth, final float minMass, final List<String> fields, final List<List<Object>> brokenBondsFilter) throws Exception {
		final Fragment fragment = new Fragment(new Precursor(getAtomContainer(smiles)));
		final Collection<Fragment> fragments = fragment.fragment(maximumTreeDepth, minMass);
		final boolean getMasses = fields.indexOf(Headers.METFRAG_MZ.name()) != -1;
		final boolean getFormulae = fields.indexOf(Headers.METFRAG_FORMULAE.name()) != -1;
		final boolean getBonds = fields.indexOf(Headers.METFRAG_BROKEN_BONDS.name()) != -1;
		
		final List<Float> masses = getMasses ? new ArrayList<>() : null;
		final List<String> formulae = getFormulae ? new ArrayList<>() : null;
		final List<Collection<Bond>> brokenBonds = getBonds ? new ArrayList<>() : null;
		
		for (final Fragment childFragment : fragments) {
			final Collection<Bond> fragBrokenBonds = childFragment.getBrokenBonds();
			
			if(filter(fragBrokenBonds, brokenBondsFilter)) {
				for (float ionMassCorrection : ION_MASS_CORRECTIONS) {
					if(masses != null) {
						masses.add(Float.valueOf(childFragment.getMonoisotopicMass() + ionMassCorrection));
					}
					
					if(formulae != null) {
						formulae.add(childFragment.getFormula());
					}
					
					if(brokenBonds != null) {
						brokenBonds.add(fragBrokenBonds);
					}
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
	 * @param smilesHeader
	 * @param fields
	 * @param maximumTreeDepth
	 * @param minMass
	 * @param maxLenSmiles
	 * @param maxRecords
	 * @throws Exception
	 */
	public static void fragment(final File inFile, final File outFile, final String smilesHeader,
			final List<String> fields, final int maximumTreeDepth, final float minMass, 
			final int maxLenSmiles, final int maxRecords) throws Exception {
		final List<List<Object>> brokenBondsFilter = null;
		
		fragment(inFile, outFile, smilesHeader, fields, maximumTreeDepth, minMass, brokenBondsFilter, maxLenSmiles, maxRecords);
	}

	/**
	 * 
	 * @param inFile
	 * @param outFile
	 * @param smilesHeader
	 * @param fields
	 * @param maximumTreeDepth
	 * @param minMass
	 * @param brokenBondsFilter
	 * @param maxLenSmiles
	 * @param maxRecords
	 * @throws Exception
	 */
	public static void fragment(final File inFile, final File outFile, final String smilesHeader,
			final List<String> fields, final int maximumTreeDepth, final float minMass, final List<List<Object>> brokenBondsFilter,
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
						final Object[] fragmentData = getFragmentData(smiles, maximumTreeDepth, minMass, fields, brokenBondsFilter);
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
	private static boolean filter(final Collection<Bond> fragBrokenBonds, final List<List<Object>> brokenBondsFilter) {
		if(brokenBondsFilter != null) {
			for(final Bond fragBrokenBond : fragBrokenBonds) {
				for(int i = 0; i < brokenBondsFilter.size(); i++) {
					final List<Object> filter = brokenBondsFilter.get(i);
					
					if(!fragBrokenBond.passesFilter(filter)) {
						return false;
					}
				}
			}
		}
		
		return true;
	}
	
	/**
	 * 
	 * @param smiles
	 * @return IAtomContainer
	 * @throws CDKException 
	 * @throws Exception
	 */
	private static IAtomContainer getAtomContainer(final String smiles) throws CDKException {
		final SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
		final IAtomContainer molecule = parser.parseSmiles(smiles);
		final Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());

		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		aromaticity.apply(molecule);

		final CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(molecule.getBuilder());

		for (int i = 0; i < molecule.getAtomCount(); i++) {
			hAdder.addImplicitHydrogens(molecule, molecule.getAtom(i));
		}

		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);

		removeHydrogens(molecule);
		return molecule;
	}

	/**
	 * 
	 * @param molecule
	 */
	private static void removeHydrogens(final IAtomContainer molecule) {
		final Collection<IAtom> hydrogenAtoms = new ArrayList<>();

		for (IAtom atom : molecule.atoms()) {
			if (atom.getSymbol().equals("H")) { //$NON-NLS-1$
				hydrogenAtoms.add(atom);
			}

			int numberHydrogens = 0;

			for (IAtom neighbour : molecule.getConnectedAtomsList(atom)) {
				if (neighbour.getSymbol().equals("H")) { //$NON-NLS-1$
					numberHydrogens++;
				}
			}

			atom.setImplicitHydrogenCount(Integer.valueOf(numberHydrogens));
		}

		for (IAtom atom : hydrogenAtoms) {
			molecule.removeAtom(atom);
		}
	}

	/**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		final File inFile = new File(new File(args[0]).getAbsolutePath());
		final File outFile = new File(new File(args[1]).getAbsolutePath());
		final String smilesHeader = args[2];
		final int maximumTreeDepth = Integer.parseInt(args[3]);
		final float minMass = Float.parseFloat(args[4]);
		final int maxLenSmiles = Integer.parseInt(args[5]);
		final int maxRecords = Integer.parseInt(args[6]);
		final List<String> fields = Arrays.asList(Arrays.copyOfRange(args, 7, args.length));
		
		// final List<List<Object>> brokenBondsFilter = null;
		final List<List<Object>> brokenBondsFilter = new ArrayList<>();
		brokenBondsFilter.add(Arrays.asList(new Object[] {null, "SINGLE", Boolean.FALSE})); //$NON-NLS-1$
		
		fragment(inFile,
				outFile,
				smilesHeader,
				fields,
				maximumTreeDepth,
				minMass,
				brokenBondsFilter,
				maxLenSmiles,
				maxRecords);
	}
}
