/**
 * 
 */
package uk.ac.liverpool.metfrag;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import org.junit.Assert;
import org.junit.Test;

/**
 * 
 * @author neilswainston
 */
public class MetFragFragmenterTest {

	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testGetFragments() throws Exception {
		final float[] expected = {
				90.97453f,
				106.94498f,
				// 115.98977f, // [C4H4ClNO-H]+
				117.98543f, // [C4H2ClFN]+
				133.95588f, // [C4H2Cl2N]+
				143.9885f, // [C5H2ClFN2]+
				144.99633f,
				// 146.00416f, // [C5H2ClFN2+H]+H+
				// 151.94645f, // [C4H2Cl2FN-H]+
				160.96678f,
				163.0069f,
				// 172.99124f, // [C6H4ClFN2O-H]+
				178.95735f, // [C5H2Cl2FN2]+
				// 178.97735f, // [C5H2Cl2N2O+2H]+H+
				// 180.97301f, // [C5H2Cl2FN2+H]+H+
				// 196.96792f, // [C5H2Cl2FN2O+H]+H+
				208.96792f, // [C6H4Cl2FN2O]+
				236.96283f // [C7H4Cl2FN2O2]+
		};

		doTest("C(C(=O)O)OC1=NC(=C(C(=C1Cl)N)Cl)F", expected, null, 3, 0.0f, false); //$NON-NLS-1$
	}
	
	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testGetFragmentsPrecursor() throws Exception {
		final float[] expected = { 14.01511f, 15.02294f, 15.022941f, 16.03077f, 17.002192f, 18.010021f, 29.0386f, 30.04643f, 31.017853f, 32.02568f, 46.041344f, 47.04917f };
		doTest("CCO", expected, null, 2, 0.0f, false); //$NON-NLS-1$
	}
	
	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testGetFragmentsFilterCC() throws Exception {
		final float[] expected = { 31.017853f, 32.02568f, 46.041344f, 47.04917f };
		
		final List<List<Object>> brokenBondsFilter = new ArrayList<>();
		brokenBondsFilter.add(Arrays.asList(new Object[] {Arrays.asList(new String[] {"C", "C"}), null, null})); //$NON-NLS-1$ //$NON-NLS-2$
		
		doTest("CCO", expected, brokenBondsFilter, 2, 20.0f, true); //$NON-NLS-1$
	}
	
	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testGetFragmentsFilterSingle() throws Exception {
		final float[] expected = { 28.030771f, 29.0386f };
		
		final List<List<Object>> brokenBondsFilter = new ArrayList<>();
		brokenBondsFilter.add(Arrays.asList(new Object[] {null, "SINGLE", null})); //$NON-NLS-1$
		
		doTest("C=C", expected, brokenBondsFilter, 2, 0.0f, true); //$NON-NLS-1$
	}

	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testGetFragmentsAromatic() throws Exception {
		final float[] expected = { 13.007279f, 14.01511f, 26.01511f, 27.02294f, 39.02294f, 40.03077f, 52.03077f, 53.038597f, 65.0386f, 66.046425f, 78.046425f, 79.05425f };
		doTest("C1=CC=CC=C1", expected, null, 2, 0f, true); //$NON-NLS-1$
	}
	
	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testGetFragmentsAromaticFiltered() throws Exception {
		final float[] expected = { 78.046425f, 79.05425f };
		
		final List<List<Object>> brokenBondsFilter = new ArrayList<>();
		brokenBondsFilter.add(Arrays.asList(new Object[] {null, null, Boolean.FALSE}));
		
		doTest("C1=CC=CC=C1", expected, brokenBondsFilter, 2, 0.0f, true); //$NON-NLS-1$
	}
	
	/**
	 * 
	 * @param smiles
	 * @param expected
	 * @param brokenBondsFilter
	 * @param maximumTreeDepth
	 * @param minMass
	 * @param biDirectional
	 * @throws Exception
	 */
	@SuppressWarnings("unchecked")
	private static void doTest(final String smiles, final float[] expected, final List<List<Object>> brokenBondsFilter, final int maximumTreeDepth, final float minMass, final boolean biDirectional) throws Exception {
		final List<String> headers = Arrays.asList(new String[] {
				MetFragFragmenter.Headers.METFRAG_MZ.name(),
				MetFragFragmenter.Headers.METFRAG_FORMULAE.name(),
				MetFragFragmenter.Headers.METFRAG_BROKEN_BONDS.name()
				});
		final Object[] data = MetFragFragmenter.getFragmentData(smiles, maximumTreeDepth, minMass, headers, brokenBondsFilter);
		final List<Float> fragments = (List<Float>)data[MetFragFragmenter.Headers.METFRAG_MZ.ordinal()];
		final Float[] uniqueFrags = new TreeSet<>(fragments).toArray(new Float[0]);
		
		final float[] returned = new float[uniqueFrags.length];
		int idx = 0;
		
		for(final Float value: uniqueFrags) {
			returned[idx++] = value.floatValue();
		}
		
		doCompare(returned, expected);
		
		if(biDirectional) {
			doCompare(expected, returned);
		}
	}
	
	/**
	 * 
	 * @param returned
	 * @param expected
	 */
	private static void doCompare(final float[] returned, final float[] expected) {
		final float epsilon = 1e-4f;

		for (double mass : expected) {
			final DoubleStream ds = IntStream.range(0, returned.length).mapToDouble(i -> returned[i]);
			Assert.assertTrue(ds.anyMatch(x -> x > mass - epsilon && x < mass + epsilon));
		}
	}
}