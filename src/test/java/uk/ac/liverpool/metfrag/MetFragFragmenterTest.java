/**
 * 
 */
package uk.ac.liverpool.metfrag;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
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
	public void testGetFragmentsPrecursor() throws Exception {
		final double[] expected = { 15.02294, 16.03077, 18.010021, 30.04643, 32.02568, 47.049168 };
		doTest("CCO", expected, null); //$NON-NLS-1$
	}
	
	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testGetFragmentsFilterCC() throws Exception {
		final double[] expected = { 47.049168, 16.03077, 32.02568 };
		
		final List<List<Object>> brokenBondsFilter = new ArrayList<>();
		brokenBondsFilter.add(Arrays.asList(new Object[] {Arrays.asList(new String[] {"C", "C"})})); //$NON-NLS-1$ //$NON-NLS-2$
		brokenBondsFilter.add(null);
		brokenBondsFilter.add(null);
		
		doTest("CCO", expected, brokenBondsFilter); //$NON-NLS-1$
	}
	
	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testGetFragmentsFilterSingle() throws Exception {
		final double[] expected = { 29.0386 };
		
		final List<List<Object>> brokenBondsFilter = new ArrayList<>();
		brokenBondsFilter.add(null);
		brokenBondsFilter.add(Arrays.asList(new Object[] {"SINGLE"})); //$NON-NLS-1$
		brokenBondsFilter.add(null);
		
		doTest("C=C", expected, brokenBondsFilter); //$NON-NLS-1$
	}

	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testGetFragmentsAromatic() throws Exception {
		final double[] expected = { 27.02294, 40.03077, 53.0386, 66.04643, 79.05426, 14.01511 };
		doTest("C1=CC=CC=C1", expected, null); //$NON-NLS-1$
	}
	
	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testGetFragmentsAromaticFiltered() throws Exception {
		final double[] expected = { 79.05426 };
		
		final List<List<Object>> brokenBondsFilter = new ArrayList<>();
		brokenBondsFilter.add(null);
		brokenBondsFilter.add(null);
		brokenBondsFilter.add(Arrays.asList(new Object[] {Boolean.FALSE}));
		
		doTest("C1=CC=CC=C1", expected, brokenBondsFilter); //$NON-NLS-1$
	}
	
	/**
	 * 
	 * @param smiles
	 * @param expected
	 * @throws Exception
	 */
	@SuppressWarnings("unchecked")
	private static void doTest(final String smiles, final double[] expected, final List<List<Object>> brokenBondsFilter) throws Exception {
		final List<String> headers = Arrays.asList(new String[] {
				MetFragFragmenter.Headers.METFRAG_MZ.name(),
				MetFragFragmenter.Headers.METFRAG_FORMULAE.name(),
				MetFragFragmenter.Headers.METFRAG_BROKEN_BONDS.name()
				});
		final Object[] data = MetFragFragmenter.getFragmentData(smiles, 2, headers, brokenBondsFilter);
		final List<Float> fragments = (List<Float>)data[MetFragFragmenter.Headers.METFRAG_MZ.ordinal()];
		final Float[] uniqueFrags = new HashSet<>(fragments).toArray(new Float[0]);
		
		final double[] returned = new double[uniqueFrags.length];
		int idx = 0;
		
		for(final Float value: uniqueFrags) {
			returned[idx++] = value.doubleValue();
		}
		
		doCompare(returned, expected);
		doCompare(expected, returned);
	}
	
	/**
	 * 
	 * @param returned
	 * @param expected
	 */
	private static void doCompare(final double[] returned, final double[] expected) {
		final float epsilon = 1e-5f;

		for (double mass : expected) {
			final DoubleStream ds = IntStream.range(0, returned.length).mapToDouble(i -> returned[i]);
			Assert.assertTrue(ds.anyMatch(x -> x > mass - epsilon && x < mass + epsilon));
		}
	}
}