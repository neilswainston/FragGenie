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
		final float[] expected = { 15.02294f, 16.03077f, 18.010021f, 30.04643f, 32.02568f, 47.049168f };
		doTest("CCO", expected, null); //$NON-NLS-1$
	}
	
	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testGetFragmentsFilterCC() throws Exception {
		final float[] expected = { 47.049168f, 16.03077f, 32.02568f };
		
		final List<List<Object>> brokenBondsFilter = new ArrayList<>();
		brokenBondsFilter.add(Arrays.asList(new Object[] {Arrays.asList(new String[] {"C", "C"}), null, null})); //$NON-NLS-1$ //$NON-NLS-2$
		
		doTest("CCO", expected, brokenBondsFilter); //$NON-NLS-1$
	}
	
	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testGetFragmentsFilterSingle() throws Exception {
		final float[] expected = { 29.0386f };
		
		final List<List<Object>> brokenBondsFilter = new ArrayList<>();
		brokenBondsFilter.add(Arrays.asList(new Object[] {null, "SINGLE", null})); //$NON-NLS-1$
		
		doTest("C=C", expected, brokenBondsFilter); //$NON-NLS-1$
	}

	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testGetFragmentsAromatic() throws Exception {
		final float[] expected = { 27.02294f, 40.03077f, 53.0386f, 66.04643f, 79.05426f, 14.01511f };
		doTest("C1=CC=CC=C1", expected, null); //$NON-NLS-1$
	}
	
	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testGetFragmentsAromaticFiltered() throws Exception {
		final float[] expected = { 79.05426f };
		
		final List<List<Object>> brokenBondsFilter = new ArrayList<>();
		brokenBondsFilter.add(Arrays.asList(new Object[] {null, null, Boolean.FALSE}));
		
		doTest("C1=CC=CC=C1", expected, brokenBondsFilter); //$NON-NLS-1$
	}
	
	/**
	 * 
	 * @param smiles
	 * @param expected
	 * @throws Exception
	 */
	@SuppressWarnings("unchecked")
	private static void doTest(final String smiles, final float[] expected, final List<List<Object>> brokenBondsFilter) throws Exception {
		final List<String> headers = Arrays.asList(new String[] {
				MetFragFragmenter.Headers.METFRAG_MZ.name(),
				MetFragFragmenter.Headers.METFRAG_FORMULAE.name(),
				MetFragFragmenter.Headers.METFRAG_BROKEN_BONDS.name()
				});
		final Object[] data = MetFragFragmenter.getFragmentData(smiles, 2, headers, brokenBondsFilter);
		final List<Float> fragments = (List<Float>)data[MetFragFragmenter.Headers.METFRAG_MZ.ordinal()];
		final Float[] uniqueFrags = new HashSet<>(fragments).toArray(new Float[0]);
		
		final float[] returned = new float[uniqueFrags.length];
		int idx = 0;
		
		for(final Float value: uniqueFrags) {
			returned[idx++] = value.floatValue();
		}
		
		doCompare(returned, expected);
		doCompare(expected, returned);
	}
	
	/**
	 * 
	 * @param returned
	 * @param expected
	 */
	private static void doCompare(final float[] returned, final float[] expected) {
		final float epsilon = 1e-5f;

		for (double mass : expected) {
			final DoubleStream ds = IntStream.range(0, returned.length).mapToDouble(i -> returned[i]);
			Assert.assertTrue(ds.anyMatch(x -> x > mass - epsilon && x < mass + epsilon));
		}
	}
}