/**
 * 
 */
package uk.ac.liverpool.metfrag;

import java.util.List;
import java.util.Map;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import org.junit.Assert;
import org.junit.Test;

/**
 * 
 * @author neilswainston
 */
public class MetFragTest {

	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testGetFragments() throws Exception {
		final double[] expected = { 90.97453, 106.94498,
				// 115.98977,
				// 117.98543,
				// 133.95588,
				// 143.9885,
				144.99633,
				// 146.00416,
				151.94645, 160.96678, 163.0069,
				// 172.99124,
				// 178.95735,
				178.97735,
				// 180.97301,
				// 196.96792,
				// 208.96792,
				// 236.96283
				};

		final float[] fragments = MetFrag.getFragments("C(C(=O)O)OC1=NC(=C(C(=C1Cl)N)Cl)F", 2); //$NON-NLS-1$
		final float epsilon = 1e-5f;

		for (double mass : expected) {
			final DoubleStream ds = IntStream.range(0, fragments.length).mapToDouble(i -> fragments[i]);
			Assert.assertTrue(ds.anyMatch(x -> x > mass - epsilon && x < mass + epsilon));
		}
	}
	
	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testGetFragmentsPrecursor() throws Exception {
		final double[] expected = {
				15.02294,
				16.03077,
				18.010021,
				30.04643,
				32.02568,
				47.049168
			};

		final float[] fragments = MetFrag.getFragments("CCO", 2); //$NON-NLS-1$
		final float epsilon = 1e-5f;

		for (double mass : expected) {
			final DoubleStream ds = IntStream.range(0, fragments.length).mapToDouble(i -> fragments[i]);
			Assert.assertTrue(ds.anyMatch(x -> x > mass - epsilon && x < mass + epsilon));
		}
	}

	/**
	 * 
	 * @throws Exception
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testMatch() throws Exception {
		final List<Map<String, Object>> results = (List<Map<String, Object>>) MetFrag.match(MetFragTestData.SMILES,
				MetFragTestData.MZ, MetFragTestData.INTEN);
		Assert.assertEquals(5, results.size());
		Assert.assertEquals("C(C(=O)O)OC1=NC(=C(C(=C1Cl)N)Cl)F", results.get(0).get("SMILES")); //$NON-NLS-1$ //$NON-NLS-2$
		Assert.assertEquals(1188.6403357016206, ((Double) results.get(0).get("FragmenterScore")).doubleValue(), 1e-6); //$NON-NLS-1$
	}
}