/**
 * 
 */
package uk.ac.liverpool.metfrag;

import java.util.List;
import java.util.Map;
import java.util.stream.DoubleStream;

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
	@Test
	public void testMatch() throws Exception {
		final List<Map<String, Object>> results = (List<Map<String, Object>>) MetFrag.match(MetFragTestData.SMILES,
				MetFragTestData.MZ, MetFragTestData.INTEN);
		Assert.assertEquals(5, results.size());
		Assert.assertEquals("C(C(=O)O)OC1=NC(=C(C(=C1Cl)N)Cl)F", results.get(0).get("SMILES")); //$NON-NLS-1$ //$NON-NLS-2$
		Assert.assertEquals(1188.6403357016206, ((Double) results.get(0).get("FragmenterScore")).doubleValue(), 1e-6); //$NON-NLS-1$
	}

	/**
	 * 
	 * @throws Exception
	 */
	@Test
	public void testGetFragments() throws Exception {
		final double[] expected = { 90.97445, 106.94476,
				// 115.98965,
				117.9854, 133.95592, 143.98846, 144.99625,
				// 146.0041,
				151.94641, 160.96668, 163.00682,
				// 172.99055,
				178.95724, 178.97725,
				// 180.97293,
				// 196.96778,
				208.9678, 236.96245 };

		final double[] fragments = MetFrag.getFragments("C(C(=O)O)OC1=NC(=C(C(=C1Cl)N)Cl)F", 2); //$NON-NLS-1$
		final double epsilon = 1e-3;

		for (double mass : expected) {
			Assert.assertTrue(DoubleStream.of(fragments).anyMatch(x -> x > mass - epsilon && x < mass + epsilon));
		}
	}
}