/**
 * 
 */
package uk.ac.liverpool.metfrag.web;

import java.io.IOException;
import java.io.StringReader;

import javax.json.Json;
import javax.json.JsonArray;
import javax.json.JsonObject;
import javax.json.JsonReader;

import org.junit.Assert;
import org.junit.Test;

import uk.ac.liverpool.metfrag.MetFragFragmenter.Headers;

/**
 * 
 * @author neilswainston
 */
public class FragmentAppEngineTest {

	/**
	 * 
	 * @throws IOException
	 */
	@SuppressWarnings("static-method")
	@Test
	public void testPost() throws IOException {
		final MockHttpServletRequest request = new MockHttpServletRequest();
		final MockHttpServletResponse response = new MockHttpServletResponse();
		new FragmentAppEngine().doPost(request, response);
		verify(response);
	}
	
	/**
	 * 
	 * @param response
	 */
	private static void verify(final MockHttpServletResponse response) {
		Assert.assertEquals("application/json", response.getContentType()); //$NON-NLS-1$
		Assert.assertEquals("UTF-8", response.getCharacterEncoding()); //$NON-NLS-1$

		final float[] expectedMz = new float[] {15.022941f, 16.03077f, 18.010021f, 30.04643f, 32.02568f, 47.04917f};
		final String[] expectedFormulae = new String[] {"CH2", "CH3", "HO", "C2H5", "CH3O", "C2H6O"}; //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$ //$NON-NLS-4$ //$NON-NLS-5$ //$NON-NLS-6$
		final int[][] expectedBrokenBonds = new int[][] {
			{12384, 12416},
			{12384},
			{12384, 12416},
			{12416},
			{12384},
			{}
		};
		
		final String resp = response.getWriterContent().toString();
		
		try(final JsonReader jsonReader = Json.createReader(new StringReader(resp))) {
			final JsonObject results = jsonReader.readObject();
			Assert.assertArrayEquals(expectedMz, toFloatArray(results.getJsonArray(Headers.METFRAG_MZ.name())), 1e-6f);
			Assert.assertArrayEquals(expectedFormulae, toStringArray(results.getJsonArray(Headers.METFRAG_FORMULAE.name())));
			Assert.assertArrayEquals(expectedBrokenBonds, toIntArray(results.getJsonArray(Headers.METFRAG_BROKEN_BONDS.name())));
		}
	}
	
	/**
	 * 
	 * @param jsonArray
	 * @return String[]
	 */
	private static String[] toStringArray(final JsonArray jsonArray) {
		final String[] array = new String[jsonArray.size()];
		
		for(int i = 0; i < jsonArray.size(); i++) {
			array[i] = jsonArray.getString(i);
		}
		
		return array;
	}
	
	/**
	 * 
	 * @param jsonArray
	 * @return double[]
	 */
	private static float[] toFloatArray(final JsonArray jsonArray) {
		final float[] array = new float[jsonArray.size()];
		
		for(int i = 0; i < jsonArray.size(); i++) {
			array[i] = (float)jsonArray.getJsonNumber(i).doubleValue();
		}
		
		return array;
	}
	
	/**
	 * 
	 * @param jsonArray
	 * @return int[][]
	 */
	private static int[][] toIntArray(final JsonArray jsonArray) {
		final int[][] outerArray = new int[jsonArray.size()][];
		
		for(int i = 0; i < jsonArray.size(); i++) {
			final JsonArray innerJsonArray = jsonArray.getJsonArray(i);
			final int[] innerArray = new int[innerJsonArray.size()];
			
			for(int j = 0; j < innerJsonArray.size(); j++) {
				innerArray[j] = innerJsonArray.getInt(j);
			}
			
			outerArray[i] = innerArray;
			
		}
		
		return outerArray;
	}
}