/**
 * 
 */
package uk.ac.liverpool.metfrag;

import java.io.IOException;
import java.io.StringReader;

import javax.json.Json;
import javax.json.JsonArray;
import javax.json.JsonReader;

import org.junit.Assert;
import org.junit.Test;

/**
 * 
 * @author neilswainston
 */
public class MetFragFragmentServletTest {

	/**
	 * 
	 * @throws IOException
	 */
	@Test
	public void testGet() throws IOException {
		final MockHttpServletResponse response = new MockHttpServletResponse();
		new MetFragFragmentServlet().doGet(null, response);
		verify(response);
	}
	
	/**
	 * 
	 * @throws IOException
	 */
	@Test
	public void testPost() throws IOException {
		final MockHttpServletResponse response = new MockHttpServletResponse();
		new MetFragFragmentServlet().doPost(new MockMetFragFragmentServletRequest(), response);
		verify(response);
	}
	
	/**
	 * 
	 * @param response
	 */
	private static void verify(final MockHttpServletResponse response) {
		Assert.assertEquals("application/json", response.getContentType()); //$NON-NLS-1$
		Assert.assertEquals("UTF-8", response.getCharacterEncoding()); //$NON-NLS-1$

		final String resp = response.getWriterContent().toString();

		try(final JsonReader jsonReader = Json.createReader(new StringReader(resp))) {
			final JsonArray results = jsonReader.readArray();
			Assert.assertEquals(5, results.size());
		}
	}
}