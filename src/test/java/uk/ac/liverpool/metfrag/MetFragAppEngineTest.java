package uk.ac.liverpool.metfrag;

import java.io.IOException;

import org.junit.Assert;
import org.junit.Test;

/**
 * 
 * @author neilswainston
 */
public class MetFragAppEngineTest {

	/**
	 * 
	 * @throws IOException
	 */
	@Test
	public void test() throws IOException {
		final MockHttpServletResponse response = new MockHttpServletResponse();
		new MetFragAppEngine().doGet(null, response);
		Assert.assertEquals("text/plain", response.getContentType());
		Assert.assertEquals("UTF-8", response.getCharacterEncoding());
		Assert.assertEquals("MetFrag match", response.getWriterContent().toString());
	}
}