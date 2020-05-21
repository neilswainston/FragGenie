/**
 * 
 */
package uk.ac.liverpool.metfrag;

import javax.json.Json;
import javax.json.JsonArrayBuilder;
import javax.json.JsonObjectBuilder;

/**
 * 
 * @author neilswainston
 */
public class MockMetFragFragmentServletRequest extends MockHttpServletRequest {
	
	/**
	 * 
	 */
	public MockMetFragFragmentServletRequest() {
		super(getQuery());
	}
	
	/**
	 * 
	 * @return String
	 */
	private static String getQuery() {
		final JsonObjectBuilder builder = Json.createObjectBuilder();
		
		final JsonArrayBuilder smilesBuilder = Json.createArrayBuilder();
		
		for(String value : MetFragTestData.SMILES) {
			smilesBuilder.add(value);
		}
		
		builder.add("smiles", smilesBuilder.build());
		builder.add("maximumTreeDepth", 2);
		
		return builder.build().toString();
	}
}