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
public class MockMetFragMatchServletRequest extends MockHttpServletRequest {
	
	/**
	 * 
	 */
	public MockMetFragMatchServletRequest() {
		super(getQuery());
	}
	
	/**
	 * 
	 * @return String
	 */
	private static String getQuery() {
		final JsonObjectBuilder builder = Json.createObjectBuilder();
		
		final JsonArrayBuilder smilesBuilder = Json.createArrayBuilder();
		final JsonArrayBuilder mzBuilder = Json.createArrayBuilder();
		final JsonArrayBuilder intenBuilder = Json.createArrayBuilder();
		
		for(String value : MetFragTestData.SMILES) {
			smilesBuilder.add(value);
		}
		
		for(double value : MetFragTestData.MZ) {
			mzBuilder.add(value);
		}
		
		for(int value : MetFragTestData.INTEN) {
			intenBuilder.add(value);
		}
		
		builder.add("smiles", smilesBuilder.build());
		builder.add("mz", mzBuilder.build());
		builder.add("inten", intenBuilder.build());
		
		return builder.build().toString();
	}
}