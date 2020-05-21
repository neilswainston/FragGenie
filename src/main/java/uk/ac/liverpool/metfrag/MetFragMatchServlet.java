/**
 * 
 */
package uk.ac.liverpool.metfrag;

import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import javax.json.Json;
import javax.json.JsonArray;
import javax.json.JsonArrayBuilder;
import javax.json.JsonObject;
import javax.json.JsonObjectBuilder;
import javax.json.JsonReader;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

/**
 * 
 * @author neilswainston
 */
@WebServlet(name = "MetFragMatchServlet", urlPatterns = { "/match" })
public class MetFragMatchServlet extends HttpServlet {

	/**
	 *
	 **/
	private static final long serialVersionUID = 1L;

	@Override
	public void doGet(final HttpServletRequest request, final HttpServletResponse response) throws IOException {
		try {
			run(MetFragTestData.SMILES, MetFragTestData.MZ, MetFragTestData.INTEN, response);
		}
		catch(Exception e) {
			MetFragUtils.handleException(e);
		}
	}
	
	@Override
	public void doPost(final HttpServletRequest request, final HttpServletResponse response) throws IOException {
		try(final JsonReader jsonReader = MetFragUtils.getReader(request)) {
			final JsonObject json = jsonReader.readObject();
			final JsonArray smiles = (JsonArray)json.get("smiles");
			final JsonArray mz = (JsonArray)json.get("mz");
			final JsonArray inten = (JsonArray)json.get("inten");
			
			try {
				run(MetFragUtils.toStringArray(smiles), MetFragUtils.toDoubleArray(mz), MetFragUtils.toIntArray(inten), response);
			}
			catch(Exception e) {
				MetFragUtils.handleException(e);
			}
		}
	}
	
	/**
	 * 
	 * @param smiles
	 * @param mz
	 * @param inten
	 * @param response
	 * @throws Exception
	 */
	private static void run(final String[] smiles, final double[] mz, final int[] inten, final HttpServletResponse response) throws Exception{
		final Collection<Map<String,Object>> results = MetFrag.match(smiles, mz, inten);
		final JsonArray json = toJson(results);
		MetFragUtils.sendJson(json, response);
	}

	/**
	 * 
	 * @param results
	 * @return JsonArray
	 */
	private static JsonArray toJson(final Collection<Map<String,Object>> results) {
		final JsonArrayBuilder arrayBuilder = Json.createArrayBuilder();

		for(Map<String,Object> result : results) {
			final JsonObjectBuilder resBuilder = Json.createObjectBuilder();
			
			for (Map.Entry<String, Object> entry : result.entrySet()) {
				final Object value = entry.getValue();
				
				if(value instanceof String) {
					resBuilder.add(entry.getKey(), (String)value);
				}
				else if(value instanceof Byte) {
					resBuilder.add(entry.getKey(), (Byte)value);
				}
				else if(value instanceof Integer) {
					resBuilder.add(entry.getKey(), (Integer)value);
				}
				else if(value instanceof Double) {
					resBuilder.add(entry.getKey(), (Double)value);
				}
			}
			
			arrayBuilder.add(resBuilder.build());
		}
		
		return arrayBuilder.build();
	}
}