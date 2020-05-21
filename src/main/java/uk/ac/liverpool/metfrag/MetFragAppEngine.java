/**
 * 
 */
package uk.ac.liverpool.metfrag;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.io.Writer;
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
@WebServlet(name = "MetFrag", urlPatterns = { "/match" })
public class MetFragAppEngine extends HttpServlet {

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
			e.printStackTrace();
			
			final Writer writer = new StringWriter();
            e.printStackTrace(new PrintWriter(writer));
			throw new IOException(writer.toString());
		}
	}
	
	@Override
	public void doPost(final HttpServletRequest request, final HttpServletResponse response) throws IOException {
		final StringBuilder builder = new StringBuilder();
		
		try(final BufferedReader reader = request.getReader()) {
			String line;
			
	        while((line = reader.readLine()) != null) {
	        	builder.append(line).append('\n');
	        }
		}

	    final String query = builder.toString();
		
		try(final JsonReader jsonReader = Json.createReader(new StringReader(query))) {
			final JsonObject json = jsonReader.readObject();
			final JsonArray smiles = (JsonArray)json.get("smiles");
			final JsonArray mz = (JsonArray)json.get("mz");
			final JsonArray inten = (JsonArray)json.get("inten");
			
			try {
				run(toStringArray(smiles), toDoubleArray(mz), toIntArray(inten), response);
			}
			catch(Exception e) {
				e.printStackTrace();
				
				final Writer writer = new StringWriter();
	            e.printStackTrace(new PrintWriter(writer));
				throw new IOException(writer.toString());
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
		response.setContentType("application/json");
		response.setCharacterEncoding("UTF-8");
		response.getWriter().print(json.toString());
	}
	
	/**
	 * 
	 * @param jsonArray
	 * @return String[]
	 */
	private static String[] toStringArray(final JsonArray jsonArray) {
		final String[] array = new String[jsonArray.size()];
		
		for(int i=0; i< jsonArray.size(); i++) {
			array[i] = jsonArray.getString(i);
		}
		
		return array;
	}
	
	/**
	 * 
	 * @param jsonArray
	 * @return double[]
	 */
	private static double[] toDoubleArray(final JsonArray jsonArray) {
		final double[] array = new double[jsonArray.size()];
		
		for(int i=0; i< jsonArray.size(); i++) {
			array[i] = jsonArray.getJsonNumber(0).doubleValue();
		}
		
		return array;
	}
	
	/**
	 * 
	 * @param jsonArray
	 * @return int[]
	 */
	private static int[] toIntArray(final JsonArray jsonArray) {
		final int[] array = new int[jsonArray.size()];
		
		for(int i=0; i< jsonArray.size(); i++) {
			array[i] = jsonArray.getInt(0);
		}
		
		return array;
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