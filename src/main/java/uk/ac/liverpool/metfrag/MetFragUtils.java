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

import javax.json.Json;
import javax.json.JsonArray;
import javax.json.JsonArrayBuilder;
import javax.json.JsonReader;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

/**
 * 
 * @author neilswainston
 */
public abstract class MetFragUtils {

	/**
	 * 
	 * @param e
	 * @throws IOException
	 */
	protected static void handleException(final Exception e) throws IOException{
		e.printStackTrace();
		
		final Writer writer = new StringWriter();
        e.printStackTrace(new PrintWriter(writer));
		throw new IOException(writer.toString());
	}
	
	/**
	 * 
	 * @param request
	 * @return JsonReader
	 * @throws IOException
	 */
	protected static JsonReader getReader(final HttpServletRequest request) throws IOException {
		final StringBuilder builder = new StringBuilder();
		
		try(final BufferedReader reader = request.getReader()) {
			String line;
			
	        while((line = reader.readLine()) != null) {
	        	builder.append(line).append('\n');
	        }
		}

	    final String query = builder.toString();
	    return Json.createReader(new StringReader(query));
	}
	
	/**
	 * 
	 * @param json
	 * @param response
	 * @throws IOException
	 */
	protected static void sendJson(final Object json, final HttpServletResponse response) throws IOException {
		response.setContentType("application/json");
		response.setCharacterEncoding("UTF-8");
		response.getWriter().print(json.toString());
	}
	
	/**
	 * 
	 * @param jsonArray
	 * @return String[]
	 */
	protected static String[] toStringArray(final JsonArray jsonArray) {
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
	protected static double[] toDoubleArray(final JsonArray jsonArray) {
		final double[] array = new double[jsonArray.size()];
		
		for(int i=0; i< jsonArray.size(); i++) {
			array[i] = jsonArray.getJsonNumber(i).doubleValue();
		}
		
		return array;
	}
	
	/**
	 * 
	 * @param jsonArray
	 * @return int[]
	 */
	protected static int[] toIntArray(final JsonArray jsonArray) {
		final int[] array = new int[jsonArray.size()];
		
		for(int i=0; i< jsonArray.size(); i++) {
			array[i] = jsonArray.getInt(i);
		}
		
		return array;
	}
	
	/**
	 * 
	 * @param array
	 * @return JsonArray
	 */
	protected static JsonArray fromDoubleArray(final double[] array) {
		final JsonArrayBuilder arrayBuilder = Json.createArrayBuilder();
		
		for(final double value : array) {
			arrayBuilder.add(value);
		}
		
		return arrayBuilder.build();
	}
}