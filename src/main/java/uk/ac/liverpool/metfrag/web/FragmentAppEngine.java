/**
 * 
 */
package uk.ac.liverpool.metfrag.web;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.io.Writer;
import java.util.Arrays;
import java.util.List;

import javax.json.Json;
import javax.json.JsonArrayBuilder;
import javax.json.JsonNumber;
import javax.json.JsonObject;
import javax.json.JsonObjectBuilder;
import javax.json.JsonReader;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import uk.ac.liverpool.metfrag.Bond;
import uk.ac.liverpool.metfrag.MetFragFragmenter;
import uk.ac.liverpool.metfrag.MetFragFragmenter.Headers;

/**
 * 
 * @author neilswainston
 */
@WebServlet(name = "liv-metfrag", urlPatterns = { "/fragment" })
public class FragmentAppEngine extends HttpServlet {

	/**
	 *
	 **/
	private static final long serialVersionUID = 1L;

	@Override
	public void doPost(final HttpServletRequest request, final HttpServletResponse response) throws IOException {
		try(final PrintWriter writer = response.getWriter()) {
			final JsonObject query = getQuery(request);
			final String smiles = query.getString("smiles"); //$NON-NLS-1$
			final int maximumTreeDepth = query.getInt("maximumTreeDepth", 3); //$NON-NLS-1$
			final JsonNumber minMassNumber = query.getJsonNumber("minMass"); //$NON-NLS-1$
			final float minMass = minMassNumber != null ? (float)minMassNumber.doubleValue() : 0.0f;
			final List<String> fields = Arrays.asList(query.getString("fields", Headers.METFRAG_MZ.name()).split(",")); //$NON-NLS-1$ //$NON-NLS-2$
			final List<List<Object>> brokenBondsFilter = null;
			final Object[] results = MetFragFragmenter.getFragmentData(smiles, maximumTreeDepth, minMass, fields, brokenBondsFilter);
			final JsonObject json = toJson(results);
			response.setContentType("application/json"); //$NON-NLS-1$
			response.setCharacterEncoding("UTF-8"); //$NON-NLS-1$
			writer.print(json.toString());
		}
		catch(Exception e) {
			e.printStackTrace();
			
			final Writer writer = new StringWriter();
            e.printStackTrace(new PrintWriter(writer));
			throw new IOException(writer.toString());
		}
	}
	
	/**
	 * 
	 * @param request
	 * @return JsonObject
	 * @throws IOException
	 */
	private static JsonObject getQuery(final HttpServletRequest request) throws IOException {
		final StringBuilder builder = new StringBuilder();
		
		try(final BufferedReader reader = request.getReader()) {
			String line;
			
	        while((line = reader.readLine()) != null) {
	        	builder.append(line).append('\n');
	        }
		}

	    final String query = builder.toString();
		
		try(final JsonReader jsonReader = Json.createReader(new StringReader(query))) {
			return jsonReader.readObject();
		}
	}

	/**
	 * 
	 * @param results
	 * @return JsonArray
	 */
	private static JsonObject toJson(final Object[] results) {
		final JsonObjectBuilder resultsBuilder = Json.createObjectBuilder();

		for(int i = 0; i < results.length; i++) {
			final Object value = results[i];
				
			if(value instanceof List) {
				final List<?> lst = (List<?>)value;
				resultsBuilder.add(Headers.values()[i].name(), getArrayBuilder(lst));
			}
		}
		
		return resultsBuilder.build();
	}
	
	/**
	 * 
	 * @param lst
	 * @return JsonArrayBuilder
	 */
	private static JsonArrayBuilder getArrayBuilder(final List<?> lst) {
		final JsonArrayBuilder resultBuilder = Json.createArrayBuilder();
		
		for(final Object value : lst) {
			if(value instanceof String) {
				resultBuilder.add((String)value);
			}
			else if(value instanceof Float) {
				resultBuilder.add(((Float)value).doubleValue());
			}
			else if(value instanceof List) {
				resultBuilder.add(getArrayBuilder((List<?>)value));
			}
			else if(value instanceof Bond) {
				final Bond bond = (Bond)value;
				resultBuilder.add(bond.encode());
			}
		}
		
		return resultBuilder;
	}
}