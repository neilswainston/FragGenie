/**
 * 
 */
package uk.ac.liverpool.metfrag;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.io.Writer;
import java.util.Collection;
import java.util.Map;
import java.util.Arrays;
import javax.json.Json;
import javax.json.JsonArrayBuilder;
import javax.json.JsonObject;
import javax.json.JsonObjectBuilder;
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
		final String[] smiles = new String[] {
				"C(C(=O)O)OC1=NC(=C(C(=C1Cl)N)Cl)F",
				"C(C(=O)OC1=NC(=C(C(=C1Cl)N)Cl)F)O",
				"C1(=C(C(=NC(=C1Cl)F)C(C(=O)O)O)Cl)N",
				"CC(=O)OOC1=NC(=C(C(=C1Cl)N)Cl)F",
				"C(C(CSCP(=O)([O-])[O-])C(=O)[O-])C(=O)[O-]"
		};
		
		final double[] mz = new double[] {
				90.97445,
				106.94476,
				110.0275,
				115.98965,
				117.9854,
				124.93547,
				124.99015,
				125.99793,
				133.95592,
				143.98846,
				144.99625,
				146.0041,
				151.94641,
				160.96668,
				163.00682,
				172.99055,
				178.95724,
				178.97725,
				180.97293,
				196.96778,
				208.9678,
				236.96245,
				254.97312};
		
		final int[] inten = new int[] {
				681,
				274,
				110,
				95,
				384,
				613,
				146,
				207,
				777,
				478,
				352,
				999,
				962,
				387,
				782,
				17,
				678,
				391,
				999,
				720,
				999,
				999,
				999};
		
		try {
			run(smiles, mz, inten, response);
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
		final String[] smiles = request.getParameterValues("smiles");
		final double[] mz = Arrays.stream(request.getParameterValues("mz")).mapToDouble(Double::parseDouble).toArray();
		final int[] inten = Arrays.stream(request.getParameterValues("inten")).mapToInt(Integer::parseInt).toArray();
		
		try {
			run(smiles, mz, inten, response);
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
	 * @param smiles
	 * @param mz
	 * @param inten
	 * @param response
	 * @throws Exception
	 */
	private static void run(final String[] smiles, final double[] mz, final int[] inten, final HttpServletResponse response) throws Exception{
		final Collection<Map<String,Object>> results = MetFrag.match(smiles, mz, inten);
		final JsonObject json = toJson(results);
		response.setContentType("application/json");
		response.setCharacterEncoding("UTF-8");
		response.getWriter().print(json.toString());
	}

	/**
	 * 
	 * @param results
	 * @return JsonObject
	 */
	private static JsonObject toJson(final Collection<Map<String,Object>> results) {
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
		
		final JsonObjectBuilder objBuilder = Json.createObjectBuilder();
		objBuilder.add("results", arrayBuilder.build());
		return objBuilder.build();
	}
}