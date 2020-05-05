/**
 * 
 */
package uk.ac.liverpool.metfrag;

import java.io.IOException;

import javax.json.Json;
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
		final MetFrag metFrag = new MetFrag();
		final String match = metFrag.match();
		final JsonObject json = toJson(match);
		response.setContentType("application/json");
		response.setCharacterEncoding("UTF-8");
		response.getWriter().print(json.toString());
	}

	/**
	 * 
	 * @param match
	 * @return JsonObject
	 */
	private static JsonObject toJson(final String match) {
		final JsonObjectBuilder jsonBuilder = Json.createObjectBuilder();
		jsonBuilder.add("match", match);
		return jsonBuilder.build();
	}
}