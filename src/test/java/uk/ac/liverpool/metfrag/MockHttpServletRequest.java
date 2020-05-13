/**
 * 
 */
package uk.ac.liverpool.metfrag;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.io.UnsupportedEncodingException;
import java.security.Principal;
import java.util.Collection;
import java.util.Enumeration;
import java.util.Locale;
import java.util.Map;

import javax.json.Json;
import javax.json.JsonArrayBuilder;
import javax.json.JsonObjectBuilder;
import javax.servlet.AsyncContext;
import javax.servlet.DispatcherType;
import javax.servlet.RequestDispatcher;
import javax.servlet.ServletContext;
import javax.servlet.ServletException;
import javax.servlet.ServletInputStream;
import javax.servlet.ServletRequest;
import javax.servlet.ServletResponse;
import javax.servlet.http.Cookie;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.servlet.http.HttpSession;
import javax.servlet.http.HttpUpgradeHandler;
import javax.servlet.http.Part;

/**
 * 
 * @author neilswainston
 */
public class MockHttpServletRequest implements HttpServletRequest {

	/**
	 * 
	 */
	private final String[] smiles = new String[] {
			"C(C(=O)O)OC1=NC(=C(C(=C1Cl)N)Cl)F",
			"C(C(=O)OC1=NC(=C(C(=C1Cl)N)Cl)F)O",
			"C1(=C(C(=NC(=C1Cl)F)C(C(=O)O)O)Cl)N",
			"CC(=O)OOC1=NC(=C(C(=C1Cl)N)Cl)F",
			"C(C(CSCP(=O)([O-])[O-])C(=O)[O-])C(=O)[O-]"
	};
	
	/**
	 * 
	 */
	private final double[] mz = new double[] {
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
	
	private final int[] inten = new int[] {
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
	
	/**
	 * 
	 */
	private final String query;
	
	/**
	 * 
	 */
	public MockHttpServletRequest() {
		final JsonObjectBuilder builder = Json.createObjectBuilder();
		
		final JsonArrayBuilder smilesBuilder = Json.createArrayBuilder();
		final JsonArrayBuilder mzBuilder = Json.createArrayBuilder();
		final JsonArrayBuilder intenBuilder = Json.createArrayBuilder();
		
		for(String value : this.smiles) {
			smilesBuilder.add(value);
		}
		
		for(double value : this.mz) {
			mzBuilder.add(value);
		}
		
		for(int value : this.inten) {
			intenBuilder.add(value);
		}
		
		builder.add("smiles", smilesBuilder.build());
		builder.add("mz", mzBuilder.build());
		builder.add("inten", intenBuilder.build());
		
		this.query = builder.build().toString();
	}

	@Override
	public Object getAttribute(String name) {
		// Empty block
		return null;
	}

	@Override
	public Enumeration<String> getAttributeNames() {
		// Empty block
		return null;
	}

	@Override
	public String getCharacterEncoding() {
		// Empty block
		return null;
	}

	@Override
	public void setCharacterEncoding(String env) throws UnsupportedEncodingException {
		// Empty block

	}

	@Override
	public int getContentLength() {
		// Empty block
		return 0;
	}

	@Override
	public long getContentLengthLong() {
		// Empty block
		return 0;
	}

	@Override
	public String getContentType() {
		// Empty block
		return null;
	}

	@Override
	public ServletInputStream getInputStream() throws IOException {
		// Empty block
		return null;
	}

	@Override
	public String getParameter(String name) {
		// Empty block
		return null;
	}

	@Override
	public Enumeration<String> getParameterNames() {
		// Empty block
		return null;
	}

	@Override
	public String[] getParameterValues(String name) {
		// Empty block
		return null;
	}

	@Override
	public Map<String, String[]> getParameterMap() {
		// Empty block
		return null;
	}

	@Override
	public String getProtocol() {
		// Empty block
		return null;
	}

	@Override
	public String getScheme() {
		// Empty block
		return null;
	}

	@Override
	public String getServerName() {
		// Empty block
		return null;
	}

	@Override
	public int getServerPort() {
		// Empty block
		return 0;
	}

	@Override
	public BufferedReader getReader() throws IOException {
		return new BufferedReader(new StringReader(this.query));
	}

	@Override
	public String getRemoteAddr() {
		// Empty block
		return null;
	}

	@Override
	public String getRemoteHost() {
		// Empty block
		return null;
	}

	@Override
	public void setAttribute(String name, Object o) {
		// Empty block

	}

	@Override
	public void removeAttribute(String name) {
		// Empty block

	}

	@Override
	public Locale getLocale() {
		// Empty block
		return null;
	}

	@Override
	public Enumeration<Locale> getLocales() {
		// Empty block
		return null;
	}

	@Override
	public boolean isSecure() {
		// Empty block
		return false;
	}

	@Override
	public RequestDispatcher getRequestDispatcher(String path) {
		// Empty block
		return null;
	}

	@Override
	public String getRealPath(String path) {
		// Empty block
		return null;
	}

	@Override
	public int getRemotePort() {
		// Empty block
		return 0;
	}

	@Override
	public String getLocalName() {
		// Empty block
		return null;
	}

	@Override
	public String getLocalAddr() {
		// Empty block
		return null;
	}

	@Override
	public int getLocalPort() {
		// Empty block
		return 0;
	}

	@Override
	public ServletContext getServletContext() {
		// Empty block
		return null;
	}

	@Override
	public AsyncContext startAsync() throws IllegalStateException {
		// Empty block
		return null;
	}

	@Override
	public AsyncContext startAsync(ServletRequest servletRequest, ServletResponse servletResponse)
			throws IllegalStateException {
		// Empty block
		return null;
	}

	@Override
	public boolean isAsyncStarted() {
		// Empty block
		return false;
	}

	@Override
	public boolean isAsyncSupported() {
		// Empty block
		return false;
	}

	@Override
	public AsyncContext getAsyncContext() {
		// Empty block
		return null;
	}

	@Override
	public DispatcherType getDispatcherType() {
		// Empty block
		return null;
	}

	@Override
	public String getAuthType() {
		// Empty block
		return null;
	}

	@Override
	public Cookie[] getCookies() {
		// Empty block
		return null;
	}

	@Override
	public long getDateHeader(String name) {
		// Empty block
		return 0;
	}

	@Override
	public String getHeader(String name) {
		// Empty block
		return null;
	}

	@Override
	public Enumeration<String> getHeaders(String name) {
		// Empty block
		return null;
	}

	@Override
	public Enumeration<String> getHeaderNames() {
		// Empty block
		return null;
	}

	@Override
	public int getIntHeader(String name) {
		// Empty block
		return 0;
	}

	@Override
	public String getMethod() {
		// Empty block
		return null;
	}

	@Override
	public String getPathInfo() {
		// Empty block
		return null;
	}

	@Override
	public String getPathTranslated() {
		// Empty block
		return null;
	}

	@Override
	public String getContextPath() {
		// Empty block
		return null;
	}

	@Override
	public String getQueryString() {
		// Empty block
		return null;
	}

	@Override
	public String getRemoteUser() {
		// Empty block
		return null;
	}

	@Override
	public boolean isUserInRole(String role) {
		// Empty block
		return false;
	}

	@Override
	public Principal getUserPrincipal() {
		// Empty block
		return null;
	}

	@Override
	public String getRequestedSessionId() {
		// Empty block
		return null;
	}

	@Override
	public String getRequestURI() {
		// Empty block
		return null;
	}

	@Override
	public StringBuffer getRequestURL() {
		// Empty block
		return null;
	}

	@Override
	public String getServletPath() {
		// Empty block
		return null;
	}

	@Override
	public HttpSession getSession(boolean create) {
		// Empty block
		return null;
	}

	@Override
	public HttpSession getSession() {
		// Empty block
		return null;
	}

	@Override
	public String changeSessionId() {
		// Empty block
		return null;
	}

	@Override
	public boolean isRequestedSessionIdValid() {
		// Empty block
		return false;
	}

	@Override
	public boolean isRequestedSessionIdFromCookie() {
		// Empty block
		return false;
	}

	@Override
	public boolean isRequestedSessionIdFromURL() {
		// Empty block
		return false;
	}

	@Override
	public boolean isRequestedSessionIdFromUrl() {
		// Empty block
		return false;
	}

	@Override
	public boolean authenticate(HttpServletResponse response) throws IOException, ServletException {
		// Empty block
		return false;
	}

	@Override
	public void login(String username, String password) throws ServletException {
		// Empty block

	}

	@Override
	public void logout() throws ServletException {
		// Empty block

	}

	@Override
	public Collection<Part> getParts() throws IOException, ServletException {
		// Empty block
		return null;
	}

	@Override
	public Part getPart(String name) throws IOException, ServletException {
		// Empty block
		return null;
	}

	@Override
	public <T extends HttpUpgradeHandler> T upgrade(Class<T> handlerClass) throws IOException, ServletException {
		// Empty block
		return null;
	}
}