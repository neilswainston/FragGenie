/**
 * 
 */
package uk.ac.liverpool.metfrag.web;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.io.UnsupportedEncodingException;
import java.security.Principal;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Enumeration;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import javax.json.Json;
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

import uk.ac.liverpool.metfrag.MetFragFragmenter.Headers;

/**
 * 
 * @author neilswainston
 */
public class MockHttpServletRequest implements HttpServletRequest {
	
	/**
	 * 
	 */
	private final String query;
	
	/**
	 * 
	 */
	public MockHttpServletRequest() {
		final JsonObjectBuilder builder = Json.createObjectBuilder();
		
		final List<String> fields = new ArrayList<>();
		
		for(final Headers header : Headers.values()) {
			fields.add(header.name());
		}
		
		builder.add("smiles", "CCO"); //$NON-NLS-1$ //$NON-NLS-2$
		builder.add("maximumTreeDepth", 3); //$NON-NLS-1$
		builder.add("minMass", 0.0); //$NON-NLS-1$
		builder.add("fields", String.join(",", fields)); //$NON-NLS-1$ //$NON-NLS-2$
		
		this.query = builder.build().toString();
	}

	@Override
	public Object getAttribute(final String name) {
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
	public void setCharacterEncoding(final String env) throws UnsupportedEncodingException {
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
	public String getParameter(final String name) {
		// Empty block
		return null;
	}

	@Override
	public Enumeration<String> getParameterNames() {
		// Empty block
		return null;
	}

	@Override
	public String[] getParameterValues(final String name) {
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
	public void setAttribute(final String name, final Object o) {
		// Empty block

	}

	@Override
	public void removeAttribute(final String name) {
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
	public RequestDispatcher getRequestDispatcher(final String path) {
		// Empty block
		return null;
	}

	@SuppressWarnings("deprecation")
	@Deprecated
	@Override
	public String getRealPath(final String path) {
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
	public AsyncContext startAsync(final ServletRequest servletRequest, final ServletResponse servletResponse)
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
	public long getDateHeader(final String name) {
		// Empty block
		return 0;
	}

	@Override
	public String getHeader(final String name) {
		// Empty block
		return null;
	}

	@Override
	public Enumeration<String> getHeaders(final String name) {
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
	public boolean isUserInRole(final String role) {
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
	public HttpSession getSession(final boolean create) {
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

	@SuppressWarnings("deprecation")
	@Deprecated
	@Override
	public boolean isRequestedSessionIdFromUrl() {
		// Empty block
		return false;
	}

	@Override
	public boolean authenticate(final HttpServletResponse response) throws IOException, ServletException {
		// Empty block
		return false;
	}

	@Override
	public void login(final String username, final String password) throws ServletException {
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
	public Part getPart(final String name) throws IOException, ServletException {
		// Empty block
		return null;
	}

	@Override
	public <T extends HttpUpgradeHandler> T upgrade(final Class<T> handlerClass) throws IOException, ServletException {
		// Empty block
		return null;
	}
}