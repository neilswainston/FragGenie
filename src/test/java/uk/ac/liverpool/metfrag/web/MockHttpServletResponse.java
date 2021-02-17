/**
 * 
 */
package uk.ac.liverpool.metfrag.web;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Collection;
import java.util.Locale;

import javax.servlet.ServletOutputStream;
import javax.servlet.http.Cookie;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletResponse;

/**
 * This mock class is created to enable basic unit testing of the
 * {@link HelloAppEngine} class. Only methods used in the unit test have a
 * non-trivial implementation.
 * 
 * Feel free to change this class or replace it using other ways for testing
 * {@link HttpServlet}s, e.g. Spring MVC Test or Mockito to suit your needs.
 */
class MockHttpServletResponse implements HttpServletResponse {

	/**
	 * 
	 */
	private String contType;

	/**
	 * 
	 */
	private String enc;

	/**
	 * 
	 */
	private StringWriter writerContent = new StringWriter();

	/**
	 * 
	 */
	private PrintWriter writer = new PrintWriter(this.writerContent);

	@Override
	public void setContentType(final String contentType) {
		this.contType = contentType;
	}

	@Override
	public String getContentType() {
		return this.contType;
	}

	@Override
	public PrintWriter getWriter() throws IOException {
		return this.writer;
	}

	public StringWriter getWriterContent() {
		return this.writerContent;
	}

	@Override
	public void flushBuffer() throws IOException {
		// Empty block
	}

	@Override
	public int getBufferSize() {
		return 0;
	}

	@Override
	public String getCharacterEncoding() {
		return this.enc;
	}

	@Override
	public Locale getLocale() {
		return null;
	}

	@Override
	public ServletOutputStream getOutputStream() throws IOException {
		return null;
	}

	@Override
	public boolean isCommitted() {
		return false;
	}

	@Override
	public void reset() {
		// Empty block
	}

	@Override
	public void resetBuffer() {
		// Empty block
	}

	@Override
	public void setBufferSize(final int size) {
		// Empty block
	}

	@Override
	public void setCharacterEncoding(final String encoding) {
		this.enc = encoding;
	}

	@Override
	public void setContentLength(final int length) {
		// Empty block
	}

	@Override
	public void setLocale(final Locale locale) {
		// Empty block
	}

	@Override
	public void addCookie(final Cookie cookie) {
		// Empty block
	}

	@Override
	public void addDateHeader(final String name, final long date) {
		// Empty block
	}

	@Override
	public void addHeader(final String name, final String value) {
		// Empty block
	}

	@Override
	public void addIntHeader(final String name, final int value) {
		// Empty block
	}

	@Override
	public boolean containsHeader(final String name) {
		return false;
	}

	@Override
	public String encodeRedirectURL(final String url) {
		return null;
	}

	@Override
	public String encodeRedirectUrl(final String url) {
		return null;
	}

	@Override
	public String encodeURL(final String url) {
		return null;
	}

	@Override
	public String encodeUrl(final String url) {
		return null;
	}

	@Override
	public void sendError(final int statusCode) throws IOException {
		// Empty block
	}

	@Override
	public void sendError(final int statusCode, final String message) throws IOException {
		// Empty block
	}

	@Override
	public void sendRedirect(final String url) throws IOException {
		// Empty block
	}

	@Override
	public void setDateHeader(final String name, final long date) {
		// Empty block
	}

	@Override
	public void setHeader(final String name, final String value) {
		// Empty block
	}

	@Override
	public void setIntHeader(final String name, final int value) {
		// Empty block
	}

	@Override
	public void setStatus(final int statusCode) {
		// Empty block
	}

	@Override
	public void setStatus(final int statusCode, final String message) {
		// Empty block
	}

	// Servlet API 3.0 and 3.1 methods
	@Override
	public void setContentLengthLong(final long length) {
		// Empty block
	}

	@Override
	public int getStatus() {
		return 0;
	}

	@Override
	public String getHeader(final String name) {
		return null;
	}

	@Override
	public Collection<String> getHeaders(final String name) {
		return null;
	}

	@Override
	public Collection<String> getHeaderNames() {
		return null;
	}
}