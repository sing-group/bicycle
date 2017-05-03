package es.cnio.bioinfo.bicycle;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * A PrintStream that redirects messages to a Logger
 * Created by lipido on 3/05/17.
 */
public class StandardStreamsToLoggerRedirector {
	private final Logger logger;
	private final Level level;
	private final MessageFilter messageFilter;
	private final ToLoggerPrintStream errToLoggerPrintStream;
	private final ToLoggerPrintStream outToLoggerPrintStream;
	private static PrintStream originalStdOut;
	private static PrintStream originalStdErr;

	static {
		originalStdOut = System.out;
		originalStdErr = System.err;
	}
	public interface MessageFilter
	{
		public abstract String filter(String msg);
	}
	
	public StandardStreamsToLoggerRedirector(Logger logger, Level level, MessageFilter filter) {
		this.logger = logger;
		this.level = level;
		this.messageFilter = filter;
		this.errToLoggerPrintStream = new ToLoggerPrintStream(this.logger, this.level, this.messageFilter, originalStdErr);
		this.outToLoggerPrintStream = new ToLoggerPrintStream(this.logger, this.level, this.messageFilter, originalStdOut);
	}
	
	public void redirectStreams() {
		System.setErr(errToLoggerPrintStream);
		System.setOut(outToLoggerPrintStream);
	}
	
	public void restoreStreams() {
		this.errToLoggerPrintStream.flush();
		this.outToLoggerPrintStream.flush();
		System.setErr(originalStdErr);
		System.setOut(originalStdOut);
	}

	class ToLoggerPrintStream extends PrintStream {
		private final Logger logger;
		private final Level level;
		private StringBuilder accumulator = new StringBuilder();

		public ToLoggerPrintStream(Logger logger, Level level, MessageFilter filter, PrintStream originalStream) {

			super( new OutputStream() {

				public boolean logging = false;
				private ByteArrayOutputStream baos = new ByteArrayOutputStream();

				@Override
				public void write(int b) throws IOException {
					if (!logging) {
						baos.write(b);
						tryDoLog();
					} else {
						originalStream.write(b);
					}
				}

				@Override
				public void write(byte[] b) throws IOException {
					if (!logging) {
						baos.write(b);
						tryDoLog();
					} else {
						originalStream.write(b);
					}
				}

				@Override
				public void write(byte[] b, int off, int len) throws IOException {
					if (!logging) {
						baos.write(b, off, len);
						tryDoLog();
					} else {
						originalStream.write(b, off, len);
					}
				}

				@Override
				public void close() throws IOException {
					if (!logging) {
						flush();
					} else {
						originalStream.close();
					}
				}

				@Override
				public void flush() throws IOException {
					if (logging) {
						originalStream.flush();
					} else {
						doLog(new String(baos.toByteArray()).length());
					}
				}

				private void tryDoLog() {
					String text = new String(baos.toByteArray());
					if (text.indexOf('\n') != -1) {
						doLog(text.indexOf('\n'));
					}
				}

				private void doLog(int end) {
					String text = new String(baos.toByteArray());
					doLog(text.substring(0, end));
					ByteArrayOutputStream newBaos = new ByteArrayOutputStream();

					if (end < text.length() - 1) {
						//System.out.println(text.length()+ " "+end);

						newBaos.write(this.baos.toByteArray(), end + 1, this.baos.toByteArray().length -
								(end+1));
					}
					this.baos = newBaos;
				}

				protected void doLog(String s) {
					if (s.length() > 0) {
						s = filter.filter(s);
						if (s.length() > 0) {
							this.logging = true;
							logger.log(level, s);
							this.logging = false;
						}
					}
				}

			});
			this.logger = logger;
			this.level = level;
		}

		protected String beforeLogging(String message) {
			return message;
		}


	}

}

