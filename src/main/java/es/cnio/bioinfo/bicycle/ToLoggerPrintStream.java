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
public class ToLoggerPrintStream extends PrintStream {
	private final Logger logger;
	private final Level level;
	private StringBuilder accumulator = new StringBuilder();

	public interface MessageFilter
	{
		public abstract String filter(String msg);
	}

	public ToLoggerPrintStream(Logger logger, Level level, MessageFilter filter) {

		super( new OutputStream() {

			private ByteArrayOutputStream baos = new ByteArrayOutputStream();

			@Override
			public void write(int b) throws IOException {

				baos.write(b);
				tryDoLog();
			}

			@Override
			public void write(byte[] b) throws IOException {

				baos.write(b);
				tryDoLog();
			}

			@Override
			public void write(byte[] b, int off, int len) throws IOException {

				baos.write(b, off, len);
				tryDoLog();
			}

			@Override
			public void close() throws IOException {
				doLog(new String(baos.toByteArray()).length());
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
						logger.log(level, s);
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
/*
	@Override
	public void print(double d) {
		accumulator.append(d);
	}

	@Override
	public void print(float f) {
		accumulator.append(f);
	}

	@Override
	public void print(int i) {
		accumulator.append(i);
	}

	@Override
	public void print(long l) {
		accumulator.append(l);
	}

	@Override
	public void println(boolean x) {
		accumulator.append(x);
		this.doLog();
	}

	@Override
	public void println(char x) {
		accumulator.append(x);
		this.doLog();
	}

	@Override
	public void println(int x) {
		accumulator.append(x);
		this.doLog();
	}

	@Override
	public void println(long x) {
		accumulator.append(x);
	}

	@Override
	public void println(float x) {
		accumulator.append(x);
	}

	@Override
	public void println(double x) {
		accumulator.append(x);
	}

	@Override
	public void write(byte[] b) throws IOException {
		accumulator.append(b);
	}

	@Override
	public void write(int b) {
		accumulator.append(b);
	}

	@Override
	public void print(boolean b) {
		accumulator.append(b);
	}

	@Override
	public void print(String s) {
		accumulator.append(s);
	}

	@Override
	public PrintStream printf(String format, Object... args) {
		return super.printf(format, args);
	}

	@Override
	public void println(String s) {
		accumulator.append(s);
		this.doLog();
	}

	@Override
	public void println() {
		this.doLog();
	}

	@Override
	public void print(Object obj) {
		accumulator.append(obj.toString());
	}

	@Override
	public void print(char c) {
		accumulator.append(c);
	}

	@Override
	public void print(char[] c) {
		accumulator.append(c);
	}

	@Override
	public void println(char[] c) {
		accumulator.append(c);
		this.doLog();
	}

	@Override
	public void println(Object o) {
		accumulator.append(o);
		this.doLog();
	}

	@Override
	public PrintStream append(char c) {
		accumulator.append(c);
		return this;
	}

	@Override
	public PrintStream append(CharSequence csq) {
		accumulator.append(csq);
		return this;
	}

	@Override
	public PrintStream append(CharSequence csq, int start, int end) {
		accumulator.append(csq, start, end);
		return this;
	}
*/



}
