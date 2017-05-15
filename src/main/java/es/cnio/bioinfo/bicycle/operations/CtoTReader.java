package es.cnio.bioinfo.bicycle.operations;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.Queue;

import es.cnio.bioinfo.bicycle.Sample;

public class CtoTReader extends ReadsReader {


	private boolean skipUnconverted = false;

	public CtoTReader(Sample sample, BufferedReader reader) {
		super(sample, reader);
	}

	public CtoTReader(Sample sample, BufferedReader reader, boolean skipUnconverted) {
		super(sample, reader);
		this.skipUnconverted = skipUnconverted;
	}

	private Queue<String> nextLines = new LinkedList<String>();

	@Override
	public String readLine() throws IOException {

		if (!nextLines.isEmpty()) {
			return nextLines.poll();
		}

		//read four lines
		String header = null;
		String sequence = null;
		String header2 = null;
		String quality = null;

		do {
			header = internalReader.readLine();

			sequence = internalReader.readLine();
			header2 = internalReader.readLine();
			quality = internalReader.readLine();

		} while (header != null && shouldSkip(header));

		if (header == null || sequence == null ||
				header2 == null || quality == null) {

			return null;
		}

		//nextLines.offer(header);
		nextLines.offer(sequence.replaceAll("C", "T"));
		nextLines.offer(header2);
		nextLines.offer(quality);

		return header.replaceAll(" ", "_") + "||" + sequence;
	}


	private boolean shouldSkip(String header) {
		return this.skipUnconverted && this.hasUnconvertedBarcode(header);
	}
}