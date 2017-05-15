package es.cnio.bioinfo.bicycle.operations;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.Queue;

public class GtoADuplicatorReader extends BufferedReader {

	private BufferedReader internalReader;

	public GtoADuplicatorReader(BufferedReader reader) {
		super(reader);
		this.internalReader = reader;
	}

	private Queue<String> nextLines = new LinkedList<String>();

	@Override
	public String readLine() throws IOException {

		if (!nextLines.isEmpty()) {
			return nextLines.poll();
		}

		//read four lines
		String header = internalReader.readLine();
		String sequence = internalReader.readLine();
		String header2 = internalReader.readLine();
		String quality = internalReader.readLine();

		if (header == null || sequence == null ||
				header2 == null || quality == null) {

			return null;
		}

		// no modifications
		//nextLines.offer(header);
		nextLines.offer(sequence);
		nextLines.offer(header2);
		nextLines.offer(quality);

		// additional GtoA
		nextLines.offer(header);
		nextLines.offer(header.split("[|][|]")[1].replaceAll("G", "A"));
		nextLines.offer(header2);
		nextLines.offer(quality);


		return header;
	}
}
