package es.cnio.bioinfo.bicycle.operations;

import java.io.BufferedReader;
import java.io.File;

import es.cnio.bioinfo.bicycle.Sample;

public class ReadsReader extends BufferedReader {
	protected BufferedReader internalReader;
	private int barcodePosition = -1;
	private Sample sample;

	public ReadsReader(Sample sample, BufferedReader baseReader) {
		super(baseReader);
		this.internalReader = baseReader;
		this.sample = sample;

	}

	private int getBarcodePosition() {
		if (barcodePosition == -1) {
			File f = sample.getReadsFiles().get(0);
			// obtengo el barcode del nombre del archivo, ejemplo: ES_LIF_s_8_TGtATT-sequence.txt
			String barcode = new String((f.getName().split("-"))[0]);
			String aux[] = barcode.split("_");
			barcode = aux[aux.length - 1];

			this.barcodePosition = -1; //inicializo
			barcodePosition = barcode.indexOf("t");
			if (barcodePosition == -1) {
				throw new IllegalArgumentException("no 't' found on barcode. barcode " + barcode);
			}
		}

		return barcodePosition;
	}

	public boolean hasUnconvertedBarcode(String readHeader) {
		String[] tokens = readHeader.split("[#]");
		if (tokens.length == 2) {
			String thisReadBarcode = tokens[1];

			if (thisReadBarcode.charAt(getBarcodePosition()) != 't' && thisReadBarcode.charAt(getBarcodePosition()) !=
					'T') {
				return true;
			}
		}
		return false;
	}
}
