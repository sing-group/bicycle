/*

Copyright 2012 Daniel Gonzalez Peña, Osvaldo Graña


This file is part of the bicycle Project. 

bicycle Project is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

bicycle Project is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser Public License for more details.

You should have received a copy of the GNU Lesser Public License
along with bicycle Project.  If not, see <http://www.gnu.org/licenses/>.
*/

package es.cnio.bioinfo.bicycle.operations;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Logger;

import es.cnio.bioinfo.bicycle.Sample;

public class BarcodeErrorComputation {
	private static final Logger logger = Logger.getLogger(BarcodeErrorComputation.class.getSimpleName());
	private Sample sample;

	public BarcodeErrorComputation(Sample sample) {
		this.sample = sample;
	}

	public double computeErrorFromBarcodes() throws IOException {
		int correctBarcodeConversion = 0;
		int failedBarcodeConversion = 0;
		double errorInSample = -1d; // inicializo a -1 (flag para indicar que NO ha podido calcularse el error)

		if (getErrorFile().exists()) {
			errorInSample = this.readErrorRate();
		} else {
			for (File inputFile : this.sample.getReadsFiles()) {

				BufferedReader br = new BufferedReader(new FileReader(inputFile));
				String line = null;
				while ((line = br.readLine()) != null) {
					if (line.startsWith("@")) {
						// obtengo el barcode del nombre del archivo, ejemplo: ES_LIF_s_8_TGtATT-sequence.txt
						String barcode = new String((inputFile.getName().split("-"))[0]);
						String aux[] = barcode.split("_");
						barcode = aux[aux.length - 1];
						// localizo la posicion de la 't' en el barcode, posicion que vendria de una c no metilada
						int barcodePosition = barcode.indexOf("t");
						String thisReadBarcode = (line.split("[#]"))[1];

						//System.out.println(barcode[i]+"--------------"+thisReadBarcode
						// +"----------------------"+thisLine);


						// AL LORO CON EL CALCULO DE ERROR RATES:
						// error Rate=n. errores/(n. errores+n. aciertos)
						// Se cuenta como error de bisulfito la presencia de 'C/c' en la posicion a mirar en el
						// barcode.
						// Se cuenta como error de secuenciación la presencia de 'A' o 'G' en la posicion a mirar del
						// barcode.
						// Se cuenta como read correctamente convertida cuando hay una 'T/t' en la posicion a mirar
						// del barcode.
						// SEGUN ESTO: el n. de errores=n. reads con error en bisulfito+n. reads con error en
						// secuenciación
						// *** de esta manera lo calculó lister, ver supplementary info página 23
						// *** en los datos que le di a Orlando de error de bisulfito (obtenidos con el script de
						// Perl) hay ligeras
						// diferencias porque consideré como error sólo la presencia de 'C' en la posición a mirar del
						// barcode, con lo que
						// si un barcode estaba mal secuenciado en esa posición (A ó G) se contó como acierto (como si
						// tuviese T) y
						// eso no es correcto
						if (thisReadBarcode.charAt(barcodePosition) == 't' || thisReadBarcode.charAt(barcodePosition)
								== 'T')
							correctBarcodeConversion++;
						else {
							failedBarcodeConversion++;

						}
					}//if(thisLine.startsWith("@"))
				} // end while	

				br.close();
			}
			errorInSample = ((double) failedBarcodeConversion) / ((double) (correctBarcodeConversion +
					failedBarcodeConversion));
		}
		this.writeErrorRate(errorInSample);
		logger.info("Error rate calculated in bisulfite conversion for " + getErrorFile() + " = " + errorInSample +
				"]");

		return errorInSample;
	}


	private double readErrorRate() throws IOException {
		final File errorFile = new File(this.sample.getProject().getOutputDirectory() + File.separator + this.sample
				.getName() + ".errorRate");

		BufferedReader reader = new BufferedReader(new FileReader(errorFile));
		String line = null;
		if ((line = reader.readLine()) != null) {
			if (line.split("=").length == 2)
				return Double.parseDouble(line.split("=")[1].trim());
		}
		throw new IllegalArgumentException("File " + errorFile + " does not contain a valid error rate");

	}

	private void writeErrorRate(double errorRate) throws IOException {
		final File errorFile = getErrorFile();

		BufferedWriter errorFileWriter = new BufferedWriter(new FileWriter(errorFile));
		errorFileWriter.write("Error rate = " + errorRate);
		errorFileWriter.flush();
		errorFileWriter.close();
	}

	private File getErrorFile() {
		return new File(this.sample.getProject().getOutputDirectory() + File.separator + this.sample.getName() + "" +
				".errorRate");

	}
}
