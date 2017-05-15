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

public class SampleBisulfitation {
	private static final Logger logger = Logger.getLogger(SampleBisulfitation.class.getSimpleName());
	private Sample sample;

	public SampleBisulfitation(Sample sample) {
		this.sample = sample;
	}

	public void computeSampleBisulfitation(boolean removeReadsWithUnconvertedBarCodes) throws IOException {
		for (File f : this.sample.getReadsFiles()) {
			logger.info("Performing CtoT in-silico bisulfitation for " + f + "...... ");

			File outputFile = getBisulfitedFile(f);
			File unconvertedReadsFile = getUnconvertedReads(f);

			if (outputFile.exists()) {
				logger.info("Removing existent file: " + outputFile);
				outputFile.delete();
			}

			BufferedWriter wr = new BufferedWriter(new FileWriter(outputFile));
			BufferedWriter unConvertedReads = new BufferedWriter(new FileWriter(unconvertedReadsFile));
			BufferedReader br = new BufferedReader(new FileReader(f));
			String lineModified;
			String thisLine = null;
			while ((thisLine = br.readLine()) != null) {
				if (thisLine.startsWith("@")) {
					// primero adjunto en la cabecera fasta la read original

					// leo la secuencia de la read
					String originalSequence = br.readLine();
					boolean useRead = true;

					if (removeReadsWithUnconvertedBarCodes) {
						// obtengo el barcode del nombre del archivo, ejemplo: ES_LIF_s_8_TGtATT-sequence.txt
						String barcode = new String((f.getName().split("-"))[0]);
						String aux[] = barcode.split("_");
						barcode = aux[aux.length - 1];
						// localizo la posicion de la 't' en el barcode, posicion que vendria de una c no metilada
						int barcodePosition = -1; //inicializo
						barcodePosition = barcode.indexOf("t");
						String[] tokens = thisLine.split("[#]");
						if (tokens.length == 2) {
							String thisReadBarcode = (thisLine.split("[#]"))[1];

							if (thisReadBarcode.charAt(barcodePosition) != 't' && thisReadBarcode.charAt
									(barcodePosition) != 'T')
								useRead = false;


							if (barcodePosition == -1) {
								String error = "[ERROR]: no barcodes defined with --barcodes";
								logger.severe(error);
								throw new RuntimeException(error);
							}
						}
					}//if(removeReadsWithUnconvertedBarCodes)

					// en caso de que quiera comprobar si la read está mal bisulfitada y eliminarla, entonces voy por
					// aquí
					 /*	if(removeUnconvertedReads){
		        	 		if(originalSequence.matches(".*[Cc][^Gg].*[Cc][^Gg].*[Cc][^Gg].*[Cc][^Gg].*")){// apunto
		        	 		esta read en el fichero de las que estan mal convertidas por el bisulfito
		        	 			unConvertedReads.write(thisLine);
		        	 			unConvertedReads.newLine();
		        	 			unConvertedReads.write(originalSequence);
		        	 			unConvertedReads.newLine();
		        	 			// por ultimo copio las dos ultimas lineas del fastQ
		        	 			unConvertedReads.write(br.readLine());
		        	 			unConvertedReads.newLine();
		        	 			unConvertedReads.write(br.readLine());
		        	 			unConvertedReads.newLine();
		        	 			
		        	 			useRead=false;
		        	 		}
		        	 		
		        	 	}*/

					if (useRead) {//que no quiere comprobar si la read está mal bisulfitada y pasa de eliminarla,
						// vamos por aquí
						//genero la nueva cabecera fasta
						lineModified = (new StringBuilder(thisLine.replace(' ', '_'))).append("||").append
								(originalSequence).toString();
						// la escribo
						wr.write(lineModified);
						wr.newLine();


						// ahora bisulfito la secuencia
						lineModified = originalSequence.replace('C', 'T');
						lineModified = lineModified.replace('c', 't');
						// escribo la nueva secuencia
						wr.write(lineModified);
						wr.newLine();
						// por ultimo copio las dos ultimas lineas del fastQ
						wr.write(br.readLine());
						wr.newLine();
						wr.write(br.readLine());
						wr.newLine();
					}
				}

			} // end while
			wr.close();
			br.close();

		}
		logger.info("[OK]");

	}

	public File getBisulfitedFile(File file) {
		if (!this.sample.getReadsFiles().contains(file)) {
			throw new IllegalArgumentException("file " + file + " is not in this sample");
		}
		return new File(this.sample.getProject().getWorkingDirectory() + File.separator + "bisulfited_CT_" +
				getPlainFileName(file));
	}

	private String getPlainFileName(File file) {
		if (this.sample.getReadsFiles().size() == 1 && this.sample.getName().equals(this.sample.getReadsFiles().get(0)
				.getName())) {
			//1-fastq sample
			return file.getName();
		} else {
			return this.sample.getName() + "_" + file.getName();
		}
	}

	public File getUnconvertedReads(File file) {
		return new File(this.sample.getProject().getWorkingDirectory() + File.separator + "unconvertedReads_" +
				getPlainFileName(file));

	}
}
