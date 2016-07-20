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

package es.cnio.bioinfo.bicycle;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.management.ManagementFactory;
import java.lang.management.RuntimeMXBean;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.SamLocusIterator;
import net.sf.picard.util.SamLocusIterator.LocusInfo;
import net.sf.picard.util.SamLocusIterator.RecordAndOffset;
import net.sf.samtools.SAMFileReader;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;

import es.cnio.bioinfo.bicycle.gatk.Strand;
import es.cnio.bioinfo.pileline.core.IntervalsIndex;
import es.cnio.bioinfo.pileline.core.IntervalsIndexFactory;
import es.cnio.bioinfo.pileline.io.LineFilter;
import es.cnio.bioinfo.pileline.refgenomeindex.GenomeIndex;
import es.cnio.bioinfo.pileline.refgenomeindex.GenomeIndexBuilder;
import es.cnio.bioinfo.pileline.refgenomeindex.GenomeIndexFactory;

/**
A pipeline for methylation analysis
@author Copyright 2010<br>Osvaldo Gra&ntilde;a, Bioinformatics Unit (CNIO)<br>Daniel Gonz&aacute;lez-Pe&ntilde;a, Higher Technical School of Computer Engineering, University of Vigo, Ourense, Spain<br><br>Distributed under the terms of the GNU General Public License. (Free Software Foundation)
@version 1.0, may 2010<br>
@deprecated This class is no longer supported
*/
public class MethylSeqOperations {

	protected static double FDR;
	private static boolean DISCARD_NON_CT=true; //do not count bases at reference cytosines which are not C or T
	protected static int DEPTH_FILTER=5; // required read depth for a cytosine to be considered in the analysis
	protected static boolean REMOVE_CLONAL; //remove clonal reads (possible PCR duplicates). keep only the best mapped read for each possible 5' position
	private static boolean CORRECT_NONCG=true; //corrects those nonCG contexts to CG where the reads have G's in the +1 nucleotide (but the reference not)
	
	
	private static double[] default_error_rates = {-1d,-1d,-1d}; //-1 means no default, compute it experimentally.
	/**
		 * updates the content of the config file of the project
		 * @param keyword
		 * @param projectFolder
		 * @param data
		 */
		public static void updateConfigFile(String keyword,StringBuilder projectFolder,String [] data){
			// almaceno el nombre de las referencias bisulfitadas en el config.txt
			StringBuilder contenido=new StringBuilder(keyword).append(":");
			
			for(int i=0;i<data.length;i++){
				if(i<data.length-1) contenido.append(data[i]).append(",");
				else contenido.append(data[i]);
			}
			
			String oldConfig=Tools.readFile(new StringBuilder(projectFolder).append("config.txt").toString());
			String[] lines=oldConfig.split("[\n]");
			//actualizo el contenido del fichero config
			StringBuilder newConfig=new StringBuilder("");
			boolean updated=false;
			for(int i=0;i<lines.length;i++){
				if(lines[i].indexOf(keyword)==-1) newConfig.append(lines).append("\n");
				else{
					newConfig.append(contenido).append("\n");
					updated=true;
				}
			}
			//si ha llegado a final de fichero y la linea no existia previamente, entonces tiene que añadirla
			if(!updated) newConfig.append(contenido).append("\n");
			Tools.saveAppending(projectFolder.toString(),"config.txt",newConfig.toString());				
				
		}
		
		/**
		 * Returns a String array with the files
		 * @param files
		 * @return
		 */
		public static String[] getFiles(StringBuilder files){
			//files recibe una cadena con el fichero o ficheros (en caso de varios ficheros irian separados por comas)
			//y los devuelve en un array
			if(files.indexOf(",")!=-1) return(files.toString().split("[,]"));
			else return(new String[]{files.toString()});
		}

		/**
		 * Performs the CtoT bisulfitation of the reference sequence
		 * @param files the list of reference files to be bisulfited, in fasta format
		 * @param path path to the reference files		
		 * @param workingDirectory where the files generated during the execution are stored
		 * @throws FileNotFoundException
		 * @return
		 */
		public static String[] performReferenceInsilicoBisulfitation_CtoT(String[] files, StringBuilder path,StringBuilder workingDirectory){	
			
			String thisLine;
			String[] bisulfited_CT_references=new String[files.length];
			for(int i=0;i<bisulfited_CT_references.length;i++) bisulfited_CT_references[i]="bisulfited_CT_"+files[i];

			for(int i=0;i<files.length;i++){
				System.out.print("Performing CtoT in-silico bisulfitation for "+files[i]+"...... ");				
				// abro buffer de escritura
				BufferedWriter wr = null;
				try {
					String inputFile=new StringBuilder(path).append(files[i]).toString();
					String outputFile=new StringBuilder(workingDirectory).append(bisulfited_CT_references[i]).toString();
					// si existe ya el fichero lo borro
					if(Tools.doesThisFileExist(outputFile)) Tools.deleteFile(outputFile);
					
					wr = new BufferedWriter(new FileWriter(outputFile));
					
					try{
						// abro buffer de lectura
						BufferedReader br = new BufferedReader(new FileReader(inputFile));
						String lineModified;
						
						while ((thisLine = br.readLine()) != null){
					         if(!thisLine.startsWith(">")){
					        	 	// bisulfito la secuencia
					        	 	lineModified=thisLine.replace('C', 'T');
					        	 	lineModified=lineModified.replace('c', 't');
					        	 	wr.write(lineModified);
					        	 	wr.newLine();
								
					         }else{
					        	 	String newFastaHeader=new StringBuilder(thisLine.replace(' ','_')).toString();					        	 	
					        	 	wr.write(newFastaHeader);
					        	 	wr.newLine();
					         }
					         
					         
					    } // end while
						
						// cierro el buffer de lectura
						br.close();
						
					}
					catch (IOException e) {
						       System.err.println("Error: " + e);
					}
					
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}finally {
					//Close the BufferedWriter
					try {
						if (wr != null) {
							wr.flush();
							wr.close();
						}
					
					} catch (IOException ex) {
						ex.printStackTrace();
					}
				}	
				
				System.out.println("[OK]");
			}//for

			return(bisulfited_CT_references);
		
		}
		
		/**
		 * Performs the GtoA bisulfitation of the reference sequence
		 * @param files the list of reference files to be bisulfited, in fasta format
		 * @param path path to the reference files	
		 * @param workingDirectory where the files generated during the execution are stored
		 * @throws FileNotFoundException 
		 * @return
		 */
		public static String[] performReferenceInsilicoBisulfitation_GtoA(String[] files, StringBuilder path,StringBuilder workingDirectory){	
			
			String thisLine;
			String[] bisulfited_GA_references=new String[files.length];
			for(int i=0;i<bisulfited_GA_references.length;i++) bisulfited_GA_references[i]="bisulfited_GA_"+files[i];
			
			for(int i=0;i<files.length;i++){
			
				System.out.print("Performing GtoA in-silico bisulfitation for "+files[i]+"...... ");
				
				// abro buffer de escritura
				BufferedWriter wr = null;
				try {
					String inputFile=new StringBuilder(path).append(files[i]).toString();
					String outputFile=new StringBuilder(workingDirectory).append(bisulfited_GA_references[i]).toString();
					// si existe ya el fichero lo borro
					if(Tools.doesThisFileExist(outputFile)) Tools.deleteFile(outputFile);
					
					wr = new BufferedWriter(new FileWriter(outputFile));
					
					try{
						// abro buffer de lectura
						BufferedReader br = new BufferedReader(new FileReader(inputFile));
						String lineModified;
						
						while ((thisLine = br.readLine()) != null){
					         if(!thisLine.startsWith(">")){
					        	 	// bisulfito la secuencia
					        	 	lineModified=thisLine.replace('G', 'A');
					        	 	lineModified=lineModified.replace('g', 'a');
					        	 	wr.write(lineModified);
					        	 	wr.newLine();
								
					         }else{
					        	 	String newFastaHeader=new StringBuilder(thisLine.replace(' ','_')).toString();
					        	 	
					        	 	wr.write(newFastaHeader);
					        	 	wr.newLine();
					         }
					         
					         
					    } // end while
						
						// cierro el buffer de lectura
						br.close();
						
					}
					catch (IOException e) {
						       System.err.println("Error: " + e);
					}
					
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}finally {
					//Close the BufferedWriter
					try {
						if (wr != null) {
							wr.flush();
							wr.close();
						}
					
					} catch (IOException ex) {
						ex.printStackTrace();
					}
				}	
				
				System.out.println("[OK]");

			}//for
			
			return(bisulfited_GA_references);
	
		}
		
		/**
		 * Performs the CtoT bisulfitation of the read sequences
		 * @param files the list of read files to be bisulfited, in fasta format
		 * @param path path to the read files	
		 * @param workingDirectory where the files generated during the execution are stored
		 * @param removeUnconvertedReads trims the reads that were not correctly bisulfited converted, like reads that contain CCCC, CHHCHHCHHCHH or CHGCHGCHGCHG
		 * @param checkForAdapters confirm to check or not to check for the presence of adapter sequences in the reads
		 * @param _5primeAdapter 5' adapter
		 * @param _3primeAdapter 3' adapter			 
		 * @throws FileNotFoundException
		 * @return
		 */
		public static String[] performReadsInsilicoBisulfitation_CtoT(String[] files, StringBuilder path,StringBuilder workingDirectory,StringBuilder outputDirectory,boolean removeUnconvertedReads,boolean removeReadsWithUnconvertedBarCodes,boolean checkForAdapters,String _5primeAdapter,String _3primeAdapter){	
			
			String thisLine;
			String[] bisulfitedReads=new String[files.length];
			for(int i=0;i<bisulfitedReads.length;i++) bisulfitedReads[i]="bisulfited_CT_"+files[i].replace('/','_');
			
			for(int i=0;i<files.length;i++){
				if(files[i].indexOf("bisulfited_")==-1){
										
					System.out.print("Performing CtoT in-silico bisulfitation for "+files[i]+"...... ");
					
					BufferedWriter wr = null;
					BufferedWriter unConvertedReads = null;
					BufferedWriter unConvertedBarcodes = null;
					String inputFile=new StringBuilder(path).append(files[i]).toString();
					String outputFile=new StringBuilder(workingDirectory).append(bisulfitedReads[i]).toString();
					String unconverted_reads=new StringBuilder(outputDirectory).append("unconvertedReads_").append(files[i].replace('/','_')).toString();
					String unconverted_barcodes=new StringBuilder(outputDirectory).append("unconvertedBarCodes_").append(files[i].replace('/','_')).toString();

					try {

						
						// si existe ya el fichero lo borro
						if(Tools.doesThisFileExist(outputFile)) Tools.deleteFile(outputFile);
						
						wr = new BufferedWriter(new FileWriter(outputFile));
						unConvertedReads = new BufferedWriter(new FileWriter(unconverted_reads));
						unConvertedBarcodes= new BufferedWriter(new FileWriter(unconverted_barcodes));
						
						
						
						try{
							BufferedReader br = new BufferedReader(new FileReader(inputFile));
							String lineModified;
							
							while ((thisLine = br.readLine()) != null){
						         if(thisLine.startsWith("@")){
						        	 	// primero adjunto en la cabecera fasta la read original
						        	 
						        	 	// leo la secuencia de la read
						        	 	String originalSequence=br.readLine();
						        	 	boolean usarRead=true;
						        	 	
						        	 	if(removeReadsWithUnconvertedBarCodes){
					        	 			// obtengo el barcode del nombre del archivo, ejemplo: ES_LIF_s_8_TGtATT-sequence.txt
						     	 			String barcode=new String((files[i].split("-"))[0]);
						     	 			String aux[]=barcode.split("_");
						     	 			barcode=aux[aux.length-1];
					        	 			// localizo la posicion de la 't' en el barcode, posicion que vendria de una c no metilada
						     	 			int barcodePosition=-1; //inicializo
					        	 			barcodePosition=barcode.indexOf("t");	
					        	 			String[] tokens = thisLine.split("[#]");
					        	 			if (tokens.length==2){
						        	 			String thisReadBarcode=(thisLine.split("[#]"))[1];
						        	 			
						        	 			if(thisReadBarcode.charAt(barcodePosition)!='t' && thisReadBarcode.charAt(barcodePosition)!='T') usarRead=false;
						        	 			
						        	 		
						        	 			if(barcodePosition==-1){
						        	 				System.err.println("[ERROR]: no barcodes defined with --barcodes\n");
						        	 				System.exit(1);
						        	 			}
					        	 			}
						        	 	}//if(removeReadsWithUnconvertedBarCodes)
						        	 	
						        	 	// en caso de que quiera comprobar si la read está mal bisulfitada y eliminarla, entonces voy por aquí
						        	 	if(removeUnconvertedReads){
						        	 		if(originalSequence.matches(".*[Cc][^Gg].*[Cc][^Gg].*[Cc][^Gg].*[Cc][^Gg].*")){// apunto esta read en el fichero de las que estan mal convertidas por el bisulfito
						        	 			unConvertedReads.write(thisLine);
						        	 			unConvertedReads.newLine();
						        	 			unConvertedReads.write(originalSequence);
						        	 			unConvertedReads.newLine();
						        	 			// por ultimo copio las dos ultimas lineas del fastQ
						        	 			unConvertedReads.write(br.readLine());
						        	 			unConvertedReads.newLine();
						        	 			unConvertedReads.write(br.readLine());
						        	 			unConvertedReads.newLine();
						        	 			
						        	 			usarRead=false;
						        	 		}
						        	 		
						        	 	}
						        	 	
						        	 	if(usarRead){//que no quiere comprobar si la read está mal bisulfitada y pasa de eliminarla, vamos por aquí
							        	 	//genero la nueva cabecera fasta
							        	 	lineModified=(new StringBuilder(thisLine.replace(' ','_'))).append("||").append(originalSequence).toString();
							        	 	// la escribo
							        	 	wr.write(lineModified);
							        	 	wr.newLine();
							        	 	
							        	 	
							        	 	// ahora bisulfito la secuencia
							        	 	lineModified=originalSequence.replace('C', 'T');
							        	 	lineModified=lineModified.replace('c', 't');
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
							
							br.close();
							
						}
						catch (IOException e) {
							       System.err.println("Error: " + e);
						}
						

					} catch (IOException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}finally {
						//Close the BufferedWriter
						try {
							if (wr != null) {
								wr.flush();
								wr.close();
							}
							if(unConvertedReads != null){
								unConvertedReads.flush();
								unConvertedReads.close();
								
								
								// como en este caso siempre se crea un fichero asociado al buffered writer (aunque esté vacío, si el usuario no pidió esta opción),
								// me aseguro de borrarlo si está vacío
								if(!removeUnconvertedReads) Tools.deleteFile(unconverted_reads);
							}
							if(unConvertedBarcodes!=null){
								unConvertedBarcodes.flush();
								unConvertedBarcodes.close();
							}
						} catch (IOException ex1) {
							ex1.printStackTrace();
						}
					}						
					
					System.out.println("[OK]");

				}//if(files[i].indexOf("bisulfited_")==-1)
			}//for(int i=0;i<files.length;i++)

			return(bisulfitedReads);
			
		}		
		
		/**
		 * Computes the error rates for each sample from the methylated barcodes in the reads
		 * @param files
		 * @param path
		 * @param workingDirectory
		 * @param outputDirectory
		 * @param barcodes
		 */
		public static void computeErrorRatesFromBarcodes(String[] files, StringBuilder path,StringBuilder workingDirectory,StringBuilder outputDirectory){
			
				String thisLine;
//				String[] bisulfitedReads=new String[files.length];
//				for(int i=0;i<bisulfitedReads.length;i++) bisulfitedReads[i]="bisulfited_CT_"+files[i];
			
				for(int i=0;i<files.length;i++){
//					if(files[i].indexOf("bisulfited_")==-1){
						
						String auxiliarFileName=files[i].replace('/','_');
					
						System.out.print("Computing barcode error rates for "+files[i]+"...... ");
						double errorInSample=-1d; // inicializo a -1 (flag para indicar que no ha podido calcularse el error)					
					
						BufferedWriter unConvertedBarcodes = null;
						try {
							String inputFile=new StringBuilder(path).append(files[i]).toString();
							
							/* en caso de que el nombre de la read venga con un subdirectorio (ej: control/ES_LIF....)
							 lo que hago es cambiar la '/' de subdirectorio por un guion bajo para mantener que es una muestra control
							 a partir de aqui							
							*/
							
							File errorFile = new File(new StringBuilder(outputDirectory).append(auxiliarFileName).append(".errorRate").toString());
							if (errorFile.exists() && errorFile.lastModified() > new File(inputFile).lastModified()){
								System.out.println("Skipping, already exists [OK]");
								continue;
							}
							
							String unconverted_barcodes=new StringBuilder(outputDirectory).append("unconvertedBarCodes_").append(auxiliarFileName).toString();
						
							unConvertedBarcodes= new BufferedWriter(new FileWriter(unconverted_barcodes));
						
							// para contar el numero de  reads con barcode correctamente convertido por el bisulfito
							// y el número de reads con barcode mal convertido, de aqui extraemos el error de conversion 
							// de bisulfito que hay en cada muestra
							int correctBarcodeConversion=0;
							int failedBarcodeConversion=0;
						
						
							try{
								BufferedReader br = new BufferedReader(new FileReader(inputFile));
															
								while ((thisLine = br.readLine()) != null){
							         if(thisLine.startsWith("@")){					        	 
							        	
							        	 // obtengo el barcode del nombre del archivo, ejemplo: ES_LIF_s_8_TGtATT-sequence.txt
						     	 		String barcode=new String((files[i].split("-"))[0]);
						     	 		String aux[]=barcode.split("_");
						     	 		barcode=aux[aux.length-1];
					        	 		// localizo la posicion de la 't' en el barcode, posicion que vendria de una c no metilada
						     	 		int barcodePosition=-1; //inicializo
					        	 		barcodePosition=barcode.indexOf("t");			 
							        	
				        	 			String thisReadBarcode=(thisLine.split("[#]"))[1];
				        	 			
				        	 			if(barcodePosition==-1){
				        	 				System.err.println("[ERROR]: no barcodes defined with --barcodes\n");
				        	 				System.exit(1);
				        	 			}
				        	 			
				        	 			//System.out.println(barcode[i]+"--------------"+thisReadBarcode+"----------------------"+thisLine);
				        	 			
				        	 			
				        	 			// AL LORO CON EL CALCULO DE ERROR RATES:
				        	 			// error Rate=n. errores/(n. errores+n. aciertos)						        	 			
				        	 			// Se cuenta como error de bisulfito la presencia de 'C/c' en la posicion a mirar en el barcode.						        	 			
				        	 			// Se cuenta como error de secuenciación la presencia de 'A' o 'G' en la posicion a mirar del barcode.
				        	 			// Se cuenta como read correctamente convertida cuando hay una 'T/t' en la posicion a mirar del barcode.
				        	 			// SEGUN ESTO: el n. de errores=n. reads con error en bisulfito+n. reads con error en secuenciación
				        	 			// *** de esta manera lo calculó lister, ver supplementary info página 23
				        	 			// *** en los datos que le di a Orlando de error de bisulfito (obtenidos con el script de Perl) hay ligeras
				        	 			// diferencias porque consideré como error sólo la presencia de 'C' en la posición a mirar del barcode, con lo que
				        	 			// si un barcode estaba mal secuenciado en esa posición (A ó G) se contó como acierto (como si tuviese T) y 
				        	 			// eso no es correcto						        	 			
				        	 			if(thisReadBarcode.charAt(barcodePosition)=='t' || thisReadBarcode.charAt(barcodePosition)=='T') correctBarcodeConversion++;
				        	 			else{
				        	 				// la tasa de
				        	 				failedBarcodeConversion++;
				        	 				unConvertedBarcodes.write(thisLine);
				        	 				unConvertedBarcodes.newLine();	
				        	 			}
							         }//if(thisLine.startsWith("@"))		
								} // end while	
								
								br.close();							
							}
							catch (IOException e) {
								System.err.println("Error: " + e);
							}
						
							// estimacion del error
							errorInSample=((double)failedBarcodeConversion)/((double)(correctBarcodeConversion+failedBarcodeConversion));
							BufferedWriter errorFileWriter = new BufferedWriter(new FileWriter(errorFile));
							errorFileWriter.write("Error rate = "+ errorInSample);
							errorFileWriter.flush();
							errorFileWriter.close();
						
						} catch (IOException e1) {
							// TODO Auto-generated catch block
							e1.printStackTrace();
						}finally {
							//Close the BufferedWriter
							try {								
								if(unConvertedBarcodes!=null){
									unConvertedBarcodes.flush();
									unConvertedBarcodes.close();
								}
							} catch (IOException ex1) {
								ex1.printStackTrace();
							}
						}						
						
						System.out.println("[OK]");
						System.out.println("[Error rate calculated in bisulfite conversion for "+auxiliarFileName+" = "+errorInSample+"]");
					
					//}//if(files[i].indexOf("bisulfited_")==-1)
				}//for(int i=0;i<files.length;i++)

		}
		
		
		/**
		* Creates the Bowtie indexes for the reference
		* @param bisulfitedRefs set of reference files
		* @param working directory directory with the references
		* @param bowtieDirectory directory with the bowtie program files
		*/
		public static void performBowtieReferenceIndexing(String[] bisulfitedRefs,StringBuilder workingDirectory,StringBuilder bowtieDirectory){	

			//lo primero indexar las referencias
			for(int i=0;i<bisulfitedRefs.length;i++){
				String ref=bisulfitedRefs[i];

				System.out.print("Bowtie: indexing "+ref+"...... ");
				String [] command=new String[]{(new StringBuilder(bowtieDirectory)).append("bowtie-build").toString(),(new StringBuilder(workingDirectory)).append(ref).toString(),(new StringBuilder(workingDirectory)).append(ref).toString()};
				
				int result=Tools.executeProcessWait(command,null,null);	
				
				if(result==0) System.out.println("[OK]");
				else{
					System.err.print("\nError: ");
					for(int j=0;j<command.length;j++) System.err.print(command[j]+" ");
					System.err.println("\n");
				}

			}//for(int i=0;i<bisulfitedRefs.length;i++)

		}

		/**
		 * Extracts the name of the different samples (subirs inside READS directory) in case that there were. Otherwise returns the names of the fastqs	
		 * @param rawReadFileNames
		 * @return
		 */
		public static String[] extractSubdirNames(String [] rawReadFileNames){
			
			//recibe como entrada un array con:
			//a) valores como "ES_LIF.txt,ES_RA.txt",
			// o bien
			//b) valores como "control/ES_LIF.txt,experiment/ES_RA.txt"
			
			// en el caso a) debe retornar en el string el mismo nombre de cada fastq
			// en el caso b) debe retornar los nombres de cada subdir
			StringBuffer nombres=new StringBuffer("");
			
			for(int i=0;i<rawReadFileNames.length;i++){
				// si tiene subdir
				if(rawReadFileNames[i].indexOf('/')!=-1){
					String sample=(rawReadFileNames[i].split("/"))[0];
					// me apunto el nombre de la muestra si no lo había incluido todavía
					if(nombres.indexOf(sample)==-1) nombres.append(sample).append(','); 
				}				
			}
			
			// hay que joderse, si NO hay subdirectorios con distintas muestras dentro de READS, devuelve los ficheros de reads del fastq tal cual entran
			// en bisulfitedReads, hay que hacerlo así.
			if(nombres.length()>0) return(nombres.toString().split(","));
			else return(rawReadFileNames);
		}
		
		/**
		 * Return the names of the error rate files for each one of the input fastqs. In the case that there are several fastqs per sample it calculates
		 * also a barcode-error-rate common to all the fastqs of the same sample.
		 */
		public static String[] getErrorRateFileNames(String [] rawReadFileNames,StringBuilder readsDirectory,StringBuilder outputDirectory,String controlGenome){
			//recibe como entrada un array con:
			//a) valores como "ES_LIF.txt,ES_RA.txt",
			// o bien
			//b) valores como "control/ES_LIF.txt,experiment/ES_RA.txt"
			
			// en el caso a) debe retornar el mismo nombre de cada fastq adjuntando el sufijo ".errorRate"	-> 	ES_LIF.txt.errorRate	
			// en el caso b) debe retornar el nombre de la muestra adjuntando el sufijo ".errorRate" -> control.errorRate
			String[] errorRateFileNames=null;
			StringBuilder auxMuestras=new StringBuilder("");
			
			/*compruebo si ya he incluido esta muestra (este subdir): como itera sobre cada fastq del subdirectorio => se añadiria el nombre
			 de la misma muestra varias veces, para evitarlo meto el siguiente control*/
			String muestrasIncorporadas=new String("");
			for(int i=0;i<rawReadFileNames.length;i++){
				//caso de que hay muestras (subdirs) de la forma "control/ES_LIF.txt,experiment/ES_RA.txt"
				if(rawReadFileNames[i].indexOf('/')!=-1){
					// si todavía no he incluido la presente muestra, entonces la incluyo
					if(muestrasIncorporadas.indexOf((rawReadFileNames[i].split("/"))[0])==-1){
						auxMuestras.append((rawReadFileNames[i].split("/"))[0]).append(".errorRate").append(",");
						muestrasIncorporadas=(new StringBuilder(rawReadFileNames[i])).append(muestrasIncorporadas).toString();
					}
				}
				// si no hay varios fastqs por muestra (subdirs), tipo "ES_LIF.txt,ES_RA.txt"
				else auxMuestras.append(rawReadFileNames[i]).append(".errorRate").append(",");
				
			}//for(int i=0;i<rawReadFileNames.length;i++)

			errorRateFileNames=auxMuestras.toString().split(",");
			
			//NOTA: en el caso de que haya varios fastqs para una misma muestra, es justo en este paso donde va a calcular
			// el error conjunto de todos los fastqs de esa muestra (como si fuera un único fastq)
			
			// para saber si hay subdirs de muestras me llega con mirar si el primero de ellos tiene barra
			// indicadora de subdirectorio '/'
			if(rawReadFileNames[0].indexOf('/')!=-1){
				
				// recupero los nombres de las muestras (nombres de los subdirs)
				StringBuffer samples=new StringBuffer("");
				for(int i=0;i<rawReadFileNames.length;i++){
					// compruebo una vez más que sea una muestra (subdir) con varios fastqs
					if(rawReadFileNames[i].indexOf('/')!=-1){
						String sample=(rawReadFileNames[i].split("/"))[0];
						// me apunto el nombre de la muestra si no lo había incluido todavía
						if(samples.indexOf(sample)==-1) samples.append(sample).append(','); 
					}				
				}
				
				// recuperados los nombres de las muestras tengo que recoger los nombres de los fastqs de cada muestra
				String [] samplesArray=samples.toString().split(",");
				
				for(int i=0;i<samplesArray.length;i++){
					StringBuilder fastqs=new StringBuilder("");
					
					// recorro las muestras
					for(int j=0;j<rawReadFileNames.length;j++){
						// apunto si el fastq es de esta muestra
						if(rawReadFileNames[j].indexOf(samplesArray[i])!=-1) fastqs.append(rawReadFileNames[j]).append(",");						
					}
					
					// finalmente calculo el 'barcode error' de toda la muestra
					// recibe como argumentos el nombre del fichero de salida que almacena el error rate común, los fastqs 
					// pertenecientes a esa muestra, el subdirectorio de las reads y el outputdir
					if(controlGenome==null) MethylSeqOperations.calculateFullBarCodeError(samplesArray[i]+".errorRate",fastqs.toString().split(","),readsDirectory,outputDirectory);
					
				}
				
			
			}//if(rawReadFileNames[0].indexOf('/')!=-1)
			
			
			return(errorRateFileNames);
		}

		// recibe como argumentos el nombre del fichero de salida que almacena el error rate común y los fastqs de esa muestra
		// el barcode de cada fastq lo obtiene del nombre del archivo
		//this.calculateFullBarCodeError(sample+".errorRate",fastqs.toString().split(","));
		/* 
		 * Calculates the 'barcode error rate' for a set of fastqs of the same sample.
		 * By now the barcode of each fastq is read directly from the name of the fastq file (waiting for a better solution):
		 * the fastq file name must be like for example : ES_LIF_s_8_TGtATT-sequence.txt
		 * 
		 */
		private static void calculateFullBarCodeError(String commonOutputFileForErrorRate,String [] fastqFiles,StringBuilder readsDirectory,StringBuilder outputDirectory){
	
			// fichero de salida con el error
			File errorFile = new File(new StringBuilder(outputDirectory).append(commonOutputFileForErrorRate).toString());
			String thisLine;
			// para contar el numero de  reads con barcode correctamente convertido por el bisulfito
			// y el número de reads con barcode mal convertido, de aqui extraemos el error de conversion 
			// de bisulfito que hay en cada muestra
			int correctBarcodeConversion=0;
			int failedBarcodeConversion=0;
			double errorInSample=-1d; // inicializo a -1 (flag para indicar que NO ha podido calcularse el error)			
			String inputFile=null;
			
			// cada posicion de fastqFiles contiene algo como esto: "control/ES_LIF_s_8_TGtATT-sequence.txt,experiment/ES_RA_s_8_TAGtTT-sequence.txt"
			for(int i=0;i<fastqFiles.length;i++){
					
					System.out.print("Computing barcode error rates for "+fastqFiles[i]+"...... ");
					BufferedWriter unConvertedBarcodes = null;
										
					try {
						inputFile=new StringBuilder(readsDirectory).append(fastqFiles[i]).toString();
						
						
						/* no tengo que comprobar si el fichero de error es más reciente que el de reads, porque lo voy a calcular
						 de cualquier manera
						if (errorFile.exists() && errorFile.lastModified() > new File(inputFile).lastModified()){
							System.out.println("Skipping, already exists [OK]");
							continue;
						}*/
						
						// el fichero de unconverted barcodes lo dejo independiente para cada fastq
						String unconverted_barcodes=new StringBuilder(outputDirectory).append("unconvertedBarCodes_").append(fastqFiles[i].replace('/','_')).toString();
					
						unConvertedBarcodes= new BufferedWriter(new FileWriter(unconverted_barcodes));
					
						try{
							BufferedReader br = new BufferedReader(new FileReader(inputFile));
														
							while ((thisLine = br.readLine()) != null){
						         if(thisLine.startsWith("@")){					        	 
						        	// obtengo el barcode del nombre del archivo, ejemplo: ES_LIF_s_8_TGtATT-sequence.txt
				     	 			String barcode=new String((fastqFiles[i].split("-"))[0]);
				     	 			String aux[]=barcode.split("_");
				     	 			barcode=aux[aux.length-1];
			        	 			// localizo la posicion de la 't' en el barcode, posicion que vendria de una c no metilada
			        	 			int barcodePosition=barcode.indexOf("t");
			        	 			String thisReadBarcode=(thisLine.split("[#]"))[1];
			        	 			
			        	 			//System.out.println(barcode[i]+"--------------"+thisReadBarcode+"----------------------"+thisLine);
			        	 			
			        	 			
			        	 			// AL LORO CON EL CALCULO DE ERROR RATES:
			        	 			// error Rate=n. errores/(n. errores+n. aciertos)						        	 			
			        	 			// Se cuenta como error de bisulfito la presencia de 'C/c' en la posicion a mirar en el barcode.						        	 			
			        	 			// Se cuenta como error de secuenciación la presencia de 'A' o 'G' en la posicion a mirar del barcode.
			        	 			// Se cuenta como read correctamente convertida cuando hay una 'T/t' en la posicion a mirar del barcode.
			        	 			// SEGUN ESTO: el n. de errores=n. reads con error en bisulfito+n. reads con error en secuenciación
			        	 			// *** de esta manera lo calculó lister, ver supplementary info página 23
			        	 			// *** en los datos que le di a Orlando de error de bisulfito (obtenidos con el script de Perl) hay ligeras
			        	 			// diferencias porque consideré como error sólo la presencia de 'C' en la posición a mirar del barcode, con lo que
			        	 			// si un barcode estaba mal secuenciado en esa posición (A ó G) se contó como acierto (como si tuviese T) y 
			        	 			// eso no es correcto						        	 			
			        	 			if(thisReadBarcode.charAt(barcodePosition)=='t' || thisReadBarcode.charAt(barcodePosition)=='T') correctBarcodeConversion++;
			        	 			else{
			        	 				// la tasa de
			        	 				failedBarcodeConversion++;
			        	 				unConvertedBarcodes.write(thisLine);
			        	 				unConvertedBarcodes.newLine();	
			        	 			}
						         }//if(thisLine.startsWith("@"))		
							} // end while	
							
							br.close();							
						}
						catch (IOException e) {
							System.err.println("Error: " + e);
						}
					

					
					} catch (IOException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}finally {
						//Close the BufferedWriter
						try {								
							if(unConvertedBarcodes!=null){
								unConvertedBarcodes.flush();
								unConvertedBarcodes.close();
							}
						} catch (IOException ex1) {
							ex1.printStackTrace();
						}
					}						
					
					System.out.println("[OK]");
					
			}//for(int i=0;i<fastqFiles.length;i++)

			// estimacion del error
			errorInSample=((double)failedBarcodeConversion)/((double)(correctBarcodeConversion+failedBarcodeConversion));
			try{
				BufferedWriter errorFileWriter = new BufferedWriter(new FileWriter(errorFile));
				errorFileWriter.write("Error rate = "+ errorInSample);
				errorFileWriter.flush();
				errorFileWriter.close();
				System.out.println("[Error rate calculated in bisulfite conversion for "+errorFile+" = "+errorInSample+"]");
			}catch(IOException e){
				System.err.println("Error: " + e);			
			}			
			
		}
		
		
		/**
		 * align with Bowtie against the CtoT reference
		 * @param alignmentOutputFile name of the bowtie output file: in case of several replicates of the same sample they are mapped together to gain read depth and there is just one output common for all of them
		 * @param bisulfitedRefs set of bisulfited reference files
		 * @param bisulfitedReads set of bisulfited read files
		 * @param workingDirectory directory directory with the references
		 * @param bowtieDirectory directory with the bowtie program files
		 * @param outputDirectory directory where the alignments are going to be saved
		 * @param bisulfitedRefs
		 * @param p (Bowtie parameter, see Bowtie Manual)
		 * @param e (Bowtie parameter, see Bowtie Manual)
		 * @param l (Bowtie parameter, see Bowtie Manual)
		 * @param n (Bowtie parameter, see Bowtie Manual)
		 * @param k (Bowtie parameter, see Bowtie Manual)
		 */
		
		public static void performBowtieAlignmentwithSAMoutput(final String alignmentOutputFile, String ref,String [] bisulfitedReads,StringBuilder workingDirectory,StringBuilder bowtieDirectory,StringBuilder outputDirectory,String p,String e,String l,String n,String k,String m,String chunk,String solexaQ){
			// ahora a alinear:
			// para cada read bisulfitada
			String reads = null; 
			StringBuilder readsBuilder = new StringBuilder();
			
			//se alinean todas las reads hacia un mismo sam
			for(int i=0;i<bisulfitedReads.length;i++){
				if (i>0){
					readsBuilder.append(",");
				}
				final String read=bisulfitedReads[i].replace('/','_');
				readsBuilder.append(workingDirectory).append(read);
			}
			reads = readsBuilder.toString();
			

			// para cada referencia bisulfitada
			//for(int j=0;j<bisulfitedRefs.length;j++){
			//	String ref=bisulfitedRefs[j];
				
			//	String alignmentFile=new StringBuilder(outputDirectory).append(read).append("__").append(ref).append(".aligned").toString();
				
				System.out.print("Bowtie: aligning ["+reads+"] against ["+ref+"]...... (see .log file)...... ");
				
				
				String [] command=new String[]{(new StringBuilder(bowtieDirectory)).append("bowtie").toString(),"-t",chunk.split(" ")[0],chunk.split(" ")[1],"-p",p,solexaQ,"-e",e,"-l",l,"-n",n,"-k",k,"-m",m,"-S","--best","--nomaqround",(new StringBuilder(workingDirectory)).append(ref).toString(),new StringBuilder(reads).toString(),alignmentOutputFile};
									
				// escribo en un fichero los logs de bowtie: salida estandar+salida error
				class SynchronizedOutputStream extends OutputStream{
					private OutputStream delegate;
					public SynchronizedOutputStream(OutputStream delegate) {
						this.delegate = delegate;
					}
					
					public synchronized void write(int b) throws IOException {
						delegate.write(b);
						
					}
					
				};
				
				String outFile=alignmentOutputFile+".log";
				
				try {
					SynchronizedOutputStream out = new SynchronizedOutputStream(new FileOutputStream(new File(outFile)));
					int result=Tools.executeProcessWait(command,out,out);
					
					if(result==0) System.out.println("[OK]");
					else{
						System.err.print("\nError: ");
						for(int h=0;h<command.length;h++) System.err.print(new StringBuilder(command[h]).append(" "));
						System.err.println("\n");
					}
					
				} catch (FileNotFoundException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}

				System.out.println("\tBowtie output file: "+alignmentOutputFile);
			//}//for(int j=0;j<refs.length;j++)


	
		}
		public static void performBowtieAlignmentwithAmbigousTaggedSAMoutput(final String alignmentOutputFileCT, final String alignmentOutputFileGA, String refCT, String refGA,String [] bisulfitedReads,final StringBuilder workingDirectory,final StringBuilder bowtieDirectory,StringBuilder outputDirectory,String p,final String e,final String l,final String n,final String k,final String m,final String chunk,final String solexaQ) throws IOException{
			
		
			int threads = Integer.parseInt(p)/2;
			if (threads == 0) threads = 1;
			
			final List<File> readsFiles = new LinkedList<File>();
			
			for(int i=0;i<bisulfitedReads.length;i++){
				final String read=bisulfitedReads[i].replace('/','_');
				readsFiles.add(new File(workingDirectory+read));
			}
			List<BufferedReader> streamsCT = FastqSplitter.splitfastq(readsFiles, threads);
			List<BufferedReader> streamsGA = FastqSplitter.splitfastq(readsFiles, threads);
			//final String read=bisulfitedReads[i].replace('/','_'); 

			// para cada referencia bisulfitada
			
	
		
			
			
			abstract class LineProcessor{
				public abstract void processLine(String line);
			}
			
			
			
			class AlignerThread extends Thread{
				String ref;
				LineProcessor out;
				private String logFileName;
				private BufferedReader readsStream;
				private boolean nohead;
				boolean shouldStop = false;
				private Strand strand;
				public AlignerThread(String ref, Strand strand, LineProcessor out, String logfileName, BufferedReader readsStream, boolean nohead){
					this.ref = ref;
					this.out = out;
					this.logFileName = logfileName;
					this.readsStream = readsStream;
					this.nohead = nohead;
					this.strand = strand;
				}
				public void run(){
					
					
					System.out.println("Bowtie: aligning "+readsFiles.toString()+" against ["+ref+"]...... (see .log file)...... ");
					
					String [] command = null;
					if (nohead){
						command=new String[]{new StringBuilder(bowtieDirectory).append("bowtie").toString(), "-t",chunk.split(" ")[0],chunk.split(" ")[1],"--mm",solexaQ,"-e",e,"-l",l,"-n",n,"-k",k,"-m",m,"-S","--sam-nohead","--sam-RG","ID:"+strand.name()+"\tSM:noname","--best","--nomaqround",(new StringBuilder(workingDirectory)).append(ref).toString(),"-"};
					}else{
						command=new String[]{new StringBuilder(bowtieDirectory).append("bowtie").toString(), "-t",chunk.split(" ")[0],chunk.split(" ")[1],solexaQ,"-e",e,"-l",l,"-n",n,"-k",k,"-m",m,"-S","--best","--sam-RG","ID:"+strand.name()+"\tSM:noname","--nomaqround",(new StringBuilder(workingDirectory)).append(ref).toString(),"-"};				
					}
					
					
					String outFile=logFileName;
					
					
					try {
						FileOutputStream outLog = new FileOutputStream(new File(outFile));
						
						final Process process=Tools.executeProcess(command,null,outLog);
						
						BufferedReader stdOut = new BufferedReader(new InputStreamReader(process.getInputStream()));
						
						String line = null;
						
						//thread to feed lines
						new Thread(){
							public void run() {			
								PrintStream ps = new PrintStream(process.getOutputStream());
								
								String readsLine = null;
								System.out.println("Start read feeding to alignment against "+ref+" log file:"+logFileName);
								//long feeded = 0;
								//long processed = 0;
								try {
									/*String[] reads = new String[]{
											"@SRR019571.1_ECKER08_55:7:1:19:1204_length=53||AGTGGGTAATATTTTTAAAAATAATGTGTAAGAGTATTATGAGAAAGAAGAGG",
											"AGTGGGTAATATTTTTAAAAATAATGTGTAAGAGTATTATGAGAAAGAAGAGG","+SRR019571.1 ECKER08_55:7:1:19:1204 length=53",
											"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII=IIIIIIIG?HIIIIBIEII",
											"@SRR020275.240_SALK_2021:7:1:3:590_length=87||GNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
											"\nGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
											"\n+SRR020275.240 SALK_2021:7:1:3:590 length=87",
											"\n#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"};


									
									for(String readsLine: reads){*/
								//	File temp = File.createTempFile("caca", "fastq");
									//temp.deleteOnExit();
									
									//PrintStream tempOut = new PrintStream(new FileOutputStream(temp));
									while ((readsLine=readsStream.readLine())!=null && !shouldStop){
										
									//	System.out.println(readsLine);
										//long startFeed = System.currentTimeMillis();
										
										
										/* DEBUG: IGNORE FEED*/
									//	if (feeded >= 360000){
											
										ps.println(readsLine);
									//	}
										
										//long endFeed = System.currentTimeMillis();
									//	timewaiting += endFeed - startFeed;
										//tempOut.println(readsLine);
										
										//feeded++;
										/*if (feeded==80000){
											System.out.println("forced feed finish, after 20000 reads");
											break;
										}*/
									/*	if (feeded % 1000 == 0){
								//			System.out.println("Feeded "+feeded/4+" reads. Time total time "+((System.currentTimeMillis()-start))+" msecs. Waiting for bowtie "+timewaiting+" msecs");
											System.out.println("Feeded "+feeded+" reads.");
										}*/
										
									}
							//		tempOut.flush();
							//		tempOut.close();
									ps.flush();
									ps.close();
									
									System.out.println("Finished read feeding to alignment against "+ref+" log file:"+logFileName);
								} catch (Exception e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
							};
						}.start();
						
						try {
							while ((line=stdOut.readLine())!=null){
								//System.out.println(line);
								out.processLine(line);
								//processed++;
								
								/*if (processed % 1000 == 0){
									System.out.println("feeded "+feeded/4+" reads. processed: "+processed+ " (not processed: "+(feeded/4-processed)+")");
								}*/
							}
						} catch (IOException e1) {
							// TODO Auto-generated catch block
							e1.printStackTrace();
						}
						
						shouldStop = true; //bowtie sends a null output, so the input feed should stop
					} catch (FileNotFoundException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}
					
				}
			}
		
			final PrintStream outCT = new PrintStream(new FileOutputStream(alignmentOutputFileCT));
			final PrintStream outGA = new PrintStream(new FileOutputStream(alignmentOutputFileGA));
			class AlignerPostprocessor{
				private String CTLine=null;
				private String GALine=null;
				
				private int tagCount = 0;
				private int mergeCount = 0;
				
				StringBuilder outputBufferCT = new StringBuilder(100000);
				StringBuilder outputBufferGA = new StringBuilder(100000);
				
				public void close(){
					flushBuffer();
					System.out.println("Both alignments finished. Ambigous reads: "+tagCount);
					
				}
				public void merge(){
					CTLine = replaceOriginalRead(CTLine).trim();
					GALine = replaceOriginalRead(GALine).trim();
					
					
					outputBufferCT.append(CTLine);
					outputBufferGA.append(GALine);
					
					if (!CTLine.startsWith("@")){
						String[] tokensCT = CTLine.split("\t",7);
						String[] tokensGA = GALine.split("\t",7);
						
						//System.out.println(tokensCT[0]+" = "+tokensGA[0]);
						if (!tokensCT[0].equals(tokensGA[0])){
							throw new RuntimeException("BUG: reading two samrecords from CT and GA alignments with are a different read\nCT:"+CTLine+"\nGA:"+GALine);
						}
						//System.out.println("tokensCT[5] "+tokensCT[5]);
						//System.out.println("tokensGA[5] "+tokensGA[5]);
						if (!tokensCT[5].equals("*") && !tokensGA[5].equals("*")){
							//adding a flag to sam indicating that this is ambiguous. Also add RG (mandatory in gatk)
							outputBufferCT.append("\tZA:A:Y\tRG:Z:"+Strand.WATSON.name());
							outputBufferGA.append("\tZA:A:Y\tRG:Z:"+Strand.CRICK.name());
							
							tagCount++;
						}else{
							outputBufferCT.append("\tRG:Z:"+Strand.WATSON.name());
							outputBufferGA.append("\tRG:Z:"+Strand.CRICK.name());
							
						}
						
						
					}
					
					outputBufferCT.append("\n");
					outputBufferGA.append("\n");
					
					mergeCount++;
					
					if (mergeCount % 1000==0){
						flushBuffer();
						System.out.println(new Date()+": Merged "+mergeCount+" alignment lines");
					}
					
					CTLine = null;
					GALine = null;
					//System.out.println("merged!");
					
				}
				
				private void flushBuffer() {
					synchronized(outCT){
						outCT.print(outputBufferCT.toString());
						outCT.flush();
						outputBufferCT.setLength(0);

						outGA.print(outputBufferGA.toString());
						outGA.flush();
						outputBufferGA.setLength(0);
					}
					
				}
				public String replaceOriginalRead(String samline){
					if(!samline.startsWith("@")){
						//System.out.println(samline);
		        	 	final String[] tokens=samline.split("[\t]");
		        	 	// recupero la secuencia inicial
		        	 	final String [] firstColumn=tokens[0].split("[|][|]");
		        	 	String originalRead=firstColumn[1];
		        	 	final StringBuilder lineModified=new StringBuilder();
		        	 	
		        	 	String strand=tokens[1];
		        	 	try{
		        	 	if((Integer.parseInt(strand) & 0x0010)==0x0010){
		        	 		originalRead=getReverseComplementary(originalRead);
		        	 	}
		        	 	}catch(NumberFormatException e){
		        	 		System.out.println(samline);
		        	 		throw e;
		        	 	}
		        	 	
		        	 	for(int j=0;j<tokens.length;j++){
		        	 		// le adjunto la cabecera original de la read
		        	 		if(j==0) lineModified.append(firstColumn[0]).append("\t");
		        	 		// le adjunto la read original en lugar de la que tenia
		        	 		else if(j==9) lineModified.append(originalRead).append("\t");
		        	 		
		        	 		//else if(j==tokens.length-1) lineModified=new StringBuilder(lineModified).append(tokens[j]);
		        	 		else lineModified.append(tokens[j]).append("\t");
		        	 		
		        	 	}
		        	 	return lineModified.toString();
					}//if(!thisLine.startsWith("@"))
					return samline;
				}
				public LineProcessor CTProcessor = new LineProcessor(){

					@Override
					public void processLine(String line) {
						synchronized(AlignerPostprocessor.this){
							while(CTLine!=null){
								try {
									AlignerPostprocessor.this.wait(10000);
									if (CTLine!=null){
										System.out.println("awaking, but CTLine is still not null (if you see this message continously, bowtie may be not responding), it is: "+CTLine+"\nProcessing new CT line: "+line);
									}
								} catch (InterruptedException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
							}
							CTLine = line;
							
							if (GALine != null){
								merge();
								AlignerPostprocessor.this.notify();
							}
						}
						
					}
					
				};
				
				public LineProcessor GAProcessor = new LineProcessor(){

					@Override
					public void processLine(String line) {
						synchronized(AlignerPostprocessor.this){
							while(GALine!=null){
								try {									
									AlignerPostprocessor.this.wait(10000);
									if (GALine!=null){
										System.out.println("awaking, but GALine is still not null, it is: "+GALine+"\nProcessing new GA line: "+line);
									}
								} catch (InterruptedException e) {
								
									e.printStackTrace();
								}
							}
							GALine = line;
							
							if(CTLine!=null){
								merge();
								AlignerPostprocessor.this.notify();
							}
						}
						
					}
					
				};
				
				
			}
			
			//write the header of the sam doing a "dummy alignment"
			AlignerPostprocessor dummyposprocessor = new AlignerPostprocessor();
			LineProcessor dummyctProcessor = dummyposprocessor.CTProcessor;
			LineProcessor dummygaProcessor = dummyposprocessor.GAProcessor;
			AlignerThread dummyThreadCT = new AlignerThread(refCT, Strand.WATSON, dummyctProcessor, alignmentOutputFileCT+"_p_head.log", new BufferedReader(new InputStreamReader(new ByteArrayInputStream(new byte[0]))), false);
			AlignerThread dummyThreadGA = new AlignerThread(refGA, Strand.CRICK, dummygaProcessor, alignmentOutputFileGA+"_p_head.log", new BufferedReader(new InputStreamReader(new ByteArrayInputStream(new byte[0]))), false);
			dummyThreadCT.start(); 
			dummyThreadGA.start(); 
			
			try {
				dummyThreadCT.join();
				dummyThreadGA.join();
			} catch (InterruptedException e2) {
				// TODO Auto-generated catch block
				e2.printStackTrace();
			}
			
			dummyposprocessor.close();
			
			List<Thread> alignerThreads = new LinkedList<Thread>();
			List<AlignerPostprocessor> postprocessors = new LinkedList<AlignerPostprocessor>();
			 
			for (int i = 0; i<streamsCT.size(); i++){
				AlignerPostprocessor postprocessor = new AlignerPostprocessor();
				postprocessors.add(postprocessor);
				LineProcessor ctProcessor = postprocessor.CTProcessor;
				LineProcessor gaProcessor = postprocessor.GAProcessor;
				BufferedReader streamCT = streamsCT.get(i);
				BufferedReader streamGA = streamsGA.get(i);
				AlignerThread threadCT = new AlignerThread(refCT, Strand.WATSON, ctProcessor, alignmentOutputFileCT+"_p_"+i+".log", streamCT, true);
				AlignerThread threadGA = new AlignerThread(refGA, Strand.CRICK, gaProcessor, alignmentOutputFileGA+"_p_"+i+".log", streamGA, true);
				
				threadCT.start();
				threadGA.start();
				
				alignerThreads.add(threadCT);
				alignerThreads.add(threadGA);
				
			}
			
			for (Thread t : alignerThreads){
				try {
					t.join();
					
				} catch (InterruptedException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
			}
			for (AlignerPostprocessor postprocessor: postprocessors){
				postprocessor.close();
			}
			outCT.flush();
			outGA.flush();
			outCT.close();
			outGA.close();
			
			
			
		}
		
		
		/**
		 * Replace each line in the Bowtie output with the corresponding original read (only the ones aligned to the CtoT reference): initially the Bowtie output contains the alignment for the bisulfited reads, so the only thing that is done here
		 * is to replace the bisulfited reads with the original ones in order to run Dani's GEMtools
		 * @param reads set of read files
		 * @param outputDirectory directory where the alignments are going to be saved
		 * @param readNames 
		 */			
		public static void replaceBowtieOutputWithOriginalReadsSAMoutput(StringBuilder outputDirectory, String[] readNames){
			// abro cada fichero de salida de bowtie y leo linea a linea para hacer el reemplazo correspondiente
			String[] reads=Tools.ls(outputDirectory.toString(),"aligned");
			for(int i=0;i<reads.length;i++){
				String inputFile=new StringBuilder(outputDirectory).append(reads[i]).toString();
				System.out.println(inputFile);
				//process only selected reads (readNames)
				boolean found = false;
				for (String readName: readNames){
					if (inputFile.contains(readName)){
						found =true;
						break;
					}
				}
				if(!found) continue;
				
				
				String outputFile=new StringBuilder(outputDirectory).append(reads[i]).append(".tmp").toString();						
				
				System.out.print(new StringBuilder("Replacing bisulfited CtoT reads in ").append(inputFile).append(" with original reads...... ").toString());
				BufferedWriter wr = null;
				boolean wasThisStepdone=false;
				
				try {

					
					// si existe ya el fichero lo borro
					if(Tools.doesThisFileExist(outputFile)) Tools.deleteFile(outputFile);
					
					wr = new BufferedWriter(new FileWriter(outputFile));
					
					try{
						BufferedReader br = new BufferedReader(new FileReader(inputFile));
						String thisLine;						
						
						while ((thisLine = br.readLine()) != null){
							if(!thisLine.startsWith("@")){
					        	 	String[] tokens=thisLine.split("[\t]");
					        	 	// recupero la secuencia inicial
					        	 	String [] firstColumn=tokens[0].split("[|][|]");
					        	 	StringBuilder lineModified=new StringBuilder("");
					        	 	String originalRead=null;
					        	 	try{
					        	 		originalRead=firstColumn[1];
					        	 		String strand=tokens[1];
					        	 		if((Integer.parseInt(strand) & 0x0010)==0x0010) originalRead=getReverseComplementary(originalRead);
					        	 	}catch(java.lang.ArrayIndexOutOfBoundsException e){
					        	 		System.out.println("\n*** Warning: it seems that the bowtie output files were already replaced with original reads. Leaving this step.");
					        	 		wasThisStepdone=true;
					        	 		break;					        	 		
					        	 	}
					        	 	
					        	 	for(int j=0;j<tokens.length;j++){
					        	 		// le adjunto la cabecera original de la read
					        	 		if(j==0) lineModified=new StringBuilder(lineModified).append(firstColumn[0]).append("\t");
					        	 		// le adjunto la read original en lugar de la que tenia
					        	 		else if(j==9) lineModified=new StringBuilder(lineModified).append(originalRead).append("\t");
					        	 		
					        	 		//else if(j==tokens.length-1) lineModified=new StringBuilder(lineModified).append(tokens[j]);
					        	 		else lineModified=new StringBuilder(lineModified).append(tokens[j]).append("\t");
					        	 		
					        	 	}
					        	 	
					        	 	wr.write(lineModified.toString());
					        	 	wr.newLine();
							}//if(!thisLine.startsWith("@"))
							else{
								wr.write(thisLine);
								wr.newLine();
							}
					    } // end while
						
						br.close();
						
					}
					catch (IOException e) {
						       System.err.println("Error: " + e);
					}
					
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}finally {
					//Close the BufferedWriter
					try {
						if (wr != null) {
							wr.flush();
							wr.close();
						}
					} catch (IOException ex1) {
						ex1.printStackTrace();
					}
				}	
				
				if(!wasThisStepdone) System.out.println("[OK]");
				
				// muevo el fichero temporal creado a su fichero original
				File from=new File(outputFile);
				File to=new File(inputFile);
				Tools.moveOneFileToAnother(from,to);
				
			}//for(int i=0;i<reads.length;i++)
		
		}
		

		/**
		 * Returns the reverse complementary sequence of a nucleotide sequence
		 * @param sequence
		 * @return
		 */
		public static String getReverseComplementary(String sequence){
			StringBuilder complementaria=new StringBuilder("");
			
			// 1ero construyo la complementaria
			for(int i=0;i<sequence.length();i++){
					switch(sequence.charAt(i)){
					case 'A': complementaria.append("T"); break;
					case 'T': complementaria.append("A"); break;
					case 'C': complementaria.append("G"); break;
					case 'G': complementaria.append("C"); break;
					default: complementaria.append(sequence.charAt(i));
				}
			}
			
			//ahora la reversa
			StringBuilder reversa=new StringBuilder("");
			for(int i=complementaria.length()-1;i>-1;i--){
				reversa.append(complementaria.charAt(i));				
			}

			return(reversa.toString());
		}
		
		/**
		 * Returns the list of read files to analyze
		 * 
		 * @param readsDirectory
		 * @return
		 */
		public static StringBuilder getListOfReadFiles(StringBuilder readsDirectory){
			StringBuilder reads=new StringBuilder("");
		
			/*lo primero que hago es buscar subdirectorios en el directorio READS, asumiendo que cada uno corresponde a una muestra
			 * y que dentro están los fastqs de cada muestra.
			 * Si no hay subdirectorios y solo hay ficheros, asumo que los ficheros son todos fastqs pertenecientes a la misma muestra
			*/
			
			//leo el directorio de reads
			File f=new File(readsDirectory.toString());
			String[] files=f.list();
			boolean areThereSubdirs=false;
			
			//compruebo si dentro del fichero de reads hay subdirectorios o directamente ficheros
			for(int i=0;i<files.length;i++){
				File f2=new File(readsDirectory.toString(),files[i]);
					
				// si es un subdirectorio
				if(f2.isDirectory()){
					areThereSubdirs=true;
					File subdir=new File(readsDirectory.toString(),files[i]);	
					String[] filesSubdir=subdir.list();
						
					for(int j=0;j<filesSubdir.length;j++){
						reads.append(new StringBuilder(files[i]).append("/").append(filesSubdir[j]).append(","));				
					}
						
				}
			}//for(int i=0;i<files.length;i++)
			
			// si no hay subdirectorios recojo los ficheros que hay en ese directorio como ficheros fastqs
			if(!areThereSubdirs) for(int i=0;i<files.length;i++) reads.append(new StringBuilder(files[i]).append(","));
			
			return(reads);
		}
		
		/**
		 * Checks that the files received exist
		 * @param file
		 * @param Directory
		 * @return String[]
		 */
		public static String[] searchForFiles(StringBuilder file,StringBuilder Directory,String prefix){
			// en caso de recibir varios ficheros como entrada, vendrian separados por comas
			System.out.println(new StringBuilder("Searching for file(s)...... ").append(file).toString());
			String[] files=MethylSeqOperations.getFiles(new StringBuilder(file));
			String [] output=new String[files.length];

			for(int i=0;i<files.length;i++){
				
				// si no ha puesto prefijo hago la búsqueda normal
				if(prefix==null){ 
					System.out.print(new StringBuilder(Directory).append(files[i]).append("...... ").toString());
					
					if(Tools.doesThisFileExist(new StringBuilder(Directory).append(files[i]).toString())){
						System.out.println("[OK]");
						//en este caso lo devuelvo con el mismo nombre
						output[i]=files[i];
					}
					else{
							System.err.println(new StringBuilder("[ERROR]: unable to find ").append(Directory).append(files[i]).append("\n").toString());
							System.exit(1);
					}
				}else{
					System.out.print(new StringBuilder(Directory).append(prefix).append(files[i]).append("...... ").toString());
					
					if(Tools.doesThisFileExist(new StringBuilder(Directory).append(prefix).append(files[i]).toString())){
						System.out.println("[OK]");
						//le pongo el prefijo
						output[i]=new StringBuilder(prefix).append(files[i]).toString();
					}
					else{
							System.err.println(new StringBuilder("[ERROR]: unable to find ").append(Directory).append(prefix).append(files[i]).append("\n").toString());
							System.exit(1);
					}				
				}
			}
			
			return(output);
		}
		
		/* replaces '/' by '_' en read file names when cases like the following:
		 * control/ES_LIF.txt -> control_ES_LIF.txt
		 */
		public static StringBuilder replaceSlashByHyphen(StringBuilder ReadFiles){
			// obtengo los nombres separados a un array (rompo por las ',')
			String [] ReadFileNames=MethylSeqOperations.getFiles(ReadFiles);
			
			StringBuilder finalNames=new StringBuilder("");
			
			for(int i=0;i<ReadFileNames.length;i++){
						finalNames.append(ReadFileNames[i].replace('/','_')).append(',');						
			}
			
			return(finalNames);
			
		}
		
		/* replaces '/' by '_' en read file names when cases like the following:
		 * control/ES_LIF.txt -> control_ES_LIF.txt
		 */
		public static String[] replaceSlashByHyphenIntoArray(StringBuilder ReadFiles){
			// obtengo los nombres separados a un array (rompo por las ',')
			String [] ReadFileNames=MethylSeqOperations.getFiles(ReadFiles);
				
			for(int i=0;i<ReadFileNames.length;i++){
						ReadFileNames[i]=ReadFileNames[i].replace('/','_');						
			}
			
			return(ReadFileNames);
			
		}
		
		/**
		 * Checks that the bisulfited references exist
		 * @param reference
		 * @param workingDirectory
		 * @return String[]
		 */
		public static String[] getBisulfitedReferences(StringBuilder reference,StringBuilder workingDirectory,String prefix){
			// si selecciona varias referencias, entonces vienen separados por comas
			String[] refs=MethylSeqOperations.getFiles(reference);
			
			String[] bisulfited_references=new String[refs.length];
			for(int i=0;i<refs.length;i++){
				bisulfited_references[i]=new StringBuilder(prefix).append(refs[i]).toString();
				if(!Tools.doesThisFileExist(new StringBuilder(workingDirectory).append(bisulfited_references[i]).toString())){
					System.err.println(new StringBuilder("\n[ERROR]: unable to find the bisulfited reference file ").append(workingDirectory).append(bisulfited_references[i]).toString());
					System.exit(1);							
				}//if						
			}//for
			
			return(bisulfited_references);
		}
		
		/**
		 * Checks that the index files for the bisulfited references exist
		 * @param reference
		 * @param workingDirectory
		 * @param prefix
		 * @return String[]
		 */
		public static String[] getIndexedReferences(StringBuilder reference,StringBuilder workingDirectory,String prefix){
			// si selecciona varias referencias, entonces vienen separados por comas
			String[] refs=MethylSeqOperations.getFiles(reference);
			
			// el nombre genérico en realidad es el nombre de referencia bisulfitada. Bowtie genera luego los ficheros de indexación tomando este como prefijo
			String [] indexedReferencesGenericName=new String[refs.length];
			
			for(int i=0;i<refs.length;i++){
				indexedReferencesGenericName[i]=new StringBuilder(prefix).append(refs[i]).toString();
				// si existen los indices de Bowtie para esa referencia bisulfitada
				if(!Tools.doesThisFileExist(new StringBuilder(workingDirectory).append(indexedReferencesGenericName[i]).append(".1.ebwt").toString())){
					System.err.println(new StringBuilder("\n[ERROR]: unable to find the bisulfited reference file ").append(workingDirectory).append(indexedReferencesGenericName[i]).append(".1.ebwt").toString());
					System.exit(1);							
				}//if						
			}//for
			
			return(indexedReferencesGenericName);			
		}
		
		/**
		 * Analyzes the methylation sites and creates the .methylation files in the output directory
		 */
		public static void analyzeMethylation(final StringBuilder outputDirectory, final StringBuilder workingDirectory, final StringBuilder referenceDirectory,final StringBuilder readsDirectory, String[] refNames, String[] readNames,final String readErrorRateFile, final String controlGenome, final StringBuilder samtoolsDirectory, final ErrorRateMode errorMode, final double[] errorRates,final StringBuilder bedFilesDirectory,String[] bedFiles){

			class AnalyzeMethylationThread extends Thread{
				
				File inputSamFile;
				boolean isPositive;
				String ref;
				String readErrorRateFile;
				Throwable error=null;
				public AnalyzeMethylationThread(File inputSamFile,
						boolean isPositive, String ref, String readErrorRateFile) {
					super();
					this.inputSamFile = inputSamFile;
					this.isPositive = isPositive;
					this.ref = ref;
					this.readErrorRateFile = readErrorRateFile;
				}
				//output attributes in this thread. Will have values once it finishes
				double[] p,cutoffs;
				File outputBamFile;
				File outputMethylation;
				
				
				public void run(){
					GenomeIndex index;
					try {
						index = getOrCreateIndex(new File(referenceDirectory+"/"+ref), new File(workingDirectory.toString()));
						File sorted1 = new File(inputSamFile.getAbsolutePath() + ".sorted.sam");
						sortSAM(inputSamFile, sorted1);
						this.outputBamFile = buildBAMAndIndex(sorted1, samtoolsDirectory);
						
						
						synchronized(this){
							System.out.println("Computing Error rates for "+(isPositive?"Watson":"Crick"));						

							if (errorMode == ErrorRateMode.from_barcodes){
								
								BufferedReader in = new BufferedReader(new FileReader(new File(outputDirectory.toString()+readErrorRateFile)));
								String line = in.readLine();
								line = line.substring(line.lastIndexOf('=')+1).trim();
								double error = Double.parseDouble(line);
								this.p= new double[]{error,error,error};
							}else if (errorMode == ErrorRateMode.from_control_genome){	
									this.p=new double[]{ getErrorRate(this.outputBamFile, index, controlGenome, "CG", isPositive), getErrorRate(this.outputBamFile, index, controlGenome, "CHG", isPositive), getErrorRate(this.outputBamFile, index, controlGenome, "CHH", isPositive)};
									
							}else{
									this.p = errorRates;
							}
							System.out.println("Error rates for "+(isPositive?"Watson":"Crick")+Arrays.toString(this.p));	
							this.notifyAll();
						}
						
						System.out.println("Calling methylCytosines for "+(isPositive?"Watson":"Crick"));
						this.outputMethylation = analyzeMethylationSAM(this.outputBamFile, index, outputDirectory, this.p[0], this.p[1], this.p[2], controlGenome,isPositive);
						System.out.println("Finished methylCytosines for "+(isPositive?"Watson":"Crick"));
						
						System.out.println("Computing p-val cutoffs for "+(isPositive?"Watson":"Crick"));	
						StringBuffer cutOffDetails = new StringBuffer();
						this.cutoffs = computePValCutoffs(this.outputMethylation,cutOffDetails);						
						System.out.println("Finsihed p-val cutoffs for "+(isPositive?"Watson":"Crick"));
						System.out.println(cutOffDetails.toString());
						
						
					} catch (IOException e) {
						this.error = e;
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						this.error = e;
					}
				
				}
				
			}
			try {
				List<IntervalsIndex> bedIndexes = new LinkedList<IntervalsIndex>();
				
				
				if(bedFiles!=null) for (String bed: bedFiles){					
					bedIndexes.add(IntervalsIndexFactory.createIntervalsIndex(bed, new File(bedFilesDirectory+bed), 1, 2, 3, 1, false));
				}
				
				for (String ref: refNames){
//					for (String reads : readsNames){
						/*File samFile1 = new File(outputDirectory.toString()+"/bisulfited_CT_"+reads+"__bisulfited_CT_"+ref+".aligned");
						File samFile2 = new File(outputDirectory.toString()+"/bisulfited_CT_"+reads+"__bisulfited_GA_"+ref+".aligned");*/
						File samFile1 = new File(readNames[0]);
						File samFile2 = new File(readNames[1]);
						
						
						
						System.out.println("Analysing methylation: "+ref+" "+readNames[0]+" (MIN_DEPTH="+DEPTH_FILTER+", "+" DISCARD_NOT_CT="+DISCARD_NON_CT+", FDR= "+FDR+", REMOVE_CLONAL= "+REMOVE_CLONAL+", CORRECT_NONCG="+CORRECT_NONCG+")");
						System.out.println("Analyzing methylation: "+ref+" "+readNames[1]+" (MIN_DEPTH="+DEPTH_FILTER+", "+" DISCARD_NOT_CT="+DISCARD_NON_CT+", FDR= "+FDR+", REMOVE_CLONAL= "+REMOVE_CLONAL+", CORRECT_NONCG="+CORRECT_NONCG+")");
						
						AnalyzeMethylationThread analysisCT = new AnalyzeMethylationThread(samFile1, true, ref, readErrorRateFile);
						AnalyzeMethylationThread analysisGA = new AnalyzeMethylationThread(samFile2, false, ref, readErrorRateFile);
						
						
						analysisCT.start();
						analysisGA.start();
						
						synchronized(analysisCT){
							if(analysisCT.p ==null){
								analysisCT.wait();
							}
						}
						synchronized(analysisGA){
							if(analysisGA.p ==null){
								analysisGA.wait();
							}
						}
						
						analysisCT.join(); if (analysisCT.error!=null) throw new RuntimeException(analysisCT.error);
						analysisGA.join(); if (analysisGA.error!=null) throw new RuntimeException(analysisGA.error);
						
						
						
						/*//CT
						File sorted1 = new File(samFile1.getAbsolutePath() + ".sorted.sam");
						sortSAM(samFile1, sorted1);
						File bamCT = buildBAMAndIndex(sorted1, samtoolsDirectory);
						
						
						double p1_CG = getErrorRate(bamCT, index, controlGenome, "CG", true);
						double p1_CHG = getErrorRate(bamCT, index, controlGenome, "CHG", true);
						double p1_CHH = getErrorRate(bamCT, index, controlGenome, "CHH", true);
						System.out.println("[OK]");
						System.out.println("Error rates for: "+samFile1.getName()+ " are CG: "+p1_CG+", CHG: "+p1_CHG+", CHH: "+p1_CHH);
						File methylationCT = analyzeMethylationSAM(bamCT, index, outputDirectory, p1_CG, p1_CHG, p1_CHH, controlGenome,true);
						double[] cutoffsCT = computePValCutoffs(methylationCT);

						//GA
						File sorted2 = new File(samFile2.getAbsolutePath() + ".sorted.sam");
						sortSAM(samFile2, sorted2);						
						File bamGA = buildBAMAndIndex(sorted2, samtoolsDirectory);

						
						double p2_CG = getErrorRate(bamGA, index, controlGenome, "CG", false);
						double p2_CHG = getErrorRate(bamGA, index, controlGenome, "CHG", false);
						double p2_CHH = getErrorRate(bamGA, index, controlGenome, "CHH", false);
						System.out.println("[OK]");
						System.out.println("Error rates for: "+samFile2.getName()+ " are CG: "+p2_CG+", CHG: "+p2_CHG+", CHH: "+p2_CHH);
						File methylationGA = analyzeMethylationSAM(bamGA, index, outputDirectory, p2_CG, p2_CHG, p2_CHH, controlGenome,false);
						double[] cutoffsGA = computePValCutoffs(methylationGA);
						*/
						
						writeMethylcytosines(analysisCT.outputBamFile, analysisGA.outputBamFile,
								analysisCT.cutoffs, analysisGA.cutoffs, analysisCT.p, analysisGA.p,outputDirectory, bedIndexes,bedFiles);
					
//					}
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
		

		/**
		 * Analyzes the methylation sites and creates the .methylation files in the output directory
		 * @param outputDirectory
		 * @param referenceDirectory
		 * @param refsNames
		 * @param readNames
		 * @throws IOException
		 * @throws InterruptedException
		 */
		public static void analyzeMethylationGATK(final StringBuilder outputDirectory, final StringBuilder workingDirectory, final StringBuilder referenceDirectory,final StringBuilder readsDirectory, String[] refsNames, String[] readNames,final String readErrorRateFile,final String controlGenome, final StringBuilder samtoolsDirectory, final ErrorRateMode errorMode, final double[] errorRates, int nThreads, boolean trim, boolean ambiguous, boolean removebad,final StringBuilder bedFilesDirectory,String[] bedFiles){
			
			try {
				for (String ref: refsNames){
					/*for (String reads : readNames){
						File samFileCT = new File(outputDirectory.toString()+"/bisulfited_CT_"+reads+"__bisulfited_CT_"+ref+".aligned");
						File samFileGA = new File(outputDirectory.toString()+"/bisulfited_CT_"+reads+"__bisulfited_GA_"+ref+".aligned");*/
						File samFileCT = new File(readNames[0]);
						File samFileGA = new File(readNames[1]);
						
						File fasta = new File(referenceDirectory+"/"+ref);
					
						File sortedCT = new File(samFileCT.getAbsolutePath() + ".sorted.sam");
						sortSAM(samFileCT, sortedCT);
						File outputBamFileCT = buildBAMAndIndex(sortedCT, samtoolsDirectory);
						
						File sortedGA = new File(samFileGA.getAbsolutePath() + ".sorted.sam");
						sortSAM(samFileGA, sortedGA);
						File outputBamFileGA = buildBAMAndIndex(sortedGA, samtoolsDirectory);
						
						RuntimeMXBean runtimemxBean = ManagementFactory.getRuntimeMXBean();
					
						
						String command = "java -Xmx1024M -cp "+runtimemxBean.getClassPath()+" org.broadinstitute.sting.gatk.CommandLineGATK -T ListerMethylation -I "+outputBamFileCT+" -I "+outputBamFileGA+" -R "+fasta+" --controlgenome "+controlGenome+" -nt "+nThreads+" --outdir "+outputDirectory+" --fdr "+FDR;
						
						if (bedFiles!=null) for (String bedfile : bedFiles){
							command+=" -annotation:"+bedfile+",bed "+bedFilesDirectory+bedfile;
						}
						
						command+=" --read_filter Lister";
						if (ambiguous){
							command+=" --removeambiguous";
						}
						if (trim){
							command+=" --trim";
						}
						if (removebad){
							command+=" --removebad";
						}
						if (errorMode == ErrorRateMode.from_barcodes){
							BufferedReader in = new BufferedReader(new FileReader(new File(outputDirectory.toString()+"/"+readErrorRateFile)));
							String line = in.readLine();
							line = line.substring(line.lastIndexOf('=')+1).trim();
							double error = Double.parseDouble(line);
							command+=" --errorrate "+error+","+error;
							
						}
						System.out.println("GATK command: "+command);
						final Process p = Runtime.getRuntime().exec(command.split(" "));
						/*Runtime.getRuntime().addShutdownHook(new Thread(){
							@Override
							public void run() {
								p.destroy();
							}
						});*/
						class StreamReader extends Thread{
							private InputStream stream;
							private PrintStream out;
							public StreamReader(InputStream stream, PrintStream out) {
								this.stream = stream;
								this.out = out;
								
							}
							
							@Override
							public void run() {
								BufferedReader reader = new BufferedReader(new InputStreamReader(this.stream));
								
								String line = null;
								
								try {
									while ((line = reader.readLine())!=null){
										out.println(line);
									}
								} catch (IOException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
							}
							
						}
						Thread stdOut = new StreamReader(p.getInputStream(), System.out);
						Thread stdErr = new StreamReader(p.getErrorStream(), System.err);
						stdOut.start();
						stdErr.start();
						p.waitFor();
						
						stdOut.join();
						stdErr.join();
						if (p.exitValue()!=0){
							throw new RuntimeException("Methylation analysis has non zero return value: "+p.exitValue());
						}
						
						
						
					
					//}
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
		
	
		private static void sortSAM(File sam, File output) throws InterruptedException, IOException{
			//sort the sam
			
			if(new File(output.getAbsolutePath()).exists() && new File(output.getAbsolutePath()).lastModified() > sam.lastModified()){
				//System.out.println("skipping. found the corresponding sorted file, older than the unsorted input file");
				return;
			}
			System.out.println("Sorting " + sam.getAbsolutePath());
			String outfile = sam.getAbsolutePath() + ".sorted.sam";
			String command = "java -Xmx1024M -jar ./lib/SortSam.jar I="+ sam.getAbsolutePath()+ " O="+ outfile+" SO=coordinate TMP_DIR="+sam.getAbsoluteFile().getParentFile().getAbsolutePath();
			if (Runtime.getRuntime().exec((command).split(" ")).waitFor()!=0){
				//error
				if (new File(outfile).exists()){
					new File(outfile).delete();
				}
				throw new RuntimeException("Sort failed! command: "+command);
			}
			
			
		}
		
		
		private static GenomeIndex getOrCreateIndex(File file, File directory) throws IOException {
			File indexFile =new File(directory+"/"+file.getName()+".genindex"); 
			if (!indexFile.exists()){
				System.out.println("creating index "+indexFile.getAbsolutePath());
				GenomeIndexBuilder.buildGenome(file, indexFile, null);
			}
			try{
				GenomeIndex index =  GenomeIndexFactory.createGenomeIndex(indexFile);
				return index;
			}catch(Exception e){
				System.err.println("[ERROR]: Genome index file seems to be corrupted: "+indexFile.getAbsolutePath()+". Re-building index.");
				indexFile.delete();
				return getOrCreateIndex(file, directory);
			}
			
			
		}
		
		private static double getErrorRate(File bamFile, GenomeIndex reference,
				String controlGenome, String type, boolean isPositive) throws IOException {
			
			
			//iterate over the sorted sam
			
			
			//do we have default error rates?
			if (type.equalsIgnoreCase("CG") && default_error_rates[0]!=-1d){
				return default_error_rates[0];
			}else if (type.equalsIgnoreCase("CHG") && default_error_rates[1]!=-1d){
				return default_error_rates[1];
			}else if (type.equalsIgnoreCase("CHH") && default_error_rates[2]!=-1d){
				return default_error_rates[2];
			}
			
			
			SAMFileReader bamReader = new SAMFileReader(bamFile);
			
			IntervalList controlInterval = new IntervalList(bamReader.getFileHeader());
			if (bamReader.getFileHeader().getSequence(controlGenome) == null){
				System.err.println("[ERROR]: Control genome not found: "+controlGenome);
				System.exit(1);
			}
			controlInterval.add(new Interval(controlGenome, 1, bamReader.getFileHeader().getSequence(controlGenome).getSequenceLength()));
			SamLocusIterator samit = new SamLocusIterator(bamReader, controlInterval, true);
			samit.setEmitUncoveredLoci(false);
			Iterator<LocusInfo> it = samit.iterator();

			int totalC=0;
			int totalCT=0;
			

	
			int positions=0;
			
			while (it.hasNext()) {
			
				
				LocusInfo locus = it.next();
				
				
				List<RecordAndOffset> bases = locus.getRecordAndPositions();
				if (REMOVE_CLONAL){
					bases = filterClonalReads(bases);
				}
				if (bases.size() < DEPTH_FILTER) continue;
				
				
				int position = locus.getPosition();
				
				
				try{
					reference.getSequenceLength(controlGenome);
				}catch(NullPointerException e){
					System.err.println("[ERROR] Sequence not found in index: "+controlGenome);
					throw e;
				}
				String _type = getContextType(reference, controlGenome, position, isPositive, it);
				
				
				if (_type==null || !_type.equals(type)){
					continue;
				}
				
				boolean isC = false;
				boolean isG = false;
				
				
				if (Character.toUpperCase(reference.getSingleNucleotide(controlGenome, position))== 'C'){
					isC = true;
				}else if (Character.toUpperCase(reference.getSingleNucleotide(controlGenome, position))== 'G'){
					isG = true;
				}
			
				positions++;
				
			
				for (RecordAndOffset s : bases) {
					if (isPositive && isC){
						if (s.getReadBase()=='C'){
								//System.err.println("+Error\t"+controlGenome+"\t"+position+"\t"+_type+"\t"+s.getRecord().getAlignmentStart()+"\t"+s.getRecord().getReadString());
								totalC++;
								totalCT++;
							
						}else if (s.getReadBase()=='T' || !DISCARD_NON_CT){
							//System.err.println("+OK\t"+controlGenome+"\t"+position+"\t"+_type+"\t"+s.getRecord().getAlignmentStart()+"\t"+s.getRecord().getReadString());
							
							totalCT++;
							
						
						}
					}else if (!isPositive && isG){
					//	if (type.equals("CHH")){System.out.println((char)s.getReadBase());}
						if (s.getReadBase()=='G'){
							//System.err.println("-Error\t"+controlGenome+"\t"+position+"\t"+_type+"\t"+s.getRecord().getAlignmentStart()+"\t"+s.getRecord().getReadString());
							totalC++;
							totalCT++;
							
						
						}else if (s.getReadBase()=='A' || !DISCARD_NON_CT){
							//System.err.println("-OK\t"+controlGenome+"\t"+position+"\t"+_type+"\t"+s.getRecord().getAlignmentStart()+"\t"+s.getRecord().getReadString());
							
							totalCT++;
							
							
						}
					}
				}
		
        	 		
					
			}
		
			//return ratiosum/(double)positions;
			return (double)totalC/(double)totalCT;

		}

		private static String getContextType(GenomeIndex reference, String sequence, int position, boolean isPositive, Iterator<LocusInfo> it) throws IOException{
			
			long sequenceLength;
			try{
				sequenceLength = reference.getSequenceLength(sequence);
			}catch(NullPointerException e){
				System.err.println("[ERROR] Sequence not found in index: "+sequence);
				throw e;
			}
			String type=null;
			if (isPositive){
				if (position<sequenceLength){ //we must be before the last nucleotide
					//System.out.println(Character.toUpperCase(reference.getSingleNucleotide(sequence, position)));
					if (Character.toUpperCase(reference.getSingleNucleotide(sequence, position))== 'C'){						 
						if (Character.toUpperCase(reference.getSingleNucleotide(sequence, position+1))=='G'){
							//CG site
							type = "CG";
						}else if (Character.toUpperCase(reference.getSingleNucleotide(sequence, position+1))!='N'){
							if (position < sequenceLength-1){
								if (Character.toUpperCase(reference.getSingleNucleotide(sequence, position + 2))=='G'){
									//CHG site
									type = "CHG";
									
								}else if (Character.toUpperCase(reference.getSingleNucleotide(sequence, position + 2))!='N'){
									//CHH site
									type = "CHH";
									
								}
							}
							/*if (CORRECT_NONCG){
								if (type!=null){
									//correct to CG
									LocusInfo possibleNext = currentIterator.lookAhead();
									if (possibleNext.getPosition()==position+1){
										int gCount = 0;
										List<RecordAndOffset> bases = possibleNext.getRecordAndPositions();
										for (RecordAndOffset record : bases){
											if (record.getReadBase()=='G'){
												gCount++;
											}
										}
										if ((double)gCount/(double)bases.size()>=CONSENSUS_MIN_PERCENT){
	//										System.err.println("corrected "+type+" to CG at "+sequence+":"+position);
											type="CG";
										}
									}
								}
							}*/
							
						}
					}
					
				}
			}else{
				if (position>1){ //we must be before the last nucleotide
					if (Character.toUpperCase(reference.getSingleNucleotide(sequence, position))== 'G'){						 
						if (Character.toUpperCase(reference.getSingleNucleotide(sequence, position-1))=='C'){
							//CG site
							type = "CG";
						}else if (Character.toUpperCase(reference.getSingleNucleotide(sequence, position-1))!='N'){
							if (position > 2){
								if (Character.toUpperCase(reference.getSingleNucleotide(sequence, position - 2))=='C'){
									//CHG site
									type = "CHG";
									
								}else if (Character.toUpperCase(reference.getSingleNucleotide(sequence, position - 2))!='N'){
									//CHH site
									type = "CHH";
									
								}
							}
							/*if(CORRECT_NONCG){
								if (type!=null){
									//correct to CG
									LocusInfo possiblePrevious= currentIterator.getPrevious();
									if (possiblePrevious!=null && possiblePrevious.getPosition()==position-1){
										int gCount = 0;
										List<RecordAndOffset> bases = possiblePrevious.getRecordAndPositions();
										for (RecordAndOffset record : bases){
											if (record.getReadBase()=='G'){
												gCount++;
											}
										}
										if ((double)gCount/(double)bases.size()>=CONSENSUS_MIN_PERCENT){
											System.err.println("corrected "+type+" to CG at "+sequence+":"+position);
											type="CG";
										}
									}
								}
							}*/
						}
					}
					
				}
			}
			
			return type;
		}
		
		private static File analyzeMethylationSAM(File sam, GenomeIndex reference, StringBuilder outputDirectory, double p_CG, double p_CHG, double p_CHH, String controlGenome, boolean isPositive) throws IOException, InterruptedException{
		
			
			
			File outFile = new File(outputDirectory.toString()+sam.getName()+".methylation");
			
			//outFile.deleteOnExit();
			PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(outFile)));
			
			
			
			out.println("#CHR\tPOS\tMETH_TYPE\tDEPTH\tC_COUNT\tC_T_DEPTH\t%C\tADJ_%C\tBASES_IN_THIS_POSITION\tP-VAL");
			
			//iterate over the sorted sam
			SAMFileReader inputSam = new SAMFileReader(sam);	
			
			
			SamLocusIterator samit = new SamLocusIterator(inputSam);
			samit.setEmitUncoveredLoci(false);
			
			Iterator<LocusInfo> it = samit.iterator();
			
			
			int positions=0;
			int depthsum=0;
			while (it.hasNext()) {
				LocusInfo locus = it.next();
				List<RecordAndOffset> bases = locus.getRecordAndPositions();
				if (REMOVE_CLONAL){
					bases = filterClonalReads(bases);
				}
				
				if (bases.size()<DEPTH_FILTER){
					continue;
				}
				
				
				String sequence = locus.getSequenceName();
				if (sequence.equalsIgnoreCase(controlGenome)){
					continue;
				}
				int position = locus.getPosition();
				
				//check the type of methyl site (CG, CHG, CHH) 						
				String type = getContextType(reference, sequence, position, isPositive, it);
				
				
				if (type !=null){ //we are in a possible methylated site				
					positions++;
					//count C's
					int totalC=0;
					int totalCT=0;
					int sumCQual=0;
					int sumCTQual=0;
					StringBuilder readedBases = new StringBuilder();
					depthsum+=bases.size();
					for (RecordAndOffset s : bases) {
						
						readedBases.append((char)s.getReadBase());
						if (isPositive){
							if (s.getReadBase()=='C'){
								sumCQual+=s.getBaseQuality();
								totalC++;
	
								totalCT++;
								sumCTQual+=s.getBaseQuality();
							}else if (s.getReadBase()=='T' || !DISCARD_NON_CT){
								totalCT++;
								sumCTQual+=s.getBaseQuality();
							}
						}else{
							if (s.getReadBase()=='G'){
								sumCQual+=s.getBaseQuality();
								totalC++;
	
								totalCT++;
								sumCTQual+=s.getBaseQuality();
							}else if (s.getReadBase()=='A' || !DISCARD_NON_CT){
								totalCT++;
								sumCTQual+=s.getBaseQuality();
							}
						}
					}							        	 	//wr.write(br.readLine());
	        	 	//wr.newLine();	
					String outLine = "";
					
					
					double p;
					if (type.equals("CG")){
						p = p_CG;
					}else if (type.equals("CHG")){
						p = p_CHG;
					}else{
						p = p_CHH;
					}
					BinomialDistribution binomial = new BinomialDistributionImpl(totalCT, p);
					
						double pval;
						try {
							pval = (totalC==0)? 1.0d: (1.0d-binomial.cumulativeProbability(totalC-1));
							//DecimalFormat format = new DecimalFormat();
							//format.setMaximumFractionDigits(100);
							outLine = 	locus.getSequenceName() + "\t" + 
							locus.getPosition() + "\t" + 
							type + "\t" + 
							bases.size() +"\t" + 
							totalC + "\t" + 
							totalCT + "\t" +
							(((double)totalC/(double)totalCT)*100) +"\t" + 
							(((double)sumCQual/(double)sumCTQual)*100)+ "\t" + 
							readedBases.toString() + "\t"+
							(pval);
							
						} catch (MathException e) {
							
							e.printStackTrace();
							throw new RuntimeException(e);
						}
						
						
						
									
					
					out.println(outLine);

					
					
					//debug
					//System.err.println(outLine);
					
				}
							
			}
			out.flush();
			out.close();
			System.out.println(" [OK]");
			System.out.println("Mean depth at reference C: "+((double)depthsum/(double)positions));
			//computing the best p-value cutoffs to 1% FDR
			return outFile;
			
			
			
			
			
			//summary
		/*	DecimalFormat format2 = new DecimalFormat();
			format2.setMaximumFractionDigits(2);
			PrintStream summaryOut = new PrintStream(new FileOutputStream(new File(outFile.getAbsolutePath()+".methylation.summary")));
			int cTotal = cCount[0]+cCount[1]+cCount[2];
			int mcTotal = mcCount[0]+mcCount[1]+mcCount[2];
			
			summaryOut.println("Error rates on "+controlGenome+" control genome [CG, CHG, CHH]: ["+p_CG+", "+p_CHG+", "+p_CHH+"]");
			summaryOut.println("Discard not CT bases: "+DISCARD_NOT_CT);
			summaryOut.println("DEPTH_FILTER: >="+DEPTH_FILTER);
			summaryOut.println("REMOVE_CLONAL: "+REMOVE_CLONAL);
			summaryOut.println("CORRECT_NONCG: "+CORRECT_NONCG);
			summaryOut.println("Total sites (methylated or not): "+cTotal);
			summaryOut.println("p-val cutoffs for "+FDR*100+"% FDR [CG, CHG, CHH]: "+Arrays.toString(cutoffs));
			summaryOut.println("\nMethylated sites: "+(mcCount[0]+mcCount[1]+mcCount[2])+" ("+format2.format(((double)(mcCount[0]+mcCount[1]+mcCount[2])/(double)cTotal)*100)+"%) ("+
					format2.format(((double)mcCount[0]/(double)mcTotal)*100)+"% CG, "+
					format2.format(((double)mcCount[1]/(double)mcTotal)*100)+"% CHG, "+
					format2.format(((double)mcCount[2]/(double)mcTotal)*100)+"% CHH)"			
			);
			summaryOut.println("\nPer-context counts:");
			summaryOut.println("CG: total: "+cCount[0]+" methylated "+mcCount[0]+" ("+format2.format(((double)mcCount[0]/(double)cCount[0])*100)+"%)");
			summaryOut.println("CHG: total: "+cCount[1]+" methylated "+mcCount[1]+" ("+format2.format(((double)mcCount[1]/(double)cCount[1])*100)+"%)");
			summaryOut.println("CHH: total: "+cCount[2]+" methylated "+mcCount[2]+" ("+format2.format(((double)mcCount[2]/(double)cCount[2])*100)+"%)");
			
			summaryOut.close();
		*/
		
			/*
			//compute adjusted p-values with FDR Benjamini Hochberg
			//sort the file
			
			File sortedFile = new File(outFile.getAbsolutePath()+".sorted");
			
			
			//sortedFile.deleteOnExit();
			System.out.print("Sorting p-values to compute adjusted p-values");
			int ret = Tools.executeProcessWait(("sort -nr -k 9,9 -T "+outputDirectory+" "+outFile.getAbsolutePath()).split(" "), new FileOutputStream(sortedFile),System.err);
			if (ret!=0){
				throw new RuntimeException("Couldn't sort");
			}
			System.out.println("[OK]");

			System.out.print("Computing adjusted p-values");
			File adjustedFile = new File(outFile.getAbsolutePath()+".methylation");
			PrintStream outAdjusted = new PrintStream(adjustedFile);
			
			outAdjusted.println("#CHR\tPOS\tMETH_TYPE\tDEPTH\tC_COUNT\tC_T_DEPTH\t%C\tADJ_%C\tP-VAL\tREADED_BASES\tADJ-PVAL");
			BufferedReader inSorted = new BufferedReader(new FileReader(sortedFile));
			String line = null;
			
			int rank = totalSites;
			int[] methylatedSites = {0,0,0}; //CG, CHG, CHH
			int[] totals = {0,0,0}; //CG, CHG, CHH
			while((line = inSorted.readLine())!=null){
				String[] tokens = line.split("\t");
				double pval = Double.parseDouble(tokens[8]);
				double qval = pval*(totalSites/rank);
				rank--;
				outAdjusted.println(line+"\t"+qval);
				
				
					if (tokens[2].equalsIgnoreCase("CG")){
						if (qval <= FDRThreshold){
							methylatedSites[0]++;
						}
						totals[0]++;
					}
					else if(tokens[2].equalsIgnoreCase("CHG")){
						if (qval <= FDRThreshold){
							methylatedSites[1]++;
						}
						totals[1]++;
					}else{
						if (qval <= FDRThreshold){
							methylatedSites[2]++;
						}
						totals[2]++;
					}
				
				
			}
			System.out.println("[OK]");
			outAdjusted.close();
*/
			
			
		}

		private static double[] computePValCutoffs(File outFile, StringBuffer computingDetails)
				throws FileNotFoundException, IOException {
			
			
			double[] positiveRate = {0d,0d,0d};
			double[] cutoffs = {1,1,1};
			boolean[] needAdjust = {true, true, true};
			int[] cCount = {0,0,0};
			int[] mcCount = {0,0,0};
			
			int iteration = 0;
			while(needAdjust[0] || needAdjust[1] || needAdjust[2]){
				iteration ++;
				cCount = new int[]{0,0,0};
				mcCount= new int[]{0,0,0};
				BufferedReader outIn = new BufferedReader(new FileReader(outFile));
				String line = null;
				while ((line = outIn.readLine())!=null){
					String[] tokens = line.split("\t");
					if (tokens[2].equals("CG")){
						cCount[0]++;
						if (Double.parseDouble(tokens[tokens.length-1])<cutoffs[0]){
							mcCount[0]++;
						}
					}else if (tokens[2].equals("CHG")){
						cCount[1]++;
						if (Double.parseDouble(tokens[tokens.length-1])<cutoffs[1]){
							mcCount[1]++;
						}
					}else if (tokens[2].equals("CHH")){
						cCount[2]++;
						if (Double.parseDouble(tokens[tokens.length-1])<cutoffs[2]){
							mcCount[2]++;
						}
					}
				}
				outIn.close();
				
				double[] _positiveRate = {((double)mcCount[0]/(double) cCount[0])*100, ((double)mcCount[1]/(double) cCount[1])*100, ((double)mcCount[2]/(double) cCount[2])*100};
				computingDetails.append("Iteration "+iteration+", M: "+Arrays.toString(cutoffs)+" %mC: "+Arrays.toString(_positiveRate)+"\n");
				if (Double.isNaN(_positiveRate[0]) || Double.isNaN(_positiveRate[1]) || Double.isNaN(_positiveRate[2])){
					computingDetails.append("Nan?? Abort\n");
					System.exit(1);
				}
				if (_positiveRate[0]!=positiveRate[0]){
					cutoffs[0] = FDR * _positiveRate[0]/(100d - _positiveRate[0]);
					positiveRate[0] = _positiveRate[0];
				}else{
					needAdjust[0]=false;
				}
				if (_positiveRate[1]!=positiveRate[1]){
					cutoffs[1] = FDR * _positiveRate[1]/(100d - _positiveRate[1]);
					positiveRate[1] = _positiveRate[1];
				}else{
					needAdjust[1]=false;
				}
				if (_positiveRate[2]!=positiveRate[2]){
					cutoffs[2] = FDR * _positiveRate[2]/(100d - _positiveRate[2]);
					positiveRate[2] = _positiveRate[2];
				}else{
					needAdjust[2]=false;
				}
				computingDetails.append("\tneed Adjust: "+Arrays.toString(needAdjust)+"\n");
			}
			
			computingDetails.append("Finished p-value adjust. Result "+Arrays.toString(cutoffs)+"\n");
			return cutoffs;
		}

		private static String getAnnotations(String seq, int pos, boolean isCT, List<IntervalsIndex> bedIndexes){
			
			StringBuilder annotation=new StringBuilder();
			
			boolean firstAnnotation = true;
			for (IntervalsIndex index : bedIndexes){
				if (!firstAnnotation) annotation.append("\t");
				else firstAnnotation = false;
				
				Iterator<es.cnio.bioinfo.pileline.core.Interval> intervals = index.getIntervals(seq, pos, 0, false);
				
				boolean firstInterval = true;
				int count=0;
				while (intervals.hasNext()){
					es.cnio.bioinfo.pileline.core.Interval interval = intervals.next();
					String[] dataSplit = interval.getData().split("\t"); 
					if ((dataSplit.length<3) || ((dataSplit[2].equals("+") && isCT) || (dataSplit[2].equals("-") && !isCT))){
						if (!firstInterval){
							annotation.append(":");
						}else{
							firstInterval = false;
						}
						annotation.append(interval.getData().replaceAll("\t","__"));
						count++;
					}
				}
				if (count==0){
					annotation.append("NULL");
				}
				
				
			}
			return annotation.toString();
			
		}
		private static void writeMethylcytosines(File bamCT, File bamGA,				
				double[] cutoffsCT, double[] cutoffsGA, double[] pCT, double[] pGA,
				StringBuilder outputDirectory, List<IntervalsIndex> bedIndexes,String[] bedFiles) throws FileNotFoundException,
				IOException {
			
			
			System.out.println("Writing methylcytosines....");
			 
			File outMC = new File(outputDirectory.toString()+bamCT.getName()+"_"+bamGA.getName()+".methylcytosines");
			PrintStream outMCs = new PrintStream(new BufferedOutputStream(new FileOutputStream(outMC)));

			
			outMCs.print("#seq\tpos\tstrand\tcontext\tmCs\treads\tLister_method_methylation_state");
			if(bedFiles!=null) for(int i=0;i<bedFiles.length;i++) outMCs.print("\t"+bedFiles[i]);
			outMCs.println();
			
			SAMFileReader samCTReader = new SAMFileReader(bamCT);
			SAMFileReader samGAReader = new SAMFileReader(bamGA);
			
			
			
			BufferedReader methylcytosinesCT = new BufferedReader(new FileReader(outputDirectory.toString()+bamCT.getName()+".methylation")); 
			BufferedReader methylcytosinesGA = new BufferedReader(new FileReader(outputDirectory.toString()+bamGA.getName()+".methylation")); 
			
			String lineCT = methylcytosinesCT.readLine();
			String lineGA = methylcytosinesGA.readLine();
			
			
			
			
			
			while(lineCT.startsWith("#")){
				lineCT = methylcytosinesCT.readLine();
			}
			while(lineGA.startsWith("#")){
				lineGA = methylcytosinesGA.readLine();
			}
			while(lineCT!=null || lineGA!=null){
			
				//which next?
				String[] tokensCT = lineCT!=null?lineCT.split("\t"):null;
				String[] tokensGA = lineGA!=null?lineGA.split("\t"):null;
				
				boolean isCT = false;
				
				if (lineCT == null || lineGA == null){
					isCT = (lineCT!=null);
					
				}else{
					if (!tokensCT[0].equals(tokensGA[0])){
						if (samCTReader.getFileHeader().getSequenceIndex(tokensCT[0])<samGAReader.getFileHeader().getSequenceIndex(tokensGA[0])){
							isCT = true;
						}
						
							
					}else{
						if (Integer.parseInt(tokensCT[1])<=Integer.parseInt(tokensGA[1])){
							isCT =true;
						}
					}
				}
			
				String processLine;
				if (isCT){
					processLine = lineCT;
					lineCT = methylcytosinesCT.readLine();
				}else{
					processLine = lineGA;
					lineGA = methylcytosinesGA.readLine();
				}
			
	
				double[] cutoffs = isCT?cutoffsCT: cutoffsGA;
				
				String[] tokens = processLine.split("\t");
				boolean methylated = false;
				if (tokens[2].equals("CG")){
					
					if (Double.parseDouble(tokens[tokens.length-1])<cutoffs[0]){
						methylated = true;
						
					}
				
				}else if (tokens[2].equals("CHG")){
					if (Double.parseDouble(tokens[tokens.length-1])<cutoffs[1]){
						methylated = true;
					}
				
				}else if (tokens[2].equals("CHH")){
					if (Double.parseDouble(tokens[tokens.length-1])<cutoffs[2]){
						methylated = true;
					}
					
				}
				
				//if (methylated){
				
					String prepend = null;
					String postpend = null;
					
					if (!tokens[2].equals("CG") && CORRECT_NONCG){
					
						
						double gRatio = 0.0d;
						double cRatio = 0.0d;
						
						
						SAMFileReader thisStrandReader = isCT?samCTReader:samGAReader;
						SAMFileReader oppositeSam = isCT?samGAReader:samCTReader;
						
						int downStreamPosition = Integer.parseInt(tokens[1])+(isCT?1:-1);
						
						List<RecordAndOffset> thisBases = Tools.getReadsForLocus(isCT?Strand.WATSON:Strand.CRICK, tokens[0], downStreamPosition, false, java.util.Arrays.asList(new SAMFileReader[]{thisStrandReader}), new LinkedList<SamRecordFilter>(),'X');
						
						int gCount=0;
						
							if (thisBases.size()>0){
								
								if (REMOVE_CLONAL){
									thisBases = filterClonalReads(thisBases);
								}
								//guanine count
								
								for (RecordAndOffset base : thisBases){
									if ((isCT && base.getReadBase()=='G') || (!isCT && base.getReadBase()=='C') ){
										gCount++;
									}
								}
								
								
								gRatio = ((double)gCount/(double)thisBases.size());
								
								//siteIterator.close(); //if using index, throws exception
							}
							
							List<RecordAndOffset> oppositeBases = Tools.getReadsForLocus(isCT?Strand.CRICK:Strand.WATSON, tokens[0], downStreamPosition, false, java.util.Arrays.asList(new SAMFileReader[]{oppositeSam}), new LinkedList<SamRecordFilter>(),'X');
							boolean oppositeNotCovered=false;
							int oppositeDepth=0;
							int oppositeCitosinesCount=0;
							if (oppositeBases.size()>0){
								//covered								
								if (REMOVE_CLONAL){
									oppositeBases = filterClonalReads(oppositeBases);
								}
								oppositeDepth = oppositeBases.size();
								//citosine count
								
								for (RecordAndOffset base : oppositeBases){
									if ((isCT && base.getReadBase()=='G') || (!isCT && base.getReadBase()=='C') ){ //we use G, because opposite bases are also represented in the forward strand
										oppositeCitosinesCount++;
									}
								}
								
								
								cRatio = ((double)oppositeCitosinesCount/(double)oppositeBases.size());
								
								//oppositeIterator.close(); //if using index, throws exception
							}
							
							else oppositeNotCovered=true;
							
							if (gRatio >0.25){ //consensus at +1 downstream is G?
							
								tokens[2]="CG****";
							}
							//System.out.println((isCT?"+":"-")+"\t"+tokens[0]+"\t"+tokens[1]+"\t"+downStreamPosition+"\tgratio: "+gRatio+", cRatio: "+cRatio+", gCount: "+gCount+" , cytosines in mC: "+tokens[4]+" opposite not covered: "+oppositeNotCovered+" opp depth "+oppositeDepth);
							else if (gRatio >= 0.2 && cRatio > 0.2){

								
								// correct 
								//System.out.println("\tCorrecting... due to gRatio and cRatio are >=20%");
								//should add in the opposite? Yes if is passes the biniomial test
								
								double pval;
								double p;
								if (tokens[2].equals("CHG")){
									p = isCT?pCT[1]:pGA[1];
								}else{
									p = isCT?pCT[2]:pGA[2];
								}
								BinomialDistribution binomial = new BinomialDistributionImpl(oppositeDepth, p);
								try {
									pval = (oppositeCitosinesCount==0)? 1.0d: (1.0d-binomial.cumulativeProbability(oppositeCitosinesCount-1));
									double cutoff=0;
									if (tokens[2].equals("CHG")){
										cutoff = isCT?cutoffsCT[1]:cutoffsGA[1];
									}else{
										cutoff = isCT?cutoffsCT[2]:cutoffsGA[2];
									}
									if (pval<cutoff){
										//yes!
									
										StringBuilder lineB = new StringBuilder();
										lineB.append(tokens[0]+"\t");
										lineB.append(downStreamPosition+"\t");
										lineB.append(((isCT)?"-":"+")+"\t");
										lineB.append("CG***\t");
										lineB.append(oppositeCitosinesCount+"\t");
										lineB.append(oppositeDepth+"\t"+(methylated?"METHYLATED":"UNMETHYLATED")+"\t");
										lineB.append((getAnnotations(tokens[0],Integer.parseInt(tokens[1]),isCT,bedIndexes)));
										if (downStreamPosition > Integer.parseInt(tokens[1])){
											postpend = lineB.toString(); //add line in the opposite strand, after
										}else{
											prepend = lineB.toString(); //add line in the opposite strand, after
										}
										
									}
								} catch (MathException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
									
								tokens[2]="CG*";
							//}else if (oppositeNotCovered && gRatio >= 0.2 && Math.abs(gCount-Integer.parseInt(tokens[4]))<=2){
							}else if (oppositeNotCovered && gRatio >= 0.2 && Math.abs(gRatio - (Double.parseDouble(tokens[4])/Double.parseDouble(tokens[3])))<=0.10){
								//correct
								
								//System.out.println("\tCorrecting due to gCount and cytosines are similar");
								tokens[2]="CG**";
							}else{
								//System.out.println("\tNo correction");
							}
						
							
							
							
					}
					if (prepend!=null){
						outMCs.println(prepend);
					}
					
					StringBuilder lineB = new StringBuilder();
					
					//outMCs.println("#seq\tpos\torient\tcontext\tmCs\treads");
					lineB.append(tokens[0]+"\t");
					lineB.append(tokens[1]+"\t");
					lineB.append(((isCT)?"+":"-")+"\t");
					lineB.append(tokens[2]+"\t");
					lineB.append(tokens[4]+"\t");
					lineB.append(tokens[5]+"\t");
					lineB.append((methylated?"METHYLATED":"UNMETHYLATED")+"\t");
					lineB.append((getAnnotations(tokens[0],Integer.parseInt(tokens[1]),isCT,bedIndexes)));
					outMCs.println(lineB.toString());
					
					if (postpend!=null){
						outMCs.println(postpend);
					}
				
			}
			outMCs.flush(); outMCs.close();
		
			writeSummary(outMC,outputDirectory.toString()+bamCT.getName()+"_"+bamGA.getName()+".summary",  pCT, pGA, cutoffsCT, cutoffsGA, bedIndexes);
			System.out.println("[OK]");
		}

		private static class Results{
			int[] totalmC_CT = {0,0,0};
			int[] totalmC_GA = {0,0,0};
			int[] totalC_CT = {0,0,0};
			int[] totalC_GA = {0,0,0};
			
			int[] totalRawmC_CT = {0,0,0};
			int[] totalRawmC_GA = {0,0,0};
			int[] totalRawC_CT = {0,0,0};
			int[] totalRawC_GA = {0,0,0};
			
			int correctionsCT = 0;
			int correctionsGA = 0;
			int addedCG_CT=0;
			int addedCG_GA=0;
		}
		
		private static Results getResults(File methylcytosines, LineFilter filter) throws IOException{
			Results toret = new Results();
			BufferedReader reader = new BufferedReader(new FileReader(methylcytosines));
			String line=null;
			
			while((line = reader.readLine())!=null){
				if (line.startsWith("#")){
					continue;
				}
				if (!filter.accept(line)){
					continue;
				}
				String[] parts = line.split("\t");
				boolean isCT = parts[2].equals("+");
				int[] totalmC = isCT?toret.totalmC_CT:toret.totalmC_GA;
				int[] totalC = isCT?toret.totalC_CT:toret.totalC_GA;
				int[] totalRawmC = isCT?toret.totalRawmC_CT:toret.totalRawmC_GA;
				int[] totalRawC = isCT?toret.totalRawC_CT:toret.totalRawC_GA;
				int context=0;
				if (parts[3].indexOf("CG")!=-1){
					context=0;
				}else if (parts[3].indexOf("CHG")!=-1){
					context=1;
				}else{
					context=2;
				}
				
				totalRawC[context]+=Integer.parseInt(parts[5]);
				totalRawmC[context]+=Integer.parseInt(parts[4]);
				
				totalC[context]++;
				if (parts[6].equals("METHYLATED")){
					totalmC[context]++;
				}
				
				if (parts[3].contains("*")){
					if (isCT) toret.correctionsCT ++;
					else toret.correctionsGA ++;
				}
				
				if (parts[3].equals("CG***")){
					if (isCT) toret.addedCG_CT++;
					else toret.addedCG_GA++;
				}
			}
		
			return toret;
		}
		
		private static String writeResults(Results res){
			StringBuilder builder = new StringBuilder();
			
			builder.append("\nnon-CG to CG orrections====\n");
			builder.append("Corrections Watson: "+res.correctionsCT+" of "+(res.totalC_CT[1]+res.totalC_CT[2]+res.correctionsCT)+" ("+(double)res.correctionsCT/(double)(res.totalC_CT[1]+res.totalC_CT[2]+res.correctionsCT)+"). Added "+res.addedCG_GA+" new mCG to the GA strand\n");
			builder.append("Corrections Crick: "+res.correctionsGA+" of "+(res.totalC_GA[1]+res.totalC_GA[2]+res.correctionsGA)+" ("+(double)res.correctionsGA/(double)(res.totalC_GA[1]+res.totalC_GA[2]+res.correctionsGA)+"). Added "+res.addedCG_GA+" new mCG to the GA strand\n");
			builder.append("Corrections both strands: "+(res.correctionsGA+res.correctionsCT)+" of "+(res.totalmC_GA[1]+res.totalmC_GA[2]+res.totalmC_CT[1]+res.totalmC_CT[2])+" ("+(double)(res.correctionsGA+res.correctionsCT)/(double)(res.totalmC_GA[1]+res.totalmC_GA[2]+res.totalmC_CT[1]+res.totalmC_CT[2])+"). Added "+(res.addedCG_GA+res.addedCG_CT)+" new mCG to the opposite strand\n");

			builder.append("\n=====RESULTS (Lister)====\n");
			int totalmCytosines_CT = (res.totalmC_CT[0]+res.totalmC_CT[1]+res.totalmC_CT[2]);
			int totalCytosines_CT = (res.totalC_CT[0]+res.totalC_CT[1]+res.totalC_CT[2]);
			int totalmCytosines_GA = (res.totalmC_GA[0]+res.totalmC_GA[1]+res.totalmC_GA[2]);
			int totalCytosines_GA = (res.totalC_GA[0]+res.totalC_GA[1]+res.totalC_GA[2]);
			
			builder.append("Watson-----\n");
			
			builder.append("Global methylation: "+totalmCytosines_CT+" of "+totalCytosines_CT+" ("+((double)totalmCytosines_CT/(double)totalCytosines_CT)+")\n");
			builder.append("Per context: CG:"+((double)res.totalmC_CT[0]/(double)totalmCytosines_CT)+", CHG:"+((double)res.totalmC_CT[1]/(double)totalmCytosines_CT)+", CHH:"+((double)res.totalmC_CT[2]/(double)totalmCytosines_CT)+"\n");
			builder.append("mCG: "+res.totalmC_CT[0]+"/"+res.totalC_CT[0]+" ("+((double)res.totalmC_CT[0]/(double)res.totalC_CT[0])+")\n");
			builder.append("mCHG: "+res.totalmC_CT[1]+"/"+res.totalC_CT[1]+" ("+((double)res.totalmC_CT[1]/(double)res.totalC_CT[1])+")\n");
			builder.append("mCHH: "+res.totalmC_CT[2]+"/"+res.totalC_CT[2]+" ("+((double)res.totalmC_CT[2]/(double)res.totalC_CT[2])+")\n");
			
			
		
			
			builder.append("\nCrick-----\n");
			builder.append("Global methylation: "+totalmCytosines_GA+" of "+totalCytosines_GA+" ("+((double)totalmCytosines_GA/(double)totalCytosines_GA)+")\n");
			builder.append("Per context: CG:"+((double)res.totalmC_GA[0]/(double)totalmCytosines_GA)+", CHG:"+((double)res.totalmC_GA[1]/(double)totalmCytosines_GA)+", CHH:"+((double)res.totalmC_GA[2]/(double)totalmCytosines_GA)+"\n");
			builder.append("mCG: "+res.totalmC_GA[0]+"/"+res.totalC_GA[0]+" ("+((double)res.totalmC_GA[0]/(double)res.totalC_GA[0])+")\n");
			builder.append("mCHG: "+res.totalmC_GA[1]+"/"+res.totalC_GA[1]+" ("+((double)res.totalmC_GA[1]/(double)res.totalC_GA[1])+")\n");
			builder.append("mCHH: "+res.totalmC_GA[2]+"/"+res.totalC_GA[2]+" ("+((double)res.totalmC_GA[2]/(double)res.totalC_GA[2])+")\n");
			
			
			
			
			builder.append("\nBoth Strands-----\n");
			builder.append("Global methylation: "+(totalmCytosines_CT+totalmCytosines_GA)+" of "+(totalCytosines_CT+totalCytosines_GA)+" ("+((double)(totalmCytosines_CT+totalmCytosines_GA)/(double)(totalCytosines_CT+totalCytosines_GA))+")\n");
			builder.append("Per context: CG:"+((double)(res.totalmC_CT[0]+res.totalmC_GA[0])/(double)(totalmCytosines_CT+totalmCytosines_GA))+", CHG:"+((double)(res.totalmC_CT[1]+res.totalmC_GA[1])/(double)(totalmCytosines_CT+totalmCytosines_GA))+", CHH:"+((double)(res.totalmC_CT[2]+res.totalmC_GA[2])/(double)(totalmCytosines_CT+totalmCytosines_GA))+"\n");
			builder.append("mCG: "+(res.totalmC_CT[0]+res.totalmC_GA[0])+"/"+(res.totalC_CT[0]+res.totalC_GA[0])+" ("+((double)(res.totalmC_CT[0]+res.totalmC_GA[0])/(double)(res.totalC_CT[0]+res.totalC_GA[0]))+")\n");
			builder.append("mCHG: "+(res.totalmC_CT[1]+res.totalmC_GA[1])+"/"+(res.totalC_CT[1]+res.totalC_GA[1])+" ("+((double)(res.totalmC_CT[1]+res.totalmC_GA[1])/(double)(res.totalC_CT[1]+res.totalC_GA[1]))+")\n");
			builder.append("mCHH: "+(res.totalmC_CT[2]+res.totalmC_GA[2])+"/"+(res.totalC_CT[2]+res.totalC_GA[2])+" ("+((double)(res.totalmC_CT[2]+res.totalmC_GA[2])/(double)(res.totalC_CT[2]+res.totalC_GA[2]))+")\n");
			
			builder.append("\n=====RESULTS (per-read) ====\n");
			int totalRawmCytosines_CT = (res.totalRawmC_CT[0]+res.totalRawmC_CT[1]+res.totalRawmC_CT[2]);
			int totalRawCytosines_CT = (res.totalRawC_CT[0]+res.totalRawC_CT[1]+res.totalRawC_CT[2]);
			int totalRawmCytosines_GA = (res.totalRawmC_GA[0]+res.totalRawmC_GA[1]+res.totalRawmC_GA[2]);
			int totalRawCytosines_GA = (res.totalRawC_GA[0]+res.totalRawC_GA[1]+res.totalRawC_GA[2]);
			
			builder.append("Watson-----\n");
			
			builder.append("Global methylation: "+totalRawmCytosines_CT+" of "+totalRawCytosines_CT+" ("+((double)totalRawmCytosines_CT/(double)totalRawCytosines_CT)+")\n");
			builder.append("Per context: CG:"+((double)res.totalRawmC_CT[0]/(double)totalRawmCytosines_CT)+", CHG:"+((double)res.totalRawmC_CT[1]/(double)totalRawmCytosines_CT)+", CHH:"+((double)res.totalRawmC_CT[2]/(double)totalRawmCytosines_CT)+"\n");
			builder.append("mCG: "+res.totalRawmC_CT[0]+"/"+res.totalRawC_CT[0]+" ("+((double)res.totalRawmC_CT[0]/(double)res.totalRawC_CT[0])+")\n");
			builder.append("mCHG: "+res.totalRawmC_CT[1]+"/"+res.totalRawC_CT[1]+" ("+((double)res.totalRawmC_CT[1]/(double)res.totalRawC_CT[1])+")\n");
			builder.append("mCHH: "+res.totalRawmC_CT[2]+"/"+res.totalRawC_CT[2]+" ("+((double)res.totalRawmC_CT[2]/(double)res.totalRawC_CT[2])+")\n");
			
			
		
			
			builder.append("\nCrick-----\n");
			builder.append("Global methylation: "+totalRawmCytosines_GA+" of "+totalRawCytosines_GA+" ("+((double)totalRawmCytosines_GA/(double)totalRawCytosines_GA)+")\n");
			builder.append("Per context: CG:"+((double)res.totalRawmC_GA[0]/(double)totalRawmCytosines_GA)+", CHG:"+((double)res.totalRawmC_GA[1]/(double)totalRawmCytosines_GA)+", CHH:"+((double)res.totalRawmC_GA[2]/(double)totalRawmCytosines_GA)+"\n");
			builder.append("mCG: "+res.totalRawmC_GA[0]+"/"+res.totalRawC_GA[0]+" ("+((double)res.totalRawmC_GA[0]/(double)res.totalRawC_GA[0])+")\n");
			builder.append("mCHG: "+res.totalRawmC_GA[1]+"/"+res.totalRawC_GA[1]+" ("+((double)res.totalRawmC_GA[1]/(double)res.totalRawC_GA[1])+")\n");
			builder.append("mCHH: "+res.totalRawmC_GA[2]+"/"+res.totalRawC_GA[2]+" ("+((double)res.totalRawmC_GA[2]/(double)res.totalRawC_GA[2])+")\n");
			
			
			
			
			builder.append("\nBoth Strands-----\n");
			builder.append("Global methylation: "+(totalRawmCytosines_CT+totalRawmCytosines_GA)+" of "+(totalRawCytosines_CT+totalRawCytosines_GA)+" ("+((double)(totalRawmCytosines_CT+totalRawmCytosines_GA)/(double)(totalRawCytosines_CT+totalRawCytosines_GA))+")\n");
			builder.append("Per context: CG:"+((double)(res.totalRawmC_CT[0]+res.totalRawmC_GA[0])/(double)(totalRawmCytosines_CT+totalRawmCytosines_GA))+", CHG:"+((double)(res.totalRawmC_CT[1]+res.totalRawmC_GA[1])/(double)(totalRawmCytosines_CT+totalRawmCytosines_GA))+", CHH:"+((double)(res.totalRawmC_CT[2]+res.totalRawmC_GA[2])/(double)(totalRawmCytosines_CT+totalRawmCytosines_GA))+"\n");
			builder.append("mCG: "+(res.totalRawmC_CT[0]+res.totalRawmC_GA[0])+"/"+(res.totalRawC_CT[0]+res.totalRawC_GA[0])+" ("+((double)(res.totalRawmC_CT[0]+res.totalRawmC_GA[0])/(double)(res.totalRawC_CT[0]+res.totalRawC_GA[0]))+")\n");
			builder.append("mCHG: "+(res.totalRawmC_CT[1]+res.totalRawmC_GA[1])+"/"+(res.totalRawC_CT[1]+res.totalRawC_GA[1])+" ("+((double)(res.totalRawmC_CT[1]+res.totalRawmC_GA[1])/(double)(res.totalRawC_CT[1]+res.totalRawC_GA[1]))+")\n");
			builder.append("mCHH: "+(res.totalRawmC_CT[2]+res.totalRawmC_GA[2])+"/"+(res.totalRawC_CT[2]+res.totalRawC_GA[2])+" ("+((double)(res.totalRawmC_CT[2]+res.totalRawmC_GA[2])/(double)(res.totalRawC_CT[2]+res.totalRawC_GA[2]))+")\n");
			
			
			
			
			return builder.toString();
		}
		
		
		private static void writeSummary(File methylcytosines, String filename, double[] pCT, double[] pGA, double[] cutoffsCT, double[] cutoffsGA, List<IntervalsIndex> beds) throws IOException{
		
			PrintStream outSummary = new PrintStream(new BufferedOutputStream(new FileOutputStream(filename)));
			
			outSummary.println("Experiment date: "+new Date());
			outSummary.println("===PARAMETERS===");
			outSummary.println("Discard non-CT bases: "+DISCARD_NON_CT);
			outSummary.println("DEPTH_FILTER: >="+DEPTH_FILTER);
			outSummary.println("REMOVE_CLONAL: "+REMOVE_CLONAL);
			outSummary.println("CORRECT_NON-CG: "+CORRECT_NONCG);
			outSummary.println("FDR: "+FDR);
			
			outSummary.println("Watson error rates for Lister: "+Arrays.toString(pCT));
			outSummary.println("Crick error rates for Lister: "+Arrays.toString(pGA));
			outSummary.println("Watson p-value cutoffs for Lister: "+Arrays.toString(cutoffsCT));
			outSummary.println("Crick p-value cutoffs for Lister: "+Arrays.toString(cutoffsGA));
			
			outSummary.println("=============================");
			outSummary.println("IGNORING ANNOTATIONS RESULTS");
			outSummary.println("=============================");
			
			Results ignoring = getResults(methylcytosines, new LineFilter(){

				@Override
				public boolean accept(String arg0) {
					return true;
				}
				
			});
			outSummary.println(writeResults(ignoring));

			if (beds.size()>0){
				
				class MyLineFilter implements LineFilter{
					private int annotation;
					public MyLineFilter(int annotation) {
						this.annotation = annotation;
					}
					@Override
					public boolean accept(String arg0) {
						String[] parts=arg0.split("\t");
						if (!parts[annotation].equals("NULL")){
							return true;
						}else{
							return false;
						}
					}
				}
				
				int annotation = 7;
				for (IntervalsIndex index: beds){
					outSummary.println("=============================");
					outSummary.println("RESULTS WITH ANNOTATION: "+index.getName());
					outSummary.println("=============================");
		
					
					Results forAnnotation = getResults(methylcytosines, new MyLineFilter(annotation));
					outSummary.println(writeResults(forAnnotation));
					
					annotation ++;
				}
				
				outSummary.println("=============================");
				outSummary.println("RESULTS WITH ALL ANNOTATIONS");
				outSummary.println("=============================");
	
				Results allAnnotations = getResults(methylcytosines, new LineFilter(){
	
					@Override
					public boolean accept(String arg0) {
						String[] parts=arg0.split("\t",8);
						if (parts[7].indexOf("NULL")==-1){
							return true;
						}else{
							return false;
						}
					}
					
				});
				outSummary.println(writeResults(allAnnotations));
				annotation++;
			}
			
			outSummary.close();

			
		}

		
		private static File buildBAMAndIndex(File samCT, StringBuilder samtoolsDirectory) {
			File bam = new File(samCT.getAbsolutePath()+".bam");
			if (bam.exists() && bam.lastModified()>samCT.lastModified()){
			//	System.out.println("bam exists and is older. Skip");
				
			}else{
				System.out.println("Building BAM for "+samCT);
				
				Tools.executeProcessWait(samtoolsDirectory.toString()+"/samtools view -S -b -o "+bam.getAbsolutePath()+" "+samCT.getAbsolutePath());
				System.out.println("BAM built for "+samCT);
			}
			
			
			File bai = new File(bam.getAbsolutePath()+".bai");
			
			if (bai.exists() && bam.exists() && bai.lastModified()>bam.lastModified()){
			//	System.out.println("bai exists and is older. Skip");
				
			}else{
				System.out.println("Building index for "+samCT);
				
				Tools.executeProcessWait(samtoolsDirectory.toString()+"/samtools index "+bam.getAbsolutePath());
				System.out.println("index built for "+samCT);
			}
			
			
			return bam;
			
			
		}

		private static List<RecordAndOffset> filterClonalReads(
				List<RecordAndOffset> bases) {
			HashMap<Integer, RecordAndOffset> forward = new HashMap<Integer, RecordAndOffset>();
			HashMap<Integer, RecordAndOffset> reverse = new HashMap<Integer, RecordAndOffset>();
			for (RecordAndOffset record: bases){
				if((record.getRecord().getFlags() & 0x0010)==0x0010){
					if (!reverse.containsKey(record.getRecord().getAlignmentEnd()) || reverse.get(record.getRecord().getAlignmentEnd()).getRecord().getMappingQuality() < record.getRecord().getMappingQuality()){
						reverse.put(record.getRecord().getAlignmentEnd(), record);
						
					}
				}else{
					if (!forward.containsKey(record.getRecord().getAlignmentStart()) || forward.get(record.getRecord().getAlignmentStart()).getRecord().getMappingQuality() < record.getRecord().getMappingQuality()){
						forward.put(record.getRecord().getAlignmentStart(), record);
					}
				}
			}
			LinkedList<RecordAndOffset> toret = new LinkedList<RecordAndOffset>();
			
			for (RecordAndOffset record : forward.values()){
				
				toret.add(record);
			}
			for (RecordAndOffset record : reverse.values()){
				toret.add(record);
			}
		
			/*
			
			System.err.println("Filter---------");
			System.err.println("Before: ");
			for (RecordAndOffset record: bases){
				System.err.print(","+((record.getRecord().getFlags() & 0x0010)==0x0010?"-"+record.getRecord().getAlignmentStart():"+"+record.getRecord().getAlignmentEnd()));
			}
			System.err.println("\nAfter: ");
			for (RecordAndOffset record: toret){
				System.err.print(","+((record.getRecord().getFlags() & 0x0010)==0x0010?"-"+record.getRecord().getAlignmentStart():"+"+record.getRecord().getAlignmentEnd()));
			}
			System.err.println("\n------------------\n");
			*/
			return toret;
//			
		}

		/*public static void main(String[] args) throws IOException, InterruptedException{
			analyzeMethylation(new StringBuilder("/home/lipido/Desktop/methyldata/output/"), new StringBuilder("/home/lipido/Desktop/methyldata/output/"), new String[]{"c1.fa"}, new String[]{"borrar.txt"});
		}*/
		
		/**
		 * Counts the number of aligned and non-aligned reads per bowtie output file and saves them in a csv file 'readCounts.csv'
		 * @throws IOException 
		 * 
		 */
		public static void count(StringBuilder outputDirectory){
			
			String [] files=Tools.ls(outputDirectory.toString());
			StringBuilder outputToCsv=new StringBuilder("");
			System.out.println("Counting reads...... ");
			
			// por cada fichero del directorio output
			for(String file : files){
				// cuento sobre los ficheros SAM producidos por bowtie
				if(file.endsWith("aligned")){
					// hago las cuentas
					String readFileName=(file.split("[_][_]"))[0];
					String referenceFileName=(file.split("[_][_]"))[1];
					System.out.print(new StringBuilder("Counting [").append(readFileName).append("] vs [").append(referenceFileName).append("]...... "));
					int uniqueAligned=0;
					int nonAligned=0;
					outputToCsv.append("\tNon-aligned reads vs ").append(referenceFileName).append("\tUnique-aligned reads vs ").append(referenceFileName).append("\tTotal\n");
					
					BufferedReader br;
				
					String thisLine;

					try {
						br = new BufferedReader(new FileReader(new StringBuilder(outputDirectory).append(file).toString()));
						while ((thisLine = br.readLine())!=null){
							/*if(!thisLine.split("\t")[5].equals("*")) uniqueAligned++;
							else if(!thisLine.startsWith("@")) nonAligned++;*/
							if(!thisLine.startsWith("@")){
								if(!thisLine.split("\t")[5].equals("*")) uniqueAligned++;
								else nonAligned++;
							}
						}
						br.close();
						
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					
					// ahora escribo los resultados
					outputToCsv.append(readFileName).append("\t").append(nonAligned).append("\t").append(uniqueAligned).append("\t").append((uniqueAligned+nonAligned)).append("\n");
					System.out.println("[OK]");
				}//if(file.endsWith("aligned"))
			}//for(String file : files)
			
			//imprimo el fichero de salida con las cuentas
			BufferedWriter wr;
			try {
				String outputFile=new StringBuilder(outputDirectory).append("readCounts.csv").toString();
				if(Tools.doesThisFileExist(outputFile)) Tools.deleteFile(outputFile);
				wr = new BufferedWriter(new FileWriter(outputFile));
				wr.write(outputToCsv.toString());
				wr.flush();
				wr.close();				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			
		}
		
		/**
		 * Eliminates reads that align to both CT and GA references (ambiguous reads)
		 * @param bisulfitedRefs
		 * @param bisulfitedReads
		 * @param outputDirectory
		 */		
		public static void eliminateAmbiguousReads(String[] bisulfitedRefs,String [] bisulfitedReads,StringBuilder workingDirectory, StringBuilder outputDirectory, boolean hasL, int trimMismatches){
			
			System.out.println("Eliminating ambiguous reads...... ");
			
			// para cada fichero de read
			for(int i=0;i<bisulfitedReads.length;i++){
				String read=bisulfitedReads[i];

				// para cada fichero de muestra
				for(int j=0;j<bisulfitedRefs.length;j++){
					String ref=bisulfitedRefs[j];
					
					// uso unicamente la referencia CT. Desde esta obtengo la referencia GA directamente
					if(ref.indexOf("_CT_")!=-1){
						String [] aux=ref.split("_CT_");
						String GAref=new StringBuilder(aux[0]).append("_GA_").append(aux[1]).toString();
						
						String CT=new StringBuilder(outputDirectory).append(read).append("__").append(ref).append(".aligned").toString();
						String GA=new StringBuilder(outputDirectory).append(read).append("__").append(GAref).append(".aligned").toString();											
						
						
						System.out.println("Searching for ambiguous reads between\n"+CT+"\nand\n"+GA+"\n====");
						
	
											
						String sortedCT=new StringBuilder(outputDirectory).append("SORTED_CT").toString();
						String sortedGA=new StringBuilder(outputDirectory).append("SORTED_GA").toString();
						String pasted=new StringBuilder(outputDirectory).append("PASTED_").append(read).append("__CT_GA").toString();
						// si existe ya el fichero lo borro	
						if(Tools.doesThisFileExist(sortedCT)) Tools.deleteFile(sortedCT);
						if(Tools.doesThisFileExist(sortedGA)) Tools.deleteFile(sortedGA);
						
						System.out.println("1.- deleting headers in both files......");
						
						BufferedWriter wrCT = null;
						BufferedWriter wrGA = null;
						BufferedWriter ambiguous = null;
						
						// ordeno la salida de los ficheros por la primera columna						
						try {
						
							String[] command=new String[] {
									"sh",
									"-c",
									"grep -v '^@'  "+CT+" | sort -k 1 -o "+sortedCT+" -T "+workingDirectory.toString()
									//"grep -v '^@'  "+CT+" > "+sortedCT
							};
							Runtime.getRuntime().exec(command).waitFor();
							
							
							command=new String[] {
									"sh",
									"-c",
									"grep -v '^@'  "+GA+" | sort -k 1 -o "+sortedGA+" -T "+workingDirectory.toString()
									//"grep -v '^@'  "+GA+" > "+sortedGA
							};
							Runtime.getRuntime().exec(command).waitFor();		         	
							
							
							
							
							
							
							
							// ahora hago el paste
							System.out.println("2.- pasting files......");
							// tengo que decirle al paste que cuando pega las dos lineas, las separe usando el caracter 1 en ASCII. Este valor nunca va a interferir con los códigos ASCII de
							//calidad de las reads.
							//****Por defecto separaba con un \t y eso me daba problemas al ejecutar las picard para crear el fichero sam ordenado
							Character sepchar = 1;
							String sep = sepchar.toString();
				         	command=new String[]{
				         			"sh",
				         			"-c",
				         			"paste -d "+sep+" "+sortedCT+" "+sortedGA+" > "+pasted				         			
				         	};
				         	Runtime.getRuntime().exec(command).waitFor();
				         	
				         	
				         	//Process child=Runtime.getRuntime().exec(command);
				         	//FileOutputStream fileoutputstream=new FileOutputStream(new File(pasted));
				         	//Tools.executeProcessWait(command, fileoutputstream, null);
				         	//fileoutputstream.flush();
				         	//fileoutputstream.close();
				         	
				         	
				         	
				         	
				         	
				         	
				         	// recogidos los dos ficheros ordenados y pegados, selecciono  las reads que alinean en ambas referencias
				         	// ¿cómo lo hago?
				         	// leo el fichero PASTED y las lineas que NO tienen dos 36M las rompo, una a un fichero y otra a otro.
				         	// las que tienen dos 36M, es decir, alinearon en las dos referencias, las guardo como ambigous reads
				         	System.out.println("3.- trimming ambiguous reads......");
				         	
				         	BufferedReader br = new BufferedReader(new FileReader(pasted));
							String thisLine;
							// REUTILIZO los SORTED:
							//aprovecho los ficheros sorted para guardar temporalmente las reads que voy a utilizar (las que no desecho)
							// *** RAZON POR LA CUAL AL FINAL TIENEN MENOS READS QUE EL PASTED
							wrCT = new BufferedWriter(new FileWriter(sortedCT));
							wrGA = new BufferedWriter(new FileWriter(sortedGA));
							String ambiguousReads=new StringBuilder(outputDirectory).append(read).append("__CT_GA.ambiguous").toString();
							ambiguous = new BufferedWriter(new FileWriter(ambiguousReads));
							
							
							while ((thisLine = br.readLine()) != null){

								String bothReads[]=thisLine.split(sep);
								// separo las 2 reads a cadenas distintas
								String readCT=bothReads[0];
								String readGA=bothReads[1];
								//ahora rompo cada una de ellas
								String [] tokensReadCT=readCT.split("\t");
								String [] tokensReadGA=readGA.split("\t");
								
								// AQUI DESCARTO READS AMBIGUAS (las que alinean contra las 2 referencias)
								
								// rompo la linea por el separador ASCII que le indiqué a la hora de hacer el PASTE
								//compruebo si han alineado las 2
								if(!tokensReadCT[5].equals("*") && !tokensReadGA[5].equals("*")){
									// apunto que es una read ambigua porque alinea contra las dos referencias
									ambiguous.write(thisLine);
						        	ambiguous.newLine();									
								}
								else{// si no es ambigua, pongo cada una (cada read) en su fichero
									if (hasL){
										if (!tokensReadCT[5].equals("*")){
											//is aligned, could be trimmed
											String trimmed = trim(trimMismatches, tokensReadCT, readCT);
											
											wrCT.write(trimmed);											
											wrCT.newLine();
										}else{
											wrCT.write(readCT);
							        	 	wrCT.newLine();
										}
										if (!tokensReadGA[5].equals("*")){
											//is aligned, could be trimmed
											String trimmed = trim(trimMismatches, tokensReadGA, readGA);
											
											wrGA.write(trimmed);											
											wrGA.newLine();
										}else{
											wrGA.write(readGA);
							        	 	wrGA.newLine();
										}
									}else{
										wrCT.write(readCT);
						        	 	wrCT.newLine();
						        	 	
						        	 	wrGA.write(readGA);
						        	 	wrGA.newLine();
									}
								}
								
								/*String[] parts = thisLine.split("\t");
						        if(!parts[5].equals("*") && !parts[19].equals("*")){
						       	 	// bisulfito la secuencia
						        	ambiguous.write(thisLine);
						        	ambiguous.newLine();
								}
						        else{// si no es ambigua, rompo la linea y pongo cada una (cada read) en su fichero					        	
						        		
						        		// rompo por el separador que le indiqué a la hora de hacer el PASTE
						        		String []tokens=thisLine.split(sep);
						        		wrCT.write(tokens[0]);
						        	 	wrCT.newLine();

						        	 	wrGA.write(tokens[1]);
						        	 	wrGA.newLine();
						        	 	
						         }*/
						         
						         
						    } // end while
							
							// cierro el buffer de lectura
							br.close();

				         	
						}catch(Exception e){							
							e.printStackTrace();
						}
						finally{
							if(wrCT!=null){
								try {
									wrCT.flush();
									wrCT.close();
								} catch (IOException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
								
							}
							if(wrGA!=null){
								try {
									wrGA.flush();
									wrGA.close();
								} catch (IOException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
																
							}
							if(ambiguous!=null){
								try {
									ambiguous.flush();
									ambiguous.close();
								} catch (IOException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
															
							}
						}

						Tools.deleteFile(pasted);
						
						
						
						
						System.out.println("4.- restoring '.align' files......");
						// ahora devuelvo las reads no eliminadas a cada uno de los ficheros align originales
						try{
				         	// primero para el CT
							String tmp=new StringBuilder(outputDirectory).append("tmp").toString();
							String[] command=new String[] {
									"sh",
									"-c",
									"grep '^@'  "+CT+" > "+tmp
							};
							Runtime.getRuntime().exec(command).waitFor();
							command=new String[] {
									"sh",
									"-c",
									"cat "+tmp+" "+sortedCT+" > "+CT
							};
							Runtime.getRuntime().exec(command).waitFor();
							
							Tools.deleteFile(tmp);
							
						}catch(Exception e){
							e.printStackTrace();
						}
						
						
						try{
				         	// ahora para el GA
							String tmp=new StringBuilder(outputDirectory).append("tmp").toString();
							String[] command=new String[] {
									"sh",
									"-c",
									"grep '^@'  "+GA+" > "+tmp
							};
							Runtime.getRuntime().exec(command).waitFor();
							command=new String[] {
									"sh",
									"-c",
									"cat "+tmp+" "+sortedGA+" > "+GA
							};
							Runtime.getRuntime().exec(command).waitFor();
							
							Tools.deleteFile(tmp);
							
						}catch(Exception e){
							e.printStackTrace();
						}
						
						Tools.deleteFile(sortedCT);
						Tools.deleteFile(sortedGA);
						
						System.out.println("///////////");
						
						

					}//if(ref.indexOf("_CT_")!=-1)
					
				}//for(int j=0;j<refs.length;j++)
				
			}//for(int i=0;i<reads.length;i++)

		}

	
		private static String trim(int trimMismatches, String[] tokensRead, String original) {
			String sequenceCT = tokensRead[9];
			String NMString = tokensRead[13].trim();
			int mismatches = Integer.parseInt(NMString.substring(5));
			String trimmed = null;
			if (mismatches>=trimMismatches){
				String MDString = tokensRead[12].substring(5);
				String[] MDStringTokens = MDString.split("[ACTGN]");
				
				
				
				String newQuality=null;
				int alignmentPos = Integer.parseInt(tokensRead[3]);
				int flags = Integer.parseInt(tokensRead[1]);
				if((flags & 0x0010)==0x0010){
					//is in the reverse
					int newStart = sequenceCT.length();
					for (int m = MDStringTokens.length-1; m>=MDStringTokens.length-trimMismatches; m--){
    					newStart-=Integer.parseInt(MDStringTokens[m])+1;
    				}
					newStart++;
					sequenceCT = sequenceCT.substring(newStart);
					
					alignmentPos+=newStart;
					newQuality=tokensRead[10].substring(newStart);
					
        	 	}else{
        	 		int newLength=0;
    				
    				for (int m = 0; m<trimMismatches; m++){
    					newLength+=Integer.parseInt(MDStringTokens[m])+1;
    				}
    				newLength--; //ignore the lastmismatch
    				
    				sequenceCT = sequenceCT.substring(0,newLength);
    				newQuality=tokensRead[10].substring(0,newLength);
        	 	}
				
				
						
				
			
				trimmed = tokensRead[0]+"\t"+tokensRead[1]+"\t"+tokensRead[2]+"\t"
				+alignmentPos+"\t"+tokensRead[4]+"\t"+(sequenceCT.length())+"M"+"\t"
				+tokensRead[6]+"\t"+tokensRead[7]+"\t"+tokensRead[8]+"\t"
				+sequenceCT+"\t"+newQuality+"\t"+tokensRead[11]+"\t"
				+tokensRead[12];
				
			/*	if((flags & 0x0010)==0x0010){
					System.err.println("----------Trim in reverse----------");
					System.err.println("From: "+unTokenize(tokensRead,"\t"));
					System.err.println("To: "+trimmed);
					System.err.println("--------------------------------");
					System.err.println("");
					
				}*/
			}
			else{
				trimmed = original;
			}
			
			return trimmed;
		}


		
}

