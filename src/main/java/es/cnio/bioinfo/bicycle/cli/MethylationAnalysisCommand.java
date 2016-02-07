/*

Copyright 2012 Daniel Gonzalez Pe��a, Osvaldo Gra��a


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

package es.cnio.bioinfo.bicycle.cli;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.util.Hashtable;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import es.cnio.bioinfo.bicycle.ErrorRateMode;
import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Reference;
import es.cnio.bioinfo.bicycle.Sample;
import es.cnio.bioinfo.bicycle.operations.MethylationAnalysis;

public class MethylationAnalysisCommand extends ProjectCommand {

	@Override
	public String getName() {
		return "analyze-methylation";
	}

	@Override
	public String getDescription() {
		return "Analyzes methylation levels over the Sam files with the GATK-based walker";
	}

	@Override
	public void executeImpl(CLIApplication app, Project project, Map<Option, String> parameters) throws Exception {
		
		
		int nThreads = Integer.parseInt(parameters.get(this.findOption("n")));
		boolean removeBad = parameters.containsKey(this.findOption("r"));
		boolean removeAmbiguous = parameters.containsKey(this.findOption("a"));		
		boolean removeClonal = parameters.containsKey(this.findOption("c"));
		boolean correctNonCG = parameters.containsKey(this.findOption("g"));
		int trimuntil = Integer.parseInt(parameters.get(this.findOption("t")));
		boolean trimreads = trimuntil!=-1;
		int mindepth = Integer.parseInt(parameters.get(this.findOption("d")));
		
		double fdr = Double.parseDouble(parameters.get(this.findOption("f")));
		
		List<File> bedFiles = new LinkedList<File>();
		
		if (parameters.containsKey(this.findOption("b"))){
			String bedFilesString = parameters.get(this.findOption("b"));
			if (bedFilesString!=null){
				String[] tokens = bedFilesString.split(",");
				
				for (String token: tokens){
					File bedFile = new File(token);
					if (!bedFile.exists()){
						throw new IllegalArgumentException("BED file not found: "+bedFile);
					}
					bedFiles.add(bedFile);
				}
			}
		}
		String errorString = parameters.get(this.findOption("e"));
		String[] errorTokens = errorString.split("[=]");
		String errorModeString = errorTokens[0];
		
		ErrorRateMode errorMode = ErrorRateMode.valueOf(errorModeString); //may illegal argument exception
		
		
		MethylationAnalysis ma = new MethylationAnalysis(project);
		if (errorMode == ErrorRateMode.from_control_genome){
			
			if (errorTokens.length<2){
				throw new IllegalArgumentException("control genome must be set");
			}
			
			String controlGenome = errorTokens[1];
			
			for (Sample sample: project.getSamples()){
				for (Reference reference: project.getReferences()){					
					ma.analyzeWithErrorFromControlGenome(reference, sample, trimreads, trimuntil, removeAmbiguous, removeBad, removeClonal, correctNonCG, mindepth, fdr, nThreads, bedFiles, controlGenome);
				}
			}
			
		}else if(errorMode == ErrorRateMode.FIXED){
			if (errorTokens.length<2){
				throw new IllegalArgumentException("error rates must be set");
			}
			String[] errorRates = errorTokens[1].split(",");
			if (errorRates.length!=2){
				throw new IllegalArgumentException("bad error rates. It must be <watson_rate>,<crick_rate>. E.g.: 0.01,0.02");
			}
			
			double watsonError = Double.parseDouble(errorRates[0]);
			double crickError = Double.parseDouble(errorRates[1]);

			for (Sample sample: project.getSamples()){
				for (Reference reference: project.getReferences()){					
					ma.analyzeWithFixedErrorRate(reference, sample, trimreads, trimuntil, removeAmbiguous, removeBad, removeClonal, correctNonCG, mindepth, fdr, nThreads, bedFiles, watsonError, crickError);
				}
			}
		}else if (errorMode == ErrorRateMode.from_barcodes){
			for (Sample sample: project.getSamples()){
				for (Reference reference: project.getReferences()){					
					ma.analyzeWithErrorFromBarcodes(reference, sample, trimreads, trimuntil, removeAmbiguous, removeBad, removeClonal, correctNonCG, mindepth, fdr, nThreads, bedFiles);
				}
			}
			
		}
		
		
		// <---- Cytosine METHYLATION PER ANNOTATED REGION ---->
		//IF BED FILES ARE AVAILABLE, then cytosine methylation per annotated region is also calculated
		//added by osvaldo, 23Jan2016
		if(!bedFiles.isEmpty()){
			List<File> methylcytosinesFile = new LinkedList<File>();
			
			//recovers all the methylcytosines files
			for (Sample sample: project.getSamples()){
				for (Reference reference: project.getReferences()){					
					methylcytosinesFile.add(new File(project.getOutputDirectory()+File.separator+sample.getName()+"_"+reference.getReferenceFile().getName()+".methylcytosines"));
				}
			}
			
			//calculates methylation for each annotated region in each bed file for each methylcytosines file
			
			//for each methylcytosines file
			for(File elem:methylcytosinesFile){
				System.out.println("[ Calculating methylation per annotated region for: "+elem.toString()+" ]");
				
				//foreach bed file
				for(File bed:bedFiles){
					System.out.println("\t\t[ Checking methylation for regions annotated in "+bed.toString()+" ]");
					BufferedWriter w=null;
					
					try{ 
						//output file that will store methylation values for each annotated region
						String outputFile=elem.toString().replace("methylcytosines","")+bed.getName().replace("bed", "METHYLATEDregions.txt");
						
						FileWriter o=new FileWriter(outputFile);
						w=new BufferedWriter(o);
						w.write("Region\tmCG WATSON\tdepthCG WATSON\tCG WATSON methylation\tmCG CRICK\tdepthCG CRICK\tCG CRICK methylation");
						w.write("\tmCG total methylation\tmCG total depth\tCG TOTAL methylation");
						w.write("\tmCHG WATSON\tdepthCHG WATSON\tCHG WATSON methylation\tmCHG CRICK\tdepthCHG CRICK\tCHG CRICK methylation");
						w.write("\tmCHG total methylation\tmCHG total depth\tCHG TOTAL methylation");
						w.write("\tmCHH WATSON\tdepthCHH WATSON\tCHH WATSON methylation\tmCHH CRICK\tdepthCHH CRICK\tCHH CRICK methylation");
						w.write("\tmCHH total methylation\tmCHH total depth\tCHH TOTAL methylation");
						w.newLine();
						w.flush();
						
						//opens the methylcytosine file
						FileReader f=new FileReader(elem);
						BufferedReader b=new BufferedReader(f);
						
						boolean firstLine=true;
						int columnOfInterest=-1;
						String previousGenomicRegion="";
						boolean firstCheck=true;
						
						//each hash table contains separated values for Watson and Crick
						Hashtable<String, Integer> mCG=new Hashtable<String, Integer>();
						Hashtable<String, Integer> mCHG=new Hashtable<String, Integer>();
						Hashtable<String, Integer> mCHH=new Hashtable<String, Integer>();
						Hashtable<String, Integer> depthCG=new Hashtable<String, Integer>();
						Hashtable<String, Integer> depthCHG=new Hashtable<String, Integer>();
						Hashtable<String, Integer> depthCHH=new Hashtable<String, Integer>();
						
						//initialization
						mCG.put("WATSON", 0);
						mCG.put("CRICK", 0);
						mCHG.put("WATSON", 0);
						mCHG.put("CRICK", 0);
						mCHH.put("WATSON", 0);
						mCHH.put("CRICK", 0);
						depthCG.put("WATSON", 0);
						depthCG.put("CRICK", 0);
						depthCHG.put("WATSON", 0);
						depthCHG.put("CRICK", 0);
						depthCHH.put("WATSON", 0);
						depthCHH.put("CRICK", 0);
						
						String line;
						
						//for each line in the methylcytosines file (starting from line 0, i.e., first line)
						while((line=b.readLine())!=null){							
							String tokens[]=line.split("\t");
							
							//if positioned in first (header) line (line 0), it finds out the column number that contains
							//the genomic annotations for the current bed file
							if(firstLine){
								//checks all first line headers
								for(int pos=0;pos<tokens.length;pos++){
									if(tokens[pos].equals(bed.getName())){
										//the current bed file is in column 'pos'
										columnOfInterest=pos;
										break;
									}
								}
								
								firstLine=false;								
								
							}//if(firstLine)
							else{// not in first line (not header line)
																
								if(!tokens[columnOfInterest].contains("N/A")){
								
									if(firstCheck){
										//initialization
										previousGenomicRegion=tokens[columnOfInterest];
										firstCheck=false;
									}
									
									//value accumulation
									if(tokens[columnOfInterest].equals(previousGenomicRegion)){
										//accumulates methylated cytosines
										String strand=tokens[2];
										String methylationContext=tokens[3];
										int depth=Integer.parseInt(tokens[4]);
										int methylation=Integer.parseInt(tokens[6]);
										
										//accumulates new values for the current annotated region
										switch(methylationContext){
										case "CG":
											mCG.put(strand, mCG.get(strand)+methylation);
											depthCG.put(strand, depthCG.get(strand)+depth);
											break;
											
										case "CHG":
											mCHG.put(strand, mCHG.get(strand)+methylation);
											depthCHG.put(strand, depthCHG.get(strand)+depth);
											break;
											
										case "CHH":
											mCHH.put(strand, mCHH.get(strand)+methylation);
											depthCHH.put(strand, depthCHH.get(strand)+depth);
											break;
											
										default:
											System.out.print("[Invalid methylation context]: ");
											System.out.println(line);
										}
										
										
									}else{
										// CYTOSINE COUNTS
										//counts cytosine methylation for the checked region
										System.out.println("\t\t\t[obtaining methylation for] "+tokens[columnOfInterest]);
										
										//mCG
										double result=(double)mCG.get("WATSON")/(double)depthCG.get("WATSON");										
										DecimalFormat df = new DecimalFormat("#.#######");										
										StringBuffer buf=new StringBuffer(tokens[columnOfInterest]);
										buf.append("\t").append(mCG.get("WATSON")).append("\t").append(depthCG.get("WATSON")).append("\t").append(df.format(result));
										result=(double)mCG.get("CRICK")/(double)depthCG.get("CRICK");
										buf.append("\t").append(mCG.get("CRICK")).append("\t").append(depthCG.get("CRICK")).append("\t").append(df.format(result));
										result=(double)(mCG.get("WATSON")+mCG.get("CRICK"))/(double)(depthCG.get("WATSON")+depthCG.get("CRICK"));
										buf.append("\t").append(mCG.get("WATSON")+mCG.get("CRICK")).append("\t").append(depthCG.get("WATSON")+depthCG.get("CRICK")).append("\t").append(df.format(result));
										
										//mCHG
										result=(double)mCHG.get("WATSON")/(double)depthCHG.get("WATSON");
										buf.append("\t").append(mCHG.get("WATSON")).append("\t").append(depthCHG.get("WATSON")).append("\t").append(df.format(result));
										result=(double)mCHG.get("CRICK")/(double)depthCHG.get("CRICK");
										buf.append("\t").append(mCHG.get("CRICK")).append("\t").append(depthCHG.get("CRICK")).append("\t").append(df.format(result));
										result=(double)(mCHG.get("WATSON")+mCHG.get("CRICK"))/(double)(depthCHG.get("WATSON")+depthCHG.get("CRICK"));
										buf.append("\t").append(mCHG.get("WATSON")+mCHG.get("CRICK")).append("\t").append(depthCHG.get("WATSON")+depthCHG.get("CRICK")).append("\t").append(df.format(result));
										
										//mCHH
										result=(double)mCHH.get("WATSON")/(double)depthCHH.get("WATSON");
										buf.append("\t").append(mCHH.get("WATSON")).append("\t").append(depthCHH.get("WATSON")).append("\t").append(df.format(result));
										result=(double)mCHH.get("CRICK")/(double)depthCHH.get("CRICK");
										buf.append("\t").append(mCHH.get("CRICK")).append("\t").append(depthCHH.get("CRICK")).append("\t").append(df.format(result));
										result=(double)(mCHH.get("WATSON")+mCHH.get("CRICK"))/(double)(depthCHH.get("WATSON")+depthCHH.get("CRICK"));
										buf.append("\t").append(mCHH.get("WATSON")+mCHH.get("CRICK")).append("\t").append(depthCHH.get("WATSON")+depthCHH.get("CRICK")).append("\t").append(df.format(result));
										
										w.write(buf.toString());
										w.newLine();
										
										//moves to the next annotated region
										previousGenomicRegion=tokens[columnOfInterest];
										
										//introduces the first methylation value for the next annotated region
										String strand=tokens[2];
										String methylationContext=tokens[3];
										int depth=Integer.parseInt(tokens[4]);
										int methylation=Integer.parseInt(tokens[6]);
										
										
										switch(methylationContext){
										case "CG":
											mCG.put(strand, methylation);
											depthCG.put(strand, depth);
											break;
											
										case "CHG":
											mCHG.put(strand, methylation);
											depthCHG.put(strand, depth);
											break;
											
										case "CHH":
											mCHH.put(strand, methylation);
											depthCHH.put(strand, depth);
											break;
											
										default:
											System.out.print("[Invalid methylation context]: ");
											System.out.println(line);
										}//switch
									
									}//else
									
								}//if(!tokens[columnOfInterest].equals("N/A")
								
								

							}//else{// not in first line (not header line)
							
						}//while((line=b.readLine())!=null)
						
						b.close();
						
					
					}catch(FileNotFoundException e){
						e.printStackTrace();
						
					}catch(IOException e){
						e.printStackTrace();
						
					}finally{
						if ( w != null ) w.close();
					}
					
					
					
				}
				
			}
		}//if(!bedFiles.isEmpty())
	}

	@Override
	protected List<Option> createOptions() {
		List<Option> toret = super.createOptions();
		
		toret.add(new DefaultValuedOption("threads", "n", 
				"number of threads to analyze",  "4"));
		
		toret.add(new Option("remove-uncorrectly-converted", "r", 
				"ignore non-correctly bisulfite-converted reads", true, false));
		
		toret.add(new Option("remove-ambiguous", "a", 
				"ignore reads aligned to both Watson and Crick strands", true, false));
		
		toret.add(new DefaultValuedOption("trim-reads", "t", 
				"Trim reads to the <t> mismatch. -1 means no trim",  "4"));
		
		toret.add(new DefaultValuedOption("min-depth", "d", 
				"Ignore positions with less than <d> reads",  "1"));
		
		toret.add(new DefaultValuedOption("fdr", "f", 
				"FDR threshold", "0.01"));
		
		toret.add(new DefaultValuedOption("error-mode", "e", 
				"Error rate computation mode. Valid options are: "+ErrorRateMode.from_control_genome+"=<control_genome_name>, "+ErrorRateMode.from_barcodes+", "+ErrorRateMode.FIXED.name()+"=<watson_error_rate,crick_error_rate>", ErrorRateMode.FIXED.name()+"=0.01,0.01"));
		
		toret.add(new Option("annotate-beds", "b", 
				"Comma-separated (with no spaces) list of BED files to annotate cytosines", true, true));
		
		toret.add(new Option("remove-clonal", "c", 
				"Remove clonal reads", true, false));
		
		toret.add(new Option("correct non-CG to CG", "g", 
				"Correct non-CG", true, false));
		
		return toret;
	}

}
