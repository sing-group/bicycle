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
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Reader;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import es.cnio.bioinfo.bicycle.FastqSplitter;
import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Reference;
import es.cnio.bioinfo.bicycle.Sample;
import es.cnio.bioinfo.bicycle.Tools;

public class BowtieAlignment {
	public enum PairedEndOrientation{
		FR("--fr"), FF("--ff"), RF("--rf");
		private String bowtieParameter;
		private PairedEndOrientation(String bowtieParameter) {
			this.bowtieParameter=bowtieParameter;
		}
		public String getBowtieParameter() {
			return bowtieParameter;
		}
	}

	public enum Strand{
		WATSON, CRICK,
	}
	private static final Logger logger = Logger.getLogger(BowtieAlignment.class.getSimpleName());
	
	
	
	private Project project;
	public enum Quals {
		AFTER_1_3("--solexa1.3-quals"), BEFORE_1_3("--solexa-quals"), PHRED_33("--phred33-quals");
		
		private String parameterValue;
		
		private Quals(String parameterValue) {
			this.parameterValue = parameterValue;
		}
		public String getParameterValue() {
			return parameterValue;
		}
	}
	public BowtieAlignment(Project p) {
		this.project = p;
	}
	
	
	public void buildBowtieIndex(Reference reference) throws IOException{
		ReferenceBisulfitation rb = new ReferenceBisulfitation(this.project);
		
		
		
		File bisulfitedReferenceCT = rb.getBisulfitedReference(ReferenceBisulfitation.Replacement.CT, reference);
		File bisulfitedReferenceGA = rb.getBisulfitedReference(ReferenceBisulfitation.Replacement.GA, reference);

		if (!bisulfitedReferenceCT.exists()){
			throw new IllegalArgumentException("Reference CtoT bisulfitation "+bisulfitedReferenceCT+" could not be found. Please, perform reference in-silico bisulfitation first");
		}
		
		if (!bisulfitedReferenceGA.exists()){
			throw new IllegalArgumentException("Reference GtoA bisulfitation "+bisulfitedReferenceGA+" could not be found. Please, perform reference in-silico bisulfitation first");
		}
		buildBowtieIndex(bisulfitedReferenceCT);
		
		buildBowtieIndex(bisulfitedReferenceGA);
		
		
	}

	private void buildBowtieIndex(File bisulfitedReference) {
		logger.info("Building Bowtie index for "+bisulfitedReference.toString().replaceAll(project
				.getWorkingDirectory()+File.separator, Project.WORKING_DIRECTORY));

		String bowtieBuildPath = (this.project.getBowtieDirectory()!=null?this.project.getBowtieDirectory().getAbsolutePath()+File
				.separator:"")+"bowtie-build";

		String [] command=new String[]{bowtieBuildPath,bisulfitedReference.getAbsolutePath(),bisulfitedReference.getAbsolutePath()};
		
		int result=Tools.executeProcessWait(command,null,null);	
		
		if(result==0) logger.info("Bowtie index build OK");
		else{
			String commandString = "";
			
			for(int j=0;j<command.length;j++)
				commandString+=command[j]+" ";
			
			throw new RuntimeException("Error during bowtie index build. Command was: "+commandString);
		}
	}
	
	public File getAlignmentOutputFile(Strand strand, Sample s, Reference r){
		return new File(this.project.getOutputDirectory()+File.separator+"bisulfited_CT_"+s.getName()+"_against_"+r.getReferenceFile().getName()+"_"+strand.name()+".sam");
	}
	public void performBowtieAlignment(
			final Sample sample,			
			final Reference reference,
			boolean skipUnconverted,
			int threadsNumber,
			
			/* bowtie params */
			final int e,
			final int l,
			final int n,
			
			
			final int chunkmbs,
			final Quals solexaQ) throws IOException{
		performBowtieAlignment(sample, reference, skipUnconverted, threadsNumber, e, l, n, chunkmbs, solexaQ, 0, 250);
	}
	public void performBowtieAlignment(
			final Sample sample,			
			final Reference reference,
			boolean skipUnconverted,
			int threadsNumber,
			
			/* bowtie params */
			final int e,
			final int l,
			final int n,
			
			
			final int chunkmbs,
			final Quals solexaQ,
			
			/*bowtie paired-end parameters*/
			final int I, 
			final int X) throws IOException{

		logger.info("Peforming alignment of sample "+sample.getName()+" against "+reference.getReferenceFile()
				.toString().replaceAll(project.getReferenceDirectory()+File.separator, ""));
		final int M = 1; //tag if there are more than one possible alignment with XM:i:>2
		final int k = 1; //report only 1 alignment
		
		ReferenceBisulfitation rb = new ReferenceBisulfitation(this.project);
		
		final File refCT = rb.getBisulfitedReference(ReferenceBisulfitation.Replacement.CT, reference); 
		final File refGA = rb.getBisulfitedReference(ReferenceBisulfitation.Replacement.GA, reference);
		
		if (!refCT.exists()){
			throw new IllegalArgumentException("cannot find in-silico CtoT bifulfited reference: "+refCT+". Perform reference bisulfitation first.");
		}
		
		if (!refGA.exists()){
			throw new IllegalArgumentException("cannot find in-silico GtoA bifulfited reference: "+refGA+". Perform reference bisulfitation first.");
		}
		
		final File alignmentOutputFileCT = getAlignmentOutputFile(Strand.WATSON, sample, reference);
		final File alignmentOutputFileGA = getAlignmentOutputFile(Strand.CRICK, sample, reference);
		
		
		int threads = threadsNumber/2;
		if (threads == 0) threads = 1;
		
		List<BufferedReader> streamsWATSON = null;
		List<BufferedReader> streamsCRICK = null;
		
		
		if (!sample.isPaired()){
			if (!sample.isDirectional()) {
				//non-directional (cokus)
				streamsWATSON = new LinkedList<BufferedReader>();
				for (BufferedReader reader : FastqSplitter.splitfastq(sample.getReadsFiles(), threads)) {
					streamsWATSON.add(new GtoADuplicatorReader(new CtoTReader(sample, reader, skipUnconverted)));
				}
				streamsCRICK = new LinkedList<BufferedReader>();
				for (BufferedReader reader : FastqSplitter.splitfastq(sample.getReadsFiles(), threads)) {
					streamsCRICK.add(new GtoADuplicatorReader(new CtoTReader(sample, reader, skipUnconverted)));
				}
				
				
			} else {
				//directional (lister)
				streamsWATSON = new LinkedList<BufferedReader>();
				for (BufferedReader reader : FastqSplitter.splitfastq(sample.getReadsFiles(), threads)) {
					streamsWATSON.add(new CtoTReader(sample, reader, skipUnconverted));
				}
				streamsCRICK = new LinkedList<BufferedReader>();
				for (BufferedReader reader : FastqSplitter.splitfastq(sample.getReadsFiles(), threads)) {
					streamsCRICK.add(new CtoTReader(sample, reader, skipUnconverted));
				}
				
			}
		}else {
			//paired end
			
			//streamsWATSON
			streamsWATSON = new LinkedList<BufferedReader>();
			List<BufferedReader> mate1Readers = FastqSplitter.splitfastq(sample.getReadsMate1Files(), threads);
			List<BufferedReader> mate2Readers = FastqSplitter.splitfastq(sample.getReadsMate2Files(), threads);
			
			for (int i = 0; i<mate1Readers.size(); i++){
				streamsWATSON.add(new PairedEndBowtieReader(sample, mate1Readers.get(i), mate2Readers.get(i), sample.isDirectional(), skipUnconverted));
			}
			
			//streamsCRICK
			streamsCRICK = new LinkedList<BufferedReader>();
			mate1Readers = FastqSplitter.splitfastq(sample.getReadsMate1Files(), threads);
			mate2Readers = FastqSplitter.splitfastq(sample.getReadsMate2Files(), threads);
			
			for (int i = 0; i<mate1Readers.size(); i++){
				streamsCRICK.add(new PairedEndBowtieReader(sample, mate1Readers.get(i), mate2Readers.get(i), sample.isDirectional(), skipUnconverted));
			}
		}
		
		abstract class LineProcessor{
			public abstract void processLine(String line);
		}
		
		
		
		class AlignerThread extends Thread{
			File ref;
			LineProcessor out;
			private String logFileName;
			private BufferedReader readsStream;
			private boolean nohead;
			boolean shouldStop = false;
			private Strand strand;
			
			public AlignerThread(File ref, Strand strand, LineProcessor out, String logfileName, BufferedReader readsStream, boolean nohead){
				this.ref = ref;
				this.out = out;
				this.logFileName = logfileName;
				this.readsStream = readsStream;
				this.nohead = nohead;
				this.strand = strand;
			}
			public void run(){
				logger.info("Aligning "+
						sample.getReadsFiles().toString().replaceAll(project.getReadsDirectory().toString()+File.separator,"")
						+" " +
						"against " +
						"["+ref.toString().replaceAll(project.getWorkingDirectory().toString()+File.separator, "")+"]...... " +
						"(see .log file)...... ");
				
				String [] command = null;

				String bowtiePath = (sample.getProject().getBowtieDirectory()!=null?sample.getProject()
						.getBowtieDirectory().getAbsolutePath()+File
						.separator:"")+"bowtie";
				if (!sample.isPaired()){
					if (nohead){
						command=new String[]{bowtiePath, "-t","--chunkmbs",""+chunkmbs,"--mm",solexaQ.getParameterValue
								(),"-e",""+e,"-l",""+l,"-n",""+n,"-k",""+k,"-M",""+M,"-S","--sam-nohead","--sam-RG","ID:"+strand.name(),"--sam-RG","SM:"+strand.name(),"--best","--nomaqround",ref.getAbsolutePath(),"-"};
					}else{
						command=new String[]{bowtiePath, "-t","--chunkmbs",""+chunkmbs,solexaQ.getParameterValue(),"-e",""+e,"-l",""+l,"-n",""+n,"-k",""+k,"-M",""+M,"-S","--best","--sam-RG","ID:"+strand.name(),"--sam-RG","SM:"+strand.name(),"--nomaqround",ref.getAbsolutePath(),"-"};
					}
				}else{ //paired
					if (nohead){
						command=new String[]{bowtiePath, "-t","--chunkmbs",""+chunkmbs,"--mm",solexaQ.getParameterValue(),"-e",""+e,"-l",""+l,"-n",""+n,"-k",""+k,"-M",""+M,"-S","--sam-nohead","--sam-RG","ID:"+strand.name(),"--sam-RG","SM:"+strand.name(),"--best","--nomaqround",ref.getAbsolutePath(),"-I", ""+I, "-X",""+X, "--fr", "--12", "-"};
					}else{
						command=new String[]{bowtiePath, "-t","--chunkmbs",""+chunkmbs,solexaQ.getParameterValue(),"-e",""+e,"-l",""+l,"-n",""+n,"-k",""+k,"-M",""+M,"-S","--best","--sam-RG","ID:"+strand.name(),"--sam-RG","SM:"+strand.name(),"--nomaqround",ref.getAbsolutePath(),"-I", ""+I, "-X",""+X, "--fr", "--12", "-"};
					}
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
							logger.info("Start read feeding to alignment against "+ref.toString().replaceAll(project
									.getWorkingDirectory()+File.separator, Project.WORKING_DIRECTORY)+". Log file: " +
									""+logFileName
									.replaceAll
									(project.getOutputDirectory()+File.separator, Project.OUTPUT_DIRECTORY));
							try {
								while ((readsLine=readsStream.readLine())!=null && !shouldStop){
									ps.println(readsLine);
									
								}
								ps.flush();
								ps.close();
								
								logger.info("Finished read feeding to alignment against "+ref.toString().replaceAll(project
										.getWorkingDirectory()+File.separator, Project.WORKING_DIRECTORY));
							} catch (Exception e) {						
								throw new RuntimeException(e);
							}
						};
					}.start();
					
					try {
						while ((line=stdOut.readLine())!=null){		
							out.processLine(line);						
						}
					} catch (IOException e1) {						
						throw new RuntimeException(e1);
					}
					
					shouldStop = true; //bowtie sends a null output, so the input feed should stop

				} catch (FileNotFoundException e1) {
					throw new RuntimeException(e1);
				}
				
			}
		}
		final PrintStream outCT = new PrintStream(new FileOutputStream(alignmentOutputFileCT));
		final PrintStream outGA = new PrintStream(new FileOutputStream(alignmentOutputFileGA));
		
		class AlignerPostprocessor{
			private int id;
			private String CTLine=null;
			private String GALine=null;
			
			private int tagCount = 0;
			private int mergeCount = 0;
			
			StringBuilder outputBufferCT = new StringBuilder(100000);
			StringBuilder outputBufferGA = new StringBuilder(100000);
			
			public AlignerPostprocessor(int id) {
				this.id = id;
			}
			public void close(){
				flushBuffer();
				logger.info("Both alignments have finished. Ambigous reads: "+tagCount);
				
			}
			
			private String directionalPreviousLineCT;
			private String directionalPreviousLineGA;
			
			private int previousEditDistance = -1;
			
			private int previousEditDistanceMate1 = -1;
			private int previousEditDistanceMate2 = -1;
			private int nonDirectionalPreviousEditDistanceMate1 = -1;
			private String directionalPreviousLineCTMate1;
			private String directionalPreviousLineGAMate1;
			private String directionalPreviousLineCTMate2;
			private String directionalPreviousLineGAMate2;
			private String nonDirectionalPreviousLineCTMate1;
			private String nonDirectionalPreviousLineGAMate1;
			
			private Pattern editDistancePattern = Pattern.compile("\\tNM:i:([^\\n\\t]+)");
			public void merge(){
				
				CTLine = replaceOriginalRead(CTLine).trim();
				GALine = replaceOriginalRead(GALine).trim();
				
				boolean ambiguous = false;
				if (!CTLine.startsWith("@")){
					String[] tokensCT = CTLine.split("\t");
					String[] tokensGA = GALine.split("\t");
					
					//System.out.println(tokensCT[0]+" = "+tokensGA[0]);
					if (!sample.isPaired() && !tokensCT[0].equals(tokensGA[0])){
						// Note: this does not happen when bowtie says "Exhausted best-first chunk memory for read"

						logger.severe("BUG: reading two samrecords from CT and GA alignments with are a " +
										"different read	CT:"+CTLine+"\nGA:"+GALine); System.exit(1);
					} else if (sample.isPaired() && !tokensCT[0].substring(0, tokensCT[0].length()-1).equals
							(tokensGA[0].substring(0, tokensGA[0].length()-1))){
						logger.severe("BUG: reading two samrecords from CT and GA alignments with are a different " +
								"read (ignoring last character)\nCT:"+CTLine+"\nGA:"+GALine); System.exit(1);
					}
					
					if (!tokensCT[5].equals("*") && !tokensGA[5].equals("*")){
						//adding a flag to sam indicating that this is ambiguous. Also add RG (mandatory in gatk)
						ambiguous = true;
						tagCount++;
					}
					
					if (sample.isDirectional()) {
						if (ambiguous) {
							outputBufferCT.append(CTLine+"\tZA:A:Y\tRG:Z:"+Strand.WATSON.name()+"\n");
							outputBufferGA.append(GALine+"\tZA:A:Y\tRG:Z:"+Strand.CRICK.name()+"\n");
						} else {
							outputBufferCT.append(CTLine+"\tRG:Z:"+Strand.WATSON.name()+"\n");
							outputBufferGA.append(GALine+"\tRG:Z:"+Strand.CRICK.name()+"\n");
						}
					} else {
						//get edit distance
						int currentEditDistance = -1;
						if (!tokensCT[5].equals("*")) {
							Matcher matcher = editDistancePattern.matcher(CTLine);
							matcher.find();
							currentEditDistance = Integer.parseInt(matcher.group(1));
						} 
						
						
						
						if (!tokensGA[5].equals("*")) {
							Matcher matcher = editDistancePattern.matcher(GALine);
							matcher.find();
							int GAEditDistance = Integer.parseInt(matcher.group(1));
							if (GAEditDistance < currentEditDistance || currentEditDistance == -1){
								currentEditDistance = GAEditDistance;
							}
						}
						
						
						if (sample.isPaired()) {
							if (directionalPreviousLineCTMate1 == null) {
								if (ambiguous) {
									directionalPreviousLineCTMate1 = CTLine+"\tZA:A:Y\tRG:Z:"+Strand.WATSON.name()+"\n";
									directionalPreviousLineGAMate1 = GALine+"\tZA:A:Y\tRG:Z:"+Strand.CRICK.name()+"\n";
								} else {
									directionalPreviousLineCTMate1 = CTLine+"\tRG:Z:"+Strand.WATSON.name()+"\n";
									directionalPreviousLineGAMate1 = GALine+"\tRG:Z:"+Strand.CRICK.name()+"\n";
								}
								previousEditDistanceMate1 = currentEditDistance;
							} else if (directionalPreviousLineCTMate2 == null) {
								if (ambiguous) {
									directionalPreviousLineCTMate2 = CTLine+"\tZA:A:Y\tRG:Z:"+Strand.WATSON.name()+"\n";
									directionalPreviousLineGAMate2 = GALine+"\tZA:A:Y\tRG:Z:"+Strand.CRICK.name()+"\n";
								} else {
									directionalPreviousLineCTMate2 = CTLine+"\tRG:Z:"+Strand.WATSON.name()+"\n";
									directionalPreviousLineGAMate2 = GALine+"\tRG:Z:"+Strand.CRICK.name()+"\n";
								}
								previousEditDistanceMate2 = currentEditDistance;
							} else if (nonDirectionalPreviousLineCTMate1 == null) {
								if (ambiguous) {
									nonDirectionalPreviousLineCTMate1 = CTLine+"\tZA:A:Y\tRG:Z:"+Strand.WATSON.name()+"\n";
									nonDirectionalPreviousLineGAMate1 = GALine+"\tZA:A:Y\tRG:Z:"+Strand.CRICK.name()+"\n";
								} else {
									nonDirectionalPreviousLineCTMate1 = CTLine+"\tRG:Z:"+Strand.WATSON.name()+"\n";
									nonDirectionalPreviousLineGAMate1 = GALine+"\tRG:Z:"+Strand.CRICK.name()+"\n";
								}
								nonDirectionalPreviousEditDistanceMate1 = currentEditDistance;
							} else {
								//which mate to print?? those with better edit distance
								int directionalEditDistance = Math.max(previousEditDistanceMate1, previousEditDistanceMate2);
								int nonDirectionalEditDistance = Math.max(nonDirectionalPreviousEditDistanceMate1, currentEditDistance);
								
								if (directionalEditDistance >= nonDirectionalEditDistance) {
									outputBufferCT.append(directionalPreviousLineCTMate1);
									outputBufferCT.append(directionalPreviousLineCTMate2);
									outputBufferGA.append(directionalPreviousLineGAMate1);
									outputBufferGA.append(directionalPreviousLineGAMate2);
								} else {
									//we have better edit distance
									outputBufferCT.append(nonDirectionalPreviousLineCTMate1);										
									outputBufferGA.append(nonDirectionalPreviousLineGAMate1);										
									if (ambiguous) {
										outputBufferCT.append(CTLine+"\tZA:A:Y\tRG:Z:"+Strand.WATSON.name()+"\n");
										outputBufferGA.append(GALine+"\tZA:A:Y\tRG:Z:"+Strand.CRICK.name()+"\n");
									} else {
										outputBufferCT.append(CTLine+"\tRG:Z:"+Strand.WATSON.name()+"\n");
										outputBufferGA.append(GALine+"\tRG:Z:"+Strand.CRICK.name()+"\n");
									}
									
									
								}
								previousEditDistanceMate1 = -1;
								previousEditDistanceMate2 = -1;
								nonDirectionalPreviousEditDistanceMate1 = -1;
								directionalPreviousLineCTMate1 = null;
								directionalPreviousLineGAMate1 = null;
								directionalPreviousLineCTMate2 = null;
								directionalPreviousLineGAMate2 = null;
								nonDirectionalPreviousLineCTMate1 = null;
								nonDirectionalPreviousLineGAMate1 = null;
							}
							
						} else {
							if (directionalPreviousLineCT == null) {
								
								if (ambiguous) {
									directionalPreviousLineCT = CTLine+"\tZA:A:Y\tRG:Z:"+Strand.WATSON.name()+"\n";
									directionalPreviousLineGA = GALine+"\tZA:A:Y\tRG:Z:"+Strand.CRICK.name()+"\n";
								} else {
									directionalPreviousLineCT = CTLine+"\tRG:Z:"+Strand.WATSON.name()+"\n";
									directionalPreviousLineGA = GALine+"\tRG:Z:"+Strand.CRICK.name()+"\n";
								}
								previousEditDistance = currentEditDistance;
							} else {
								// which of the two pairs should be written, the one with the best alignment. If none has alignment
								// keep the first one in the output
								if (currentEditDistance == -1) {
									// no matches in the non-directional reads, so do not output this line, only the preivous
									outputBufferCT.append(directionalPreviousLineCT);
									outputBufferGA.append(directionalPreviousLineGA);
								} else {
									if ((currentEditDistance < previousEditDistance && previousEditDistance != -1) ||
										previousEditDistance == -1 && currentEditDistance >=0 ) {
										// we have better edit distance, ignore previous, print the non-directional instead
										if (ambiguous) {
											outputBufferCT.append(CTLine+"\tZA:A:Y\tRG:Z:"+Strand.WATSON.name()+"\n");
											outputBufferGA.append(GALine+"\tZA:A:Y\tRG:Z:"+Strand.CRICK.name()+"\n");
										} else {
											outputBufferCT.append(CTLine+"\tRG:Z:"+Strand.WATSON.name()+"\n");
											outputBufferGA.append(GALine+"\tRG:Z:"+Strand.CRICK.name()+"\n");
										}
									} else {
										// no better edit distance, print only the previous
										outputBufferCT.append(directionalPreviousLineCT);
										outputBufferGA.append(directionalPreviousLineGA);
									}
								}
								
								directionalPreviousLineCT = null;
								directionalPreviousLineGA = null;
								previousEditDistance = -1;
								
							}
						}
					}
					
					
				} else {
					outputBufferCT.append(CTLine+"\n");
					outputBufferGA.append(GALine+"\n");
				}

				
				
				
				
				mergeCount++;
				
				if (mergeCount % 10000==0){
					flushBuffer();
					logger.info("Aligner thread "+this.id+": "+mergeCount+" reads processed");
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
			
			public LineProcessor CTProcessor = new LineProcessor(){

				@Override
				public void processLine(String line) {
					synchronized(AlignerPostprocessor.this){
						while(CTLine!=null){
							try {
								AlignerPostprocessor.this.wait(10000);
								if (CTLine!=null){
									logger.info("awaking, but CTLine is still not null (if you see this message continously, bowtie may be not responding), it is: "+CTLine+"\nProcessing new CT line: "+line);
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
									logger.info("awaking, but GALine is still not null, it is: "+GALine+"\nProcessing new GA line: "+line);
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
		AlignerPostprocessor dummyposprocessor = new AlignerPostprocessor(0);
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
		 
		for (int i = 0; i<streamsWATSON.size(); i++){
			AlignerPostprocessor postprocessor = new AlignerPostprocessor(i+1);
			postprocessors.add(postprocessor);
			LineProcessor ctProcessor = postprocessor.CTProcessor;
			LineProcessor gaProcessor = postprocessor.GAProcessor;
			BufferedReader streamCT = streamsWATSON.get(i);
			BufferedReader streamGA = streamsCRICK.get(i);
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

		logger.info("Alignment of sample "+sample.getName()+" OK");
		
	}
	
	private String getReverseComplementary(String sequence){
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
	private String replaceOriginalRead(String samline){
		if(!samline.startsWith("@")){
			//System.out.println(samline);
    	 	final String[] tokens=samline.split("[\t]");
    	 	// recupero la secuencia inicial
    	 	final String [] firstColumn=tokens[0].split("[|][|]");
    	 	String originalRead=null;
    	 	
    	 	int flag = Integer.parseInt(tokens[1]);
    	 	boolean mate1 = false;
    	 	boolean paired = false;
    	 	if ((flag & 0x0001)==0x0001){
    	 		//paired!
    	 		paired=true;
    	 		if((flag & 0x0040)==0x0040){
    	 			//mate1
    	 			mate1=true;
    	 			originalRead = firstColumn[1];
    	 		}else if((flag & 0x0080)==0x0080){
    	 			//mate2
    	 			mate1=false;
    	 			originalRead = firstColumn[2];
    	 		}else{
    	 			throw new RuntimeException("Malformed FLAG in SAM. It says that is a paired read, but it is not the first nor the second pair");
    	 		}
    	 	}else{    	 		
    	 		originalRead = firstColumn[1];
    	 	}
    	 	final StringBuilder lineModified=new StringBuilder();
    	 	
    	 	if((flag & 0x0010)==0x0010){
    	 		originalRead=getReverseComplementary(originalRead);
    	 	}
    	 	
    	 	
    	 	for(int j=0;j<tokens.length;j++){
    	 		// le adjunto la cabecera original de la read
    	 		if(j==0){
    	 			lineModified.append(firstColumn[0]);
    	 			if (paired){
    	 				if (mate1){
    	 					lineModified.append("/1");
    	 				}else{
    	 					lineModified.append("/2");
    	 				}
    	 			}
    	 			lineModified.append("\t");
    	 		}
    	 		// le adjunto la read original en lugar de la que tenia
    	 		else if(j==9) lineModified.append(originalRead).append("\t");
    	 		
    	 		//else if(j==tokens.length-1) lineModified=new StringBuilder(lineModified).append(tokens[j]);
    	 		else lineModified.append(tokens[j]).append("\t");
    	 		
    	 	}
    	 	return lineModified.toString();
		}//if(!thisLine.startsWith("@"))
		return samline;
	}
	
	
}


