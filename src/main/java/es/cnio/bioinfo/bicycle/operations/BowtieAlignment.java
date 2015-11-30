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
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Reader;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Logger;

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
	private static final Logger logger = Logger.getLogger(BowtieAlignment.class.getName());
	
	
	
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
		logger.info("Bowtie: indexing "+bisulfitedReference+"...... ");
		String [] command=new String[]{this.project.getBowtieDirectory().getAbsolutePath()+File.separator+("bowtie-build"),bisulfitedReference.getAbsolutePath(),bisulfitedReference.getAbsolutePath()};
		
		int result=Tools.executeProcessWait(command,null,null);	
		
		if(result==0) logger.info("[OK]");
		else{
			String commandString = "";
			
			for(int j=0;j<command.length;j++)
				commandString+=command[j]+" ";
			
			logger.severe("\nError: "+commandString);
			throw new RuntimeException("Error during bowtie index build. Command was: "+commandString);
		}
	}
	
	public File getAlignmentOutputFile(Strand strand, Sample s, Reference r){
		return new File(this.project.getOutputDirectory()+File.separator+"bisulfited_CT_"+s.getName()+"_against_"+r.getReferenceFile().getName()+"_"+strand.name()+".sam");
	}
	public void performBowtieAlignment(
			final Sample sample,			
			final Reference reference,			
			int threadsNumber,
			
			/* bowtie params */
			final int e,
			final int l,
			final int n,
			
			
			final int chunkmbs,
			final Quals solexaQ) throws IOException{
		performBowtieAlignment(sample, reference, threadsNumber, e, l, n, chunkmbs, solexaQ, 0, 250);
	}
	public void performBowtieAlignment(
			final Sample sample,			
			final Reference reference,			
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
		
		SampleBisulfitation sb = new SampleBisulfitation(sample);		
		
		List<BufferedReader> streamsWATSON = null;
		List<BufferedReader> streamsCRICK = null;
		
		
		if (!sample.isPaired()){
			final List<File> readsFiles = new LinkedList<File>();
			for(File readFile : sample.getReadsFiles()){
				File bisulfitedRead = sb.getBisulfitedFile(readFile);
				if (!bisulfitedRead.exists()){
					throw new IllegalArgumentException("cannot find in-silico bifulfited read file: "+bisulfitedRead+". Perform read bisulfitation firsrt");
				}
				readsFiles.add(bisulfitedRead);
			}
			streamsWATSON = FastqSplitter.splitfastq(readsFiles, threads);
			streamsCRICK = FastqSplitter.splitfastq(readsFiles, threads);
		}else {
			//paired end
			final List<File> readsFilesM1 = new LinkedList<File>();
			final List<File> readsFilesM2 = new LinkedList<File>();
			for(File readFile : sample.getReadsMate1Files()){
				File bisulfitedRead = sb.getBisulfitedFile(readFile);
				if (!bisulfitedRead.exists()){
					throw new IllegalArgumentException("cannot find in-silico bifulfited read file: "+bisulfitedRead+". Perform read bisulfitation firsrt");
				}
				readsFilesM1.add(bisulfitedRead);
			}
			for(File readFile : sample.getReadsMate2Files()){
				File bisulfitedRead = sb.getBisulfitedFile(readFile);
				if (!bisulfitedRead.exists()){
					throw new IllegalArgumentException("cannot find in-silico bifulfited read file: "+bisulfitedRead+". Perform read bisulfitation firsrt");
				}
				readsFilesM2.add(bisulfitedRead);
			}
			
			//streamsWATSON
			streamsWATSON = new LinkedList<BufferedReader>();
			List<BufferedReader> mate1Readers = FastqSplitter.splitfastq(readsFilesM1, threads);
			List<BufferedReader> mate2Readers = FastqSplitter.splitfastq(readsFilesM2, threads);
			
			for (int i = 0; i<mate1Readers.size(); i++){
				streamsWATSON.add(new PairedEndBowtieReader(mate1Readers.get(i), mate2Readers.get(i), Strand.WATSON));
			}
			
			//streamsCRICK
			streamsCRICK = new LinkedList<BufferedReader>();
			mate1Readers = FastqSplitter.splitfastq(readsFilesM1, threads);
			mate2Readers = FastqSplitter.splitfastq(readsFilesM2, threads);
			
			for (int i = 0; i<mate1Readers.size(); i++){
				streamsCRICK.add(new PairedEndBowtieReader(mate1Readers.get(i), mate2Readers.get(i), Strand.CRICK));
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
				
				
				logger.info("Bowtie: aligning "+sample.getReadsFiles()+" against ["+ref+"]...... (see .log file)...... ");
				
				String [] command = null;
				
				if (!sample.isPaired()){
					if (nohead){
						command=new String[]{sample.getProject().getBowtieDirectory().getAbsolutePath()+File.separator+"bowtie", "-t","--chunkmbs",""+chunkmbs,"--mm",solexaQ.getParameterValue(),"-e",""+e,"-l",""+l,"-n",""+n,"-k",""+k,"-M",""+M,"-S","--sam-nohead","--sam-RG","ID:"+strand.name(),"--sam-RG","SM:"+strand.name(),"--best","--nomaqround",ref.getAbsolutePath(),"-"};
					}else{
						command=new String[]{sample.getProject().getBowtieDirectory().getAbsolutePath()+File.separator+"bowtie", "-t","--chunkmbs",""+chunkmbs,solexaQ.getParameterValue(),"-e",""+e,"-l",""+l,"-n",""+n,"-k",""+k,"-M",""+M,"-S","--best","--sam-RG","ID:"+strand.name(),"--sam-RG","SM:"+strand.name(),"--nomaqround",ref.getAbsolutePath(),"-"};				
					}
				}else{ //paired
					if (nohead){
						command=new String[]{sample.getProject().getBowtieDirectory().getAbsolutePath()+File.separator+"bowtie", "-t","--chunkmbs",""+chunkmbs,"--mm",solexaQ.getParameterValue(),"-e",""+e,"-l",""+l,"-n",""+n,"-k",""+k,"-M",""+M,"-S","--sam-nohead","--sam-RG","ID:"+strand.name(),"--sam-RG","SM:"+strand.name(),"--best","--nomaqround",ref.getAbsolutePath(),"-I", ""+I, "-X",""+X, "--fr", "--12", "-"};
					}else{
						command=new String[]{sample.getProject().getBowtieDirectory().getAbsolutePath()+File.separator+"bowtie", "-t","--chunkmbs",""+chunkmbs,solexaQ.getParameterValue(),"-e",""+e,"-l",""+l,"-n",""+n,"-k",""+k,"-M",""+M,"-S","--best","--sam-RG","ID:"+strand.name(),"--sam-RG","SM:"+strand.name(),"--nomaqround",ref.getAbsolutePath(),"-I", ""+I, "-X",""+X, "--fr", "--12", "-"};	
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
							logger.info("Start read feeding to alignment against "+ref+" log file:"+logFileName);
							try {
									while ((readsLine=readsStream.readLine())!=null && !shouldStop){
									
									ps.println(readsLine);
									
								}
								ps.flush();
								ps.close();
								
								logger.info("Finished read feeding to alignment against "+ref+" log file:"+logFileName);
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
				logger.info("Both alignments finished. Ambigous reads: "+tagCount);
				
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
					if (!sample.isPaired() && !tokensCT[0].equals(tokensGA[0])){
						throw new RuntimeException("BUG: reading two samrecords from CT and GA alignments with are a different read\nCT:"+CTLine+"\nGA:"+GALine);
					}else if (sample.isPaired() && !tokensCT[0].substring(0, tokensCT[0].length()-1).equals(tokensGA[0].substring(0, tokensGA[0].length()-1))){
						//throw new RuntimeException("BUG: reading two samrecords from CT and GA alignments with are a different read (ignoring last character)\nCT:"+CTLine+"\nGA:"+GALine);
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


