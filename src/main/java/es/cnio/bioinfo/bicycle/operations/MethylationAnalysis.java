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
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.management.ManagementFactory;
import java.lang.management.RuntimeMXBean;
import java.security.Permission;
import java.util.List;
import java.util.logging.Logger;

import es.cnio.bioinfo.bicycle.ErrorRateMode;
import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Reference;
import es.cnio.bioinfo.bicycle.Sample;
import es.cnio.bioinfo.bicycle.Tools;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment.Strand;

public class MethylationAnalysis {
	private static final Logger logger = Logger.getLogger(MethylationAnalysis.class.getName());
	
	private Project project;
	
	public MethylationAnalysis(Project p) {
		this.project = p;
	}
	

	
	public void analyzeWithErrorFromBarcodes(Reference reference, 
			Sample sample, 
			boolean trimreads,
			int trimuntil,
			boolean removeAmbiguous, 
			boolean removeBad,
			boolean removeClonal,
			boolean correctNonCG,
			int mindepth,
			double fdr,
			int nThreads,			
			List<File> bedFiles) throws IOException, InterruptedException{
		
		this.analyze(reference, sample, trimreads, trimuntil, removeAmbiguous, removeBad, removeClonal, correctNonCG, mindepth, fdr, nThreads, ErrorRateMode.from_barcodes, "", 0d, 0d, bedFiles);
		
	}
	
	public void analyzeWithErrorFromControlGenome(Reference reference, 
			Sample sample, 
			boolean trimreads,
			int trimuntil,
			boolean removeAmbiguous, 
			boolean removeBad,
			boolean removeClonal,
			boolean correctNonCG,
			int mindepth,
			double fdr,
			int nThreads,			
			List<File> bedFiles,
			String controlGenome) throws IOException, InterruptedException{
		
		this.analyze(reference, sample, trimreads, trimuntil, removeAmbiguous, removeBad, removeClonal, correctNonCG, mindepth, fdr, nThreads, ErrorRateMode.from_control_genome, controlGenome, 0d, 0d, bedFiles);
		
	}
	
	public void analyzeWithFixedErrorRate(Reference reference, 
			Sample sample, 
			boolean trimreads,
			int trimuntil,
			boolean removeAmbiguous, 
			boolean removeBad,
			boolean removeClonal,
			boolean correctNonCG,
			int mindepth,
			double fdr,
			int nThreads,			
			List<File> bedFiles,
			double watsonError, double crickError) throws IOException, InterruptedException{
		
		this.analyze(reference, sample, trimreads, trimuntil, removeAmbiguous, removeBad, removeClonal, correctNonCG, mindepth, fdr, nThreads, ErrorRateMode.FIXED, "", watsonError, crickError, bedFiles);
		
	}
	
	public File getMethylcytosinesFile(Reference reference, Sample sample){
		return new File(this.project.getOutputDirectory()+File.separator+sample.getName()+"_"+reference.getReferenceFile().getName()+".methylcytosines");
	}
	public File getMethylcytosinesVCFFile(Reference reference, Sample sample){
		return new File(this.project.getOutputDirectory()+File.separator+sample.getName()+"_"+reference.getReferenceFile().getName()+".methylcytosines.vcf");
	}
	private void analyze(Reference reference, 
			Sample sample, 
			boolean trimreads,
			int trimuntil,
			boolean removeAmbiguous, 
			boolean removeBad,
			boolean removeClonal,
			boolean correctNonCG,
			int mindepth,
			double fdr,
			int nThreads,
			ErrorRateMode errorMode,
			String controlGenome,
			double watsonError,
			double crickError,
			List<File> bedFiles) throws IOException, InterruptedException{
		
		BowtieAlignment ba = new BowtieAlignment(this.project);
		File samFileCT = ba.getAlignmentOutputFile(Strand.WATSON, sample, reference);
		File samFileGA = ba.getAlignmentOutputFile(Strand.CRICK, sample, reference);
		File fasta = reference.getReferenceFile();
	
		File sortedCT = new File(samFileCT.getAbsolutePath() + ".sorted.sam");
		
		sortSAM(samFileCT, sortedCT);
		File outputBamFileCT = buildBAMAndIndex(sortedCT, this.project.getSamtoolsDirectory());
		
		File sortedGA = new File(samFileGA.getAbsolutePath() + ".sorted.sam");
		sortSAM(samFileGA, sortedGA);
		File outputBamFileGA = buildBAMAndIndex(sortedGA, this.project.getSamtoolsDirectory());
		
		RuntimeMXBean runtimemxBean = ManagementFactory.getRuntimeMXBean();
	
		//String command = "java -Xmx1024M -cp "+runtimemxBean.getClassPath()+" org.broadinstitute.sting.gatk.CommandLineGATK ";
		String command = "-T ListerMethylation -I "+outputBamFileCT+" -I "+outputBamFileGA+" -R "+fasta+" -nt "+nThreads+" --outdir "+project.getOutputDirectory()+" --fdr "+fdr;
		
		command+=" --methylcytosinesfile "+getMethylcytosinesFile(reference, sample).getAbsolutePath();
		command+=" --methylcytosinesvcffile "+getMethylcytosinesVCFFile(reference, sample).getAbsolutePath();
		command+=" --summaryfile "+getSummaryFile(reference, sample).getAbsolutePath();
		command+=" --methylationwatsonfile "+getMethylationFile(Strand.WATSON, reference, sample).getAbsolutePath();
		command+=" --methylationcrickfile "+getMethylationFile(Strand.CRICK, reference, sample).getAbsolutePath();
		if (removeClonal){
			command+=" --removeclonal";
		}
		
		if (bedFiles!=null) for (File bedfile : bedFiles){
			command+=" -annotation:"+bedfile.getName()+",bed "+bedfile.getAbsolutePath();
		}
		
		if (errorMode == ErrorRateMode.from_barcodes){
			BarcodeErrorComputation bec = new BarcodeErrorComputation(sample);			
			double error = bec.computeErrorFromBarcodes();
			watsonError = error;
			crickError = error;
			command+=" --errorrate "+watsonError+","+crickError;
			
		}else if (errorMode == ErrorRateMode.from_control_genome){
			command+=" --controlgenome "+controlGenome;
		}else{
			command+=" --errorrate "+watsonError+","+crickError;
		}
		
		if (correctNonCG){
			command+=" --correctnoncg";
		}
		
		command+=" --mindepth "+mindepth;
		
		if (trimreads){
			command+=" --trim";
		}
		
		command+=" --read_filter Lister";
		if (removeAmbiguous){
			command+=" --removeambiguous";
		}
		if (trimreads){
			command+=" --trimuntil "+trimuntil;
		}
		if (removeBad){
			command+=" --removebad";
		}
		
		//logger.info("GATK command: "+command);
		//final Process p = Runtime.getRuntime().exec(command.split(" "));
		try{  
			forbidSystemExitCall() ;
//			System.setProperty("snappy.disable", "true");
			
			org.broadinstitute.sting.gatk.CommandLineGATK.main(command.split(" "));
		}catch( ExitTrappedException e ) {
			//SortSam seems to have a System.exit(0) at the end!
	    } finally {
	        enableSystemExitCall() ;
	    }
		/*Runtime.getRuntime().addShutdownHook(new Thread(){
			@Override
			public void run() {
				p.destroy();
			}
		});*/
		/*class StreamReader extends Thread{
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
		}*/
	}

	public File getMethylationFile(Strand strand, Reference reference,
			Sample sample) {
		return new File(this.project.getOutputDirectory()+File.separator+sample.getName()+"_"+reference.getReferenceFile().getName()+"_"+strand.name()+".methylation");
	}

	public File getMethylationVCFFile(Strand strand, Reference reference,
			Sample sample) {
		return new File(this.project.getOutputDirectory()+File.separator+sample.getName()+"_"+reference.getReferenceFile().getName()+"_"+strand.name()+".methylation.vcf");
	}

	public File getSummaryFile(Reference reference, Sample sample) {
		return new File(this.project.getOutputDirectory()+File.separator+sample.getName()+"_"+reference.getReferenceFile().getName()+".summary");
	}



	private static void sortSAM(File sam, File output) throws InterruptedException, IOException{
		//sort the sam
		
		if(new File(output.getAbsolutePath()).exists() && new File(output.getAbsolutePath()).lastModified() > sam.lastModified()){
			//System.out.println("skipping. found the corresponding sorted file, older than the unsorted input file");
			return;
		}
		System.out.println("Sorting " + sam.getAbsolutePath());
		String outfile = sam.getAbsolutePath() + ".sorted.sam";
		
		//String classLocation = MethylationAnalysis.class.getProtectionDomain().getCodeSource().getLocation().getFile();
		//String command = "java -Xmx1024M -Dsnappy.disable=true -jar "+classLocation+"../lib/SortSam.jar I="+ sam.getAbsolutePath()+ " O="+ outfile+" SO=coordinate TMP_DIR="+sam.getAbsoluteFile().getParentFile().getAbsolutePath();
		
		/*if (Runtime.getRuntime().exec((command).split(" ")).waitFor()!=0){
			//error
			if (new File(outfile).exists()){
				new File(outfile).delete();
			}
			throw new RuntimeException("Sort failed! command: "+command);
		}*/
		

		  
		try{  
			forbidSystemExitCall() ;
//			System.setProperty("snappy.disable", "true");
			net.sf.picard.sam.SortSam.main(new String[]{"I="+sam.getAbsolutePath(),"O="+outfile, "SO=coordinate", "TMP_DIR="+sam.getAbsoluteFile().getParentFile().getAbsolutePath()});
		}catch( ExitTrappedException e ) {
			//SortSam seems to have a System.exit(0) at the end!
	    } finally {
	        enableSystemExitCall() ;	       
	        
	    }
		
	}
	private static PrintStream err;
	private static void disableSystemErr(){
		if (err==null) err = System.err;
		System.setErr(new PrintStream(new OutputStream(){

			@Override
			public void write(int arg0) throws IOException {
				
			}
			
		}));
	}
	private static void enableSystemErr(){
		if (err!=null) System.setErr(err);
		err=null;
		
	}
	private static class ExitTrappedException extends SecurityException { }
	private static void forbidSystemExitCall() {
	    final SecurityManager securityManager = new SecurityManager() {
	    	private boolean hasExited=false;
	    	public void checkPermission(Permission perm, Object context) {
	    		//System.err.println(perm);
	    		
	    	};
	      public void checkPermission( Permission permission ) {
	    	  if( permission.getName().startsWith("exitVM") ) {
	    		  hasExited = true;
	    		  if (Integer.parseInt((permission.getName().split("\\.")[1]))!=0 && !hasExited){
	    			  
	    			  throw new RuntimeException("SortSam exited with status: "+Integer.parseInt((permission.getName().split("\\.")[1]))+". Please check if you have sufficient space in file system.");
	    		  }
	    		  disableSystemErr();
	    		  throw new ExitTrappedException();
	    	  }
	      }
	    } ;
	    
	    System.setSecurityManager( securityManager ) ;
	  }

	  private static void enableSystemExitCall() {
		
	    System.setSecurityManager( null ) ;
	    enableSystemErr();
	  }
	private static File buildBAMAndIndex(File samCT, File samtoolsDirectory) {
		File bam = new File(samCT.getAbsolutePath()+".bam");
		if (bam.exists() && bam.lastModified()>samCT.lastModified()){
		//	System.out.println("bam exists and is older. Skip");
			
		}else{
			logger.info("Building BAM for "+samCT);
			
			Tools.executeProcessWait(samtoolsDirectory.getAbsolutePath()+File.separator+"samtools view -S -b -o "+bam.getAbsolutePath()+" "+samCT.getAbsolutePath());
			logger.info("BAM built for "+samCT);
		}
		
		
		File bai = new File(bam.getAbsolutePath()+".bai");
		
		if (bai.exists() && bam.exists() && bai.lastModified()>bam.lastModified()){
		//	System.out.println("bai exists and is older. Skip");
			
		}else{
			logger.info("Building index for "+samCT);
			
			Tools.executeProcessWait(samtoolsDirectory.toString()+"/samtools index "+bam.getAbsolutePath());
			logger.info("index built for "+samCT);
		}
		
		
		return bam;
		
		
	}
	
		
}


