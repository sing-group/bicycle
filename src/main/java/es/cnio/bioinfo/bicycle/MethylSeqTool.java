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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.RuntimeMXBean;
import java.util.Arrays;
import java.util.Date;
import java.util.List;

/**
A pipeline for methylation analysis
@author Copyright 2010<br>Daniel Gonz&aacute;lez-Pe&ntilde;a, Higher Technical School of Computer Engineering, University of Vigo, Ourense, Spain<br>Osvaldo Gra&ntilde;a, Bioinformatics Unit (CNIO)<br><br>Distributed under the terms of the GNU General Public License. (Free Software Foundation)
@version 1.0, may 2010<br>
@deprecated Use es.cnio.bioinfo.bicycle.Main class instead
*/ 
public class MethylSeqTool {
	private StringBuilder bowtieDirectory;
	private StringBuilder samtoolsDirectory;
	private StringBuilder workingDirectory,outputDirectory;
	private boolean a,b,c,d,e,f,h,j,l,m,n,t,u;
	private StringBuilder options,project,projectDirectory,projectFolder;
	private StringBuilder referenceDirectory,readsDirectory,reference,reads,bedFilesDirectory;
	private boolean create,run,delete;
	private String pBowtie,eBowtie,lBowtie,nBowtie,kBowtie,mBowtie,chunksBowtie,control_genome,solexaQ;
	private int trimMismatches;
	private StringBuilder bedFile;
	private ErrorRateMode errorMode;
	private double[] fixedErrorRates = {-1,-1,-1};
	
	private boolean calculateErrorFromBarcodes;
	
	
	
	public MethylSeqTool() {
		super();
		
		this.bowtieDirectory=null;	
		this.samtoolsDirectory=null;
		this.a=false;
		this.b=false;
		this.c=false;
		this.d=false;
		this.e=false;
		this.f=false;
		this.h=false;
		this.j=false;
		this.l=false;
		this.m=false;
		this.n=false;
		this.t=false;
		this.u=false;
		this.options=null;
		this.project=null;
		this.projectDirectory=null;
		this.create=false;
		this.run=false;
		this.delete=false;
		this.referenceDirectory=null;
		this.readsDirectory=null;
		this.bedFilesDirectory=null;
		this.reference=null;
		this.reads=null;
		this.projectFolder=null;
		this.workingDirectory=null;
		this.outputDirectory=null;
		this.pBowtie="1";//obligatorio correrlo así para correr 2 bowties monohilo en el paso P.h
		this.eBowtie="140";
		this.lBowtie="20";
		this.nBowtie="0";
		this.kBowtie="10";
		this.mBowtie="1"; // por defecto elimina las reads no únicas, es decir, que reportan más de m alineamientos		
		this.trimMismatches=4;
		this.chunksBowtie="--chunkmbs 64";
		this.control_genome=null;
		this.solexaQ="--solexa1.3-quals";
		this.calculateErrorFromBarcodes=false;
		this.errorMode = ErrorRateMode.from_control_genome;
		this.bedFile=null;
		
		
		
	}


	public static void main(String[] args) {
		// lo primero imprimo la linea ejecutada a un fichero de log que almaceno en el tmp
		RuntimeMXBean RuntimemxBean = ManagementFactory.getRuntimeMXBean();
		List<String> arguments = RuntimemxBean.getInputArguments();
		Date date=new Date();
		try{
			BufferedWriter executionLog=new BufferedWriter(new FileWriter("/tmp/execution.log",true));
			executionLog.write(date.toString());
			executionLog.newLine();
			executionLog.write("java -cp "+RuntimemxBean.getClassPath()+" ");
			for(int i=0;i<arguments.size();i++) executionLog.write(arguments.get(i)+" ");
			executionLog.write(" "+MethylSeqTool.class.getCanonicalName()+" ");
			for(int i=0;i<args.length;i++) executionLog.write(args[i]+" ");
			executionLog.newLine();
			executionLog.newLine();
			executionLog.flush();
			executionLog.close();
		}
		catch(Exception E){
				E.printStackTrace();
				System.exit(1);
		}


		// limpia pantalla
		System.out.println ((char)27 + "[2J");
		System.out.println("**************************************************************************************************************");
		System.out.println("* MethylSeqTool                                                                                              *");
		System.out.println("*                                                                                                            *");
		System.out.println("* Copyright 2011                                                                                             *");
		System.out.println("* Daniel González-Peña, Higher Technical School of Computer Engineering, University of Vigo                  *");
		System.out.println("* Osvaldo Graña, Bioinformatics Unit (CNIO)                                                                  *");
		System.out.println("* Distributed under the terms of the GNU General Public License. (Free Software Foundation)                  *");
		System.out.println("**************************************************************************************************************");
		MethylSeqTool P=new MethylSeqTool();
		MethylSeqOperations.FDR=0.01; // valor por defecto mientras el usuario no indique otro
		MethylSeqOperations.REMOVE_CLONAL=false;
		MethylSeqOperations.DEPTH_FILTER=1; // valor por defecto
		
		// me aseguro de no tener que disponer de un virtual frame buffer
		System.setProperty("java.awt.headless","true");
		
		for(int i=0;i<args.length;i++){
			if(args[i].equals("--create-project")) P.create=true;
			if(args[i].equals("--project-name") && args.length>i+1) P.project=new StringBuilder(args[i+1]);
			if(args[i].equals("--project-directory") && args.length>i+1) P.projectDirectory=new StringBuilder(args[i+1]);
			if(args[i].equals("--reference-directory") && args.length>i+1) P.referenceDirectory=new StringBuilder(args[i+1]);
			if(args[i].equals("--reads-directory") && args.length>i+1) P.readsDirectory=new StringBuilder(args[i+1]);
			if(args[i].equals("--bedFiles-directory") && args.length>i+1) P.bedFilesDirectory=new StringBuilder(args[i+1]);
			if(args[i].equals("--run")) P.run=true;
			if(args[i].equals("--reference") && args.length>i+1) P.reference=new StringBuilder(args[i+1]);
			if(args[i].equals("--reads") && args.length>i+1) P.reads=new StringBuilder(args[i+1]);
			if(args[i].equals("--bowtie-directory") && args.length>i+1) P.bowtieDirectory=new StringBuilder(args[i+1]);
			if(args[i].equals("--samtools-directory") && args.length>i+1) P.samtoolsDirectory=new StringBuilder(args[i+1]);
			if(args[i].equals("--bedFile") && args.length>i+1) P.bedFile=new StringBuilder(args[i+1]);
			if(args[i].equals("--options") && args.length>i+1) P.options=new StringBuilder(args[i+1]);
			if(args[i].equals("-p") && args.length>i+1) P.pBowtie=new StringBuilder("").append(args[i+1]).toString();
			if(args[i].equals("-e") && args.length>i+1) P.eBowtie=new StringBuilder("-e ").append(args[i+1]).toString();
			if(args[i].equals("-l") && args.length>i+1) P.lBowtie=new StringBuilder("-l ").append(args[i+1]).toString();
			if(args[i].equals("-n") && args.length>i+1) P.nBowtie=new StringBuilder("-n ").append(args[i+1]).toString();
			if(args[i].equals("-k") && args.length>i+1) P.kBowtie=new StringBuilder("-k ").append(args[i+1]).toString();
			if(args[i].equals("-m") && args.length>i+1) P.mBowtie=new StringBuilder("-m ").append(args[i+1]).toString();
			if(args[i].equals("--depth-filter") && args.length>i+1) MethylSeqOperations.DEPTH_FILTER=(new Integer(args[i+1])).intValue();
			if(args[i].equals("--FDR") && args.length>i+1){
					Double d=new Double(args[i+1]);
					MethylSeqOperations.FDR=d.doubleValue();					
			}
			if(args[i].equals("--errorRates") && args.length>i+1){
				if(args[i+1].indexOf(',')!=-1){
					P.errorMode = ErrorRateMode.FIXED;
					String[] numbers = args[i+1].split(",");
					if (numbers.length!=3){
						System.err.println("Error parsing: "+args[i+1]+". You need 3 error rates");
						System.exit(1);
					}
					try{
						P.fixedErrorRates[0] = Double.parseDouble(numbers[0]);
						P.fixedErrorRates[1] = Double.parseDouble(numbers[1]);
						P.fixedErrorRates[2] = Double.parseDouble(numbers[2]);
					}catch(NumberFormatException e){
						System.err.println("Unrecoginized numbers in "+args[i+1]);
						System.exit(1);
					}
				}else{
					try{
						P.errorMode =  ErrorRateMode.valueOf(args[i+1]);
						
					}catch(IllegalArgumentException e){
						System.err.println("Unrecoginized error rate calculation mode: "+args[i+1]);
						System.exit(1);
					}
				}
				
			}
			if(args[i].equals("--removeClonal")) MethylSeqOperations.REMOVE_CLONAL=true; 
			if(args[i].equals("--delete") && args.length>i+1) P.delete=true;
			if(args[i].equals("--control-genome") && args.length>i+1) P.control_genome=args[i+1];
			

		}
		
		// compruebo que la ruta incluye '/' al final
		if(P.projectDirectory!=null && P.projectDirectory.charAt(P.projectDirectory.length()-1)!='/') P.projectDirectory.append("/");
		if(P.referenceDirectory!=null && P.referenceDirectory.charAt(P.referenceDirectory.length()-1)!='/') P.referenceDirectory.append("/");
		if(P.readsDirectory!=null && P.readsDirectory.charAt(P.readsDirectory.length()-1)!='/') P.readsDirectory.append("/");
		if(P.bowtieDirectory!=null && P.bowtieDirectory.charAt(P.bowtieDirectory.length()-1)!='/') P.bowtieDirectory.append("/");
		if(P.samtoolsDirectory!=null && P.samtoolsDirectory.charAt(P.samtoolsDirectory.length()-1)!='/') P.samtoolsDirectory.append("/");
		if(P.bedFilesDirectory!=null && P.bedFilesDirectory.charAt(P.bedFilesDirectory.length()-1)!='/') P.bedFilesDirectory.append("/");
			

		// recojo las opciones seleccionadas
		if(P.options!=null && P.options.indexOf("a")!=-1) P.a=true;
		if(P.options!=null && P.options.indexOf("b")!=-1) P.b=true;
		if(P.options!=null && P.options.indexOf("c")!=-1) P.c=true;
		if(P.options!=null && P.options.indexOf("d")!=-1) P.d=true;
		if(P.options!=null && P.options.indexOf("e")!=-1) P.e=true;
		if(P.options!=null && P.options.indexOf("f")!=-1) P.f=true;
		if(P.options!=null && P.options.indexOf("h")!=-1) P.h=true;
		if(P.options!=null && P.options.indexOf("j")!=-1) P.j=true;
		if(P.options!=null && P.options.indexOf("l")!=-1) P.l=true;
		if(P.options!=null && P.options.indexOf("m")!=-1) P.m=true;
		if(P.options!=null && P.options.indexOf("n")!=-1) P.n=true;
		if(P.options!=null && P.options.indexOf("t")!=-1) P.t=true;
		if(P.options!=null && P.options.indexOf("u")!=-1){
			P.u=true;
			
		}
		
		if(P.options!=null && P.options.indexOf("s")!=-1) P.solexaQ="--solexa-quals";
		else if(P.options!=null && P.options.indexOf("p33")!=-1) P.solexaQ="--phred33-quals";
		
		// en caso de que el usuario quiere calcular el error a partir de los barcodes sin realizar de nuevo el paso 'c' (performs reads in-silico bisulfitation (replacing Cs with Ts)
		if(P.errorMode == ErrorRateMode.from_barcodes){
			P.calculateErrorFromBarcodes=true;			
		}		

		
		if(args.length==0 || (!P.create && !P.run && !P.delete) || (P.create && (P.projectDirectory==null || P.project==null || P.referenceDirectory==null || P.readsDirectory==null || P.samtoolsDirectory==null || P.bowtieDirectory==null)) || (P.run && (P.project==null || P.projectDirectory==null || P.reference==null || /*P.reads==null ||*/ P.errorMode==null || P.options==null)) || (P.delete && (P.projectDirectory==null || P.project==null))){
			
			System.out.println("\n\tArguments:");
			System.out.println("\t\t--create-project : creates a new project");
			System.out.println("\t\t--project-name : name of the project to create");
			System.out.println("\t\t--project-directory : path where to create the project");
			System.out.println("\t\t--reference-directory : directory that contains the reference files");
			System.out.println("\t\t--reads-directory : directory that contains the read files");
			System.out.println("\t\t--bedFiles-directory: directory that contains the annotation files");
			System.out.println("\t\t--reference : reference file **to use with --run (in case of several files, use commas to separate them: --reference chr1.fa,chr2.fa)");			
			//System.out.println("\t\t--reads : read file name **to use with --run (in case of several files, use commas to separate them: --reads s_1_sequence.txt,s_2_sequence.txt)");
			System.out.println("\t\t--control-genome : name of the genome used to estimate the error level for the binomial distribution (written as it appears in the fasta reference file)  **to use with --run");
			System.out.println("\t\t--FDR : FDR value (default value = 0.01)  **to use with --run");
			System.out.println("\t\t--removeClonal : removes clonal reads  **to use with --run");
			System.out.println("\t\t--errorRates : [<CG_rate>,<GHG_rate>,<GHH_rate> | from_control_genome | from_barcodes]\n\t\t\t\tSelects how to compute the error rate (1. Pre-fixed, 2. Cytosines at reference cytosines in a control genome, 3. Ratio of unconverted barcodes)");
			System.out.println("\t\t--bedFile : bed file with annotations **to use with --run (***they must be in the same directory as the read files.\n\t\t\t\tIn case of several files, use commas to separate them: --bedFile file1.bed,file2.bed)");
			System.out.println("\t\t--depth-filter: required read depth for a cytosine to be considered in the analysis (default 1)");
			System.out.println("\n\t\t--run : aligns reads and performs the analysis");
			System.out.println("\t\t--options : options to execute (for example: --options abcdefghij)");
			System.out.println("\t\t\t'a' - performs reference in-silico bisulfitation (CtoT and GtoA)");			
			System.out.println("\t\t\t'b' - performs reads in-silico bisulfitation (CtoT)");
			System.out.println("\t\t\t'c' - remove unconverted reads in step b (those were the bisulfite treatment failed, following the rule applied in Lister et al., Nature 2009)");
			System.out.println("\t\t\t'd' - remove reads with unconverted barcodes in step b");
			System.out.println("\t\t\t'e' - tells Bowtie to build indexes for both references, CtoT and GtoA");
			//System.out.println("\t\t\t'f' - aligns with Bowtie against both references (CtoT and GtoA)");			
			System.out.println("\t\t\t'h' - aligns with Bowtie against both references using multiple bowties (CtoT and GtoA)");
			System.out.println("\t\t\t'j' - analyses methylation levels over the Sam files");
			System.out.println("\t\t\t('u' - analyses methylation levels over the Sam files with the GATK-based walker)");
			System.out.println("\t\t\t'l' - trims reads to base preceding the 'x' mismatch");
			System.out.println("\t\t\t'm' - counts the number of aligned and non-aligned reads per bowtie output file and saves them in a csv file 'readCounts.csv'");
			System.out.println("\t\t\t'n' - eliminate reads that align to both CT and GA references");
			System.out.println("\t\t\t's1.3' - [--solexa1.3-quals] (default parameter for bowtie)");
			System.out.println("\t\t\t's' - [--solexa-quals]");
			System.out.println("\t\t\t'p33' - [--phred33-quals]\n");
			
			System.out.println("\t\t--delete : deletes an existing project");
			System.out.println("\t\t--samtools-directory : directory that contains the samtools programs (for example: /ALIGNERS/samtools-0.1.7a)");
			System.out.println("\t\t--bowtie-directory : directory that contains the bowtie aligner program (for example: /ALIGNERS/bowtie-0.12.5)");
			System.out.println("\t\t--bowtie default parameters : -t -p 1 --solexa1.3-quals -e 140 -l 20 -n 0 -k 10 --best --nomaqround -m 1");
			System.out.println("\t\t\tin order to change these values, type a new value for each parameter you want to change (for example: -l 25 -n 1)\n");
			
			System.out.println("\nExecution examples:");
			System.out.println("\n1.- Creating a project: java -cp ./lib/pileline.jar:./bin:./lib/apache-math/commons-math-2.1.jar:./lib/picard-tools-1.19/picard-1.19.jar:./lib/picard-tools-1.19/SortSam.jar -Dgenomeindex.headersize=10000000 es.cnio.bioinfo.methylationpipeline.MethylSeqTool --project-name myProjectName --project-directory /local/projects/ --samtools-directory /pathToSamTools/samtools-0.1.7a --bowtie-directory /pathToBowtie/bowtie-0.12.5 --reads-directory /pathToMyReadFiles/ --reference-directory /pathToMyReferenceFile/ --create-project\n");
			System.out.println("2.- Analyzing project data");
			System.out.println("2a) java -cp ./lib/pileline.jar:./bin:./lib/apache-math/commons-math-2.1.jar:./lib/picard-tools-1.19/picard-1.19.jar:./lib/picard-tools-1.19/SortSam.jar -Dgenomeindex.headersize=10000000  es.cnio.bioinfo.methylationpipeline.MethylSeqTool --project-name myProjectName --project-directory /local/projects/ --control-genome Ecoli --reference hg18.fa --errorRates from_barcodes --reads ES.txt,MEFs.txt,Lung.txt --bedFile UCSC_Longer_Isoform.bed,UCSC_CpG_islands.bed,super_RR.bed --run --options abckdefghj");
			System.out.println("2b) java -cp ./lib/pileline.jar:./bin:./lib/apache-math/commons-math-2.1.jar:./lib/picard-tools-1.19/picard-1.19.jar:./lib/picard-tools-1.19/SortSam.jar -Dgenomeindex.headersize=10000000  es.cnio.bioinfo.methylationpipeline.MethylSeqTool --project-name myProjectName --project-directory /local/projects/ --control-genome Ecoli --reference hg18_plus_Ecoli.fa --errorRates from_control_genome --reads ES.txt,MEFs.txt,Lung.txt --bedFile UCSC_Longer_Isoform.bed,UCSC_CpG_islands.bed,super_RR.bed --run --options abckdefghj");
			System.out.println("\n3.- Deleting a project:\njava -cp ./lib/pileline.jar:./bin:./lib/apache-math/commons-math-2.1.jar:./lib/picard-tools-1.19/picard-1.19.jar:./lib/picard-tools-1.19/SortSam.jar -Dgenomeindex.headersize=10000000  es.cnio.bioinfo.methylationpipeline.MethylSeqTool --project-name myProjectName --delete --projectDirectory /local/projects/ --project-name MyProject\n");
			
			
			
			
			
			System.out.println("\nCheck Parameter values:\n");
			System.out.println("--create-project: "+P.create);
			System.out.println("--project-name: "+P.project);
			System.out.println("--project-directory: "+P.projectDirectory);
			System.out.println("--reads-directory: "+P.readsDirectory);
			System.out.println("--reference-directory: "+P.referenceDirectory);
			System.out.println("--samtools-directory: "+P.samtoolsDirectory);
			System.out.println("--bowtie-directory: "+P.bowtieDirectory);
			
			
			
			System.exit(1);
			
			
			
			
			
		}

		//compruebo que el fichero (ó los ficheros) de reads existen en la ruta indicada
		// primero leo el fichero de configuracion
		P.projectFolder=new StringBuilder(P.projectDirectory).append(P.project).append("/");
		P.outputDirectory=new StringBuilder(P.projectFolder).append("output/");
		P.workingDirectory=new StringBuilder(P.projectFolder).append("workingDirectory/");
		
		

		if(P.create){// SI LO QUE QUIERE ES CREAR UN NUEVO PROYECTO			
		
			// creo el directorio de trabajo de este proyecto
			System.out.print(new StringBuilder("Creating the project folder: ").append(P.projectFolder).append("...... ").toString());
			// si ya esta creado: le aviso para que lo borre (menudo mamoncete)
			if(Tools.doesThisFileExist(P.projectFolder.toString())){
				System.err.println("\n[ERROR]: The project folder already exists, please remove it first\n");
				System.exit(1);
			}else{
				// lo creo
				Tools.mkdirs(P.projectFolder.toString());
				if(!Tools.doesThisFileExist(P.projectFolder.toString())){
					System.err.println(new StringBuilder("\n[ERROR]: unable to create ").append(P.projectFolder).toString());
					System.exit(1);
				}
				// creo tambien los subdirectorios				
				Tools.mkdirs(P.workingDirectory.toString());
				Tools.mkdirs(P.outputDirectory.toString());
				
				System.out.println("[OK]");
			}
			
			// creo el fichero de configuracion
			BufferedWriter wr=null;
			try {
				wr = new BufferedWriter(new FileWriter(new StringBuilder(P.projectFolder).append("config.txt").toString()));
				wr.write(new StringBuilder("project_name:").append(P.project).toString());
				wr.newLine();
				wr.write(new StringBuilder("project_directory:").append(P.projectDirectory).toString());
				wr.newLine();
			
				//compruebo que el directorio de referencias existe en la ruta indicada
				System.out.print(new StringBuilder("Looking for the reference(s) directory ").append(P.referenceDirectory).append("...... ").toString());
				if(Tools.doesThisFileExist(new StringBuilder(P.referenceDirectory).toString())){					
					System.out.println("[OK]");
					wr.write(new StringBuilder("reference_directory:").append(P.referenceDirectory).toString());
					wr.newLine();											
				}
				else{
						System.err.println(new StringBuilder("\n[ERROR]: unable to find ").append(P.referenceDirectory).append("\n").toString());
						System.exit(1);
				}
				
			
				//compruebo que el directorio de reads existe en la ruta indicada
				System.out.print(new StringBuilder("Looking for the read(s) directory ").append(P.readsDirectory).append("...... ").toString());
				if(Tools.doesThisFileExist(new StringBuilder(P.readsDirectory).toString())){
					System.out.println("[OK]");
					wr.write(new StringBuilder("reads_directory:").append(P.readsDirectory).toString());
					wr.newLine();
					
					
					P.reads=MethylSeqOperations.getListOfReadFiles(P.readsDirectory);
					String[] readFiles=MethylSeqOperations.searchForFiles(P.reads,P.readsDirectory,null);
					

				}
				else{
						System.err.println(new StringBuilder("\n[ERROR]: unable to find ").append(P.readsDirectory).append("\n").toString());
						System.exit(1);
				}
				
			
				wr.write(new StringBuilder("bowtie_directory:").append(P.bowtieDirectory).toString());
				wr.newLine();
				wr.write(new StringBuilder("samtools_directory:").append(P.samtoolsDirectory).toString());
				wr.newLine();
				//ahora escribo la linea ejecutada
				wr.write("java -cp "+RuntimemxBean.getClassPath()+" ");
				for(int i=0;i<arguments.size();i++) wr.write(arguments.get(i)+" ");
				wr.write(" "+MethylSeqTool.class.getCanonicalName()+" ");
				for(int i=0;i<args.length;i++) wr.write(args[i]+" ");
				wr.newLine();
			} catch (IOException e){
				System.err.println("[ERROR]: cannot create "+wr.toString()+"\n");				
				e.printStackTrace();
				System.exit(1);
			}finally{
				//cierro el BufferedWriter
				try {
					if (wr != null) {
						wr.flush();
						wr.close();
					}
				
				} catch (IOException ex) {
					System.err.println("[ERROR]: cannot close "+wr.toString()+"\n");
					ex.printStackTrace();
					System.exit(1);
				}
			}		
			
		}//if(P.create)
		
		else if(P.run){// si lo que quiere es bisulfitar, alinear, etc
			
			System.out.println("\n*** Using ["+P.solexaQ+"] with Bowtie");
			System.out.println("*** calculation of error rates: "+P.errorMode);
			System.out.println("*** Control genome: "+P.control_genome);
			System.out.println();
			
			if(P.errorMode==ErrorRateMode.from_control_genome && P.control_genome==null){
				System.out.println("[ERROR] You have selected to otbain the error rates from the control genome, but there is no control genome\n");
				System.exit(1);
			}
			
			// compruebo 1ero que existe el fichero config.txt para leer parametros, sino me paro y le aviso de ello
			File f=new File(P.projectFolder.toString(),"config.txt");
			if(f.exists()){
				String [] content=null;
				if(f.exists()){
					content=(Tools.readFile(new StringBuilder(P.projectFolder).append("config.txt").toString())).split("[\n]");
					for(int i=0;i<content.length;i++){						
						if(content[i].indexOf("project_directory")!=-1) P.projectDirectory=new StringBuilder((content[i].split("[:]"))[1]);
						if(content[i].indexOf("reference_directory")!=-1) P.referenceDirectory=new StringBuilder((content[i].split("[:]"))[1]);
						if(content[i].indexOf("reads_directory")!=-1) P.readsDirectory=new StringBuilder((content[i].split("[:]"))[1]);
						if(content[i].indexOf("bowtie_directory")!=-1) P.bowtieDirectory=new StringBuilder((content[i].split("[:]"))[1]);
						if(content[i].indexOf("samtools_directory")!=-1) P.samtoolsDirectory=new StringBuilder((content[i].split("[:]"))[1]);
					}
				}
				

				P.reads=MethylSeqOperations.getListOfReadFiles(P.readsDirectory);
				String[] readFiles=MethylSeqOperations.searchForFiles(P.reads,P.readsDirectory,null);
				
				// calcula el error a partir de los barcodes
				if(P.calculateErrorFromBarcodes){
					MethylSeqOperations.computeErrorRatesFromBarcodes(readFiles,P.readsDirectory,P.workingDirectory,P.outputDirectory);
				}

				//actualizo el contenido del fichero de configuración con la referencia/s utilizada/s, las reads analizadas y las opciones
				// de ejecucion pedidas
				BufferedWriter wr=null;
				try {
					wr = new BufferedWriter(new FileWriter(new StringBuilder(P.projectFolder).append("config.txt").toString(),true));					
					wr.write(new StringBuilder("reference:").append(P.reference).toString());
					wr.newLine();
					wr.write(new StringBuilder("reads:").append(P.reads).toString());
					wr.newLine();
					
					//ahora escribo la linea ejecutada
					wr.write("java -cp "+RuntimemxBean.getClassPath()+" ");
					for(int i=0;i<arguments.size();i++) wr.write(arguments.get(i)+" ");
					wr.write(" "+MethylSeqTool.class.getCanonicalName()+" ");
					for(int i=0;i<args.length;i++) wr.write(args[i]+" ");
					wr.newLine();
				} catch (IOException e){
					System.err.println("[ERROR]: cannot append information to "+wr.toString()+"\n");				
					e.printStackTrace();
					System.exit(1);
				}finally{
					//cierro el BufferedWriter
					try {
						if (wr != null) {
							wr.flush();
							wr.close();
						}
					
					} catch (IOException ex) {
						System.err.println("[ERROR]: cannot close "+wr.toString()+"\n");
						ex.printStackTrace();
						System.exit(1);
					}
				}
				
				
				
				
				// conversión de la referencia CtoT y GtoA
				String[] bisulfited_CT_references=null;
				String[] bisulfited_GA_references=null;
				String[] refs=null;
				if(P.a){
					//compruebo que el fichero (ó ficheros - separados por comas) de referencia existen en la ruta indicada
					refs=MethylSeqOperations.searchForFiles(P.reference,P.referenceDirectory,null);
					bisulfited_CT_references=MethylSeqOperations.performReferenceInsilicoBisulfitation_CtoT(refs,P.referenceDirectory,P.workingDirectory);
					bisulfited_GA_references=MethylSeqOperations.performReferenceInsilicoBisulfitation_GtoA(refs,P.referenceDirectory,P.workingDirectory);
					
					// actualizo cambios en el config si esta linea de informacion ya estaba presente, sino simplemente la añado
					//MethylSeqTools.updateConfigFile("bisulfited_CT_reference",P.projectFolder,bisulfited_CT_references);
				}
				
								
				// conversión de las reads CtoT
				String[] bisulfitedReads=null;
				if(P.b){
					//compruebo que el fichero (ó los ficheros) de reads existen en la ruta indicada
					bisulfitedReads=MethylSeqOperations.performReadsInsilicoBisulfitation_CtoT(readFiles,P.readsDirectory,P.workingDirectory,P.outputDirectory,P.c,P.d,false,null,null);
					
					// actualizo cambios en el config si esta linea de informacion ya estaba presente, sino simplemente la añado
					//MethylSeqTools.updateConfigFile("bisulfited_reads",P.projectFolder,bisulfitedReads);
				}				
				
				if(P.e){//INDEXACION de referencias CtoT y GtoA

					//si en esta ejecución no ha solicitado bisulfitar las referencias antes de indexar (a través de los pasos anteriores), tengo que asegurarme de que ya existen las referencias bisulfitadas
					if(bisulfited_CT_references==null) bisulfited_CT_references=MethylSeqOperations.getBisulfitedReferences(P.reference,P.workingDirectory,"bisulfited_CT_");
						
					// una vez comprobado que existen las referencias bisulfitadas, entonces a indexar....
					MethylSeqOperations.performBowtieReferenceIndexing(bisulfited_CT_references,P.workingDirectory,P.bowtieDirectory);
					
					//si en esta ejecución no ha solicitado bisulfitar las referencias antes de indexar (a través de los pasos anteriores), tengo que asegurarme de que ya existen las referencias bisulfitadas
					if(bisulfited_GA_references==null) bisulfited_GA_references=MethylSeqOperations.getBisulfitedReferences(P.reference,P.workingDirectory,"bisulfited_GA_");
					
					// si ha localizado todas las referencias bisulfitadas, entonces ha indexar
					MethylSeqOperations.performBowtieReferenceIndexing(bisulfited_GA_references,P.workingDirectory,P.bowtieDirectory);
				}
				
				//*****DEPRECATED METHOD****
				// alineamiento contra las reads contra ambas referencias (CtoT y GtoA)
				if(P.f){

					// antes de alinear debo comprobar que existen los ficheros indexados de referencias bisulfitadas
					String[] CT_indexedReferencesGenericName=MethylSeqOperations.getIndexedReferences(P.reference,P.workingDirectory,"bisulfited_CT_");
					String[] GA_indexedReferencesGenericName=MethylSeqOperations.getIndexedReferences(P.reference,P.workingDirectory,"bisulfited_GA_");
					

					
					// compruebo que existen las reads bisulfitadas					
					P.reads=MethylSeqOperations.getListOfReadFiles(P.readsDirectory);
					StringBuilder copy_P_reads=P.reads;// para quedarme una copia donde aun no he cambiado la barra de subdir '/' por '_'
					P.reads=MethylSeqOperations.replaceSlashByHyphen(P.reads);					
					
					if(bisulfitedReads==null) bisulfitedReads=MethylSeqOperations.searchForFiles(P.reads,P.workingDirectory,"bisulfited_CT_");
				
					// una vez que ya tengo los fastqs de reads del workingDirectory con su correspondiente prefijo 'bisulfited_CT_' o GA,
					// se me pueden presentar estas dos situaciones en el array bisulfitedReads:					
					// a) el directorio de READS NO contiene subdirectorios de muestras, todos los fastqs estan juntos, bisulfitedReads contiene
					// algo como 'ES_LIF.txt,ES_RA.txt' separados por comas.
					// b) el directorio de READS contiene subdirectorios de muestras, bisulfitedReads contiene algo como:
					// 	'control/ES_LIF.txt','experiment/ES_RA.txt','experiment/MEF.txt'
					// qué paso previo hay que hacer ahora para poder llamar al método performBowtieAlignmentwithSAMoutput?
					// en el caso a) hay que llamarlo una vez para cada fastq de entrada (antes se le pasaban todos juntos en un array), poniendo como
					// valor para alignmentOutputFile el nombre del fastq de entrada, que a su vez se pasa tambien como valor a String [] bisulfitedReads,
					// de uno en uno cada vez.
					// en el caso b) hay que llamar al metodo de alineamiento pasando como valor para 'String [] bisulfitedReads' todas las reads pertenecientes a la misma
					// muestra (mismo subdir), separadas por comas como siempre, y en el parametro alignmentOutputFile hay que pasar el nombre del fichero de salida, que
					// va a ser unico para todas las reads de la muestra (se alinean todas juntas y su salida va en un fichero comun). El nombre del fichero de salida
					// va a ser el nombre del subdiretorio donde estan las reads
					
					//extrae los nombres de los subirectorios (muestras) en caso de que haya. Si no hay simplemente contiene los fastqs del directorio READS,
					// vamos que es igual al propio bisulfitedReads.
					
					String [] samples=MethylSeqOperations.extractSubdirNames(copy_P_reads.toString().split(","));
					
					// ahora para llamar a bowtie:
					//a) en caso de que todas los fastqs estén en el directorio READs, le paso como argumento 'alignmentOutputFile' uno de esos ficheros fastqs, y como
					// bisulfited reads un array que en este caso contiene un unico elemento cuyo valor es ese mismo fichero.
					//b) en caso de que haya muestras (subdirs), el array 'samples' contiene las muestras, entonces en cada iteración al alinear cada muestra
					// le paso como argumento 'alignmentOutputFile' el nombre de la muestra en cuestion y como bisulfitedReads, el conjunto de ficheros de reads de 
					// esa muestra (subdir) separados por comas, de esta manera bowtie junta todos los ficheros de una misma muestra en uno solo y asi gana profundidad
					// de read
					
					for(int y=0;y<samples.length;y++){
						
						StringBuilder subconjuntoFastqsDeUnaMismaMuestra=new StringBuilder("");
						
						// si es un fastq perteneciente a ese sample en cuestion (deberia llevar el indicativo de sample incluido en el nombre, recordar- > experiment/MEFs_s_8_TtGATT-sequence.txt,experiment/Liver_s_8_TtAGTT-sequence.txt)
						for(int r=0;r<bisulfitedReads.length;r++) if(bisulfitedReads[r].indexOf(samples[y])!=-1) subconjuntoFastqsDeUnaMismaMuestra.append(bisulfitedReads[r]).append(",");
						
						String[] subsetBisulfitedReads=subconjuntoFastqsDeUnaMismaMuestra.toString().split(",");
						
						//Al loro, si hay muestras (subdirs) le adjunto al nombre el prefijo 'bisulfited_CT_', para que aparezca asi en el fichero de salida
						String alignmentOutputFileCT=samples[y];
						String alignmentOutputFileGA=samples[y];
						if(alignmentOutputFileCT.indexOf("bisulfited_CT_")==-1) alignmentOutputFileCT="bisulfited_CT_"+samples[y]+"_against_CTref";
						if(alignmentOutputFileGA.indexOf("bisulfited_CT_")==-1) alignmentOutputFileGA="bisulfited_CT_"+samples[y]+"_against_GAref";
						alignmentOutputFileCT=new StringBuilder(P.outputDirectory).append(alignmentOutputFileCT).append(".aligned").toString();
						alignmentOutputFileGA=new StringBuilder(P.outputDirectory).append(alignmentOutputFileGA).append(".aligned").toString();

						for (String ref : CT_indexedReferencesGenericName){
							MethylSeqOperations.performBowtieAlignmentwithSAMoutput(alignmentOutputFileCT,ref,subsetBisulfitedReads,P.workingDirectory,P.bowtieDirectory,P.outputDirectory,P.pBowtie,P.eBowtie,P.lBowtie,P.nBowtie,P.kBowtie,P.mBowtie,P.chunksBowtie,P.solexaQ);
						}
						for (String ref : GA_indexedReferencesGenericName){							
							MethylSeqOperations.performBowtieAlignmentwithSAMoutput(alignmentOutputFileGA,ref,subsetBisulfitedReads,P.workingDirectory,P.bowtieDirectory,P.outputDirectory,P.pBowtie,P.eBowtie,P.lBowtie,P.nBowtie,P.kBowtie,P.mBowtie,P.chunksBowtie,P.solexaQ);
						}
						MethylSeqOperations.replaceBowtieOutputWithOriginalReadsSAMoutput(P.outputDirectory,new String[]{alignmentOutputFileCT});
						MethylSeqOperations.replaceBowtieOutputWithOriginalReadsSAMoutput(P.outputDirectory,new String[]{alignmentOutputFileGA});
					}
					
					
				}

				// alineamiento contra las reads contra ambas referencias (CtoT y GtoA) usando la nueva version
				if(P.h){

					// antes de alinear debo comprobar que existen los ficheros indexados de referencias bisulfitadas
					String[] CT_indexedReferencesGenericName=MethylSeqOperations.getIndexedReferences(P.reference,P.workingDirectory,"bisulfited_CT_");
					String[] GA_indexedReferencesGenericName=MethylSeqOperations.getIndexedReferences(P.reference,P.workingDirectory,"bisulfited_GA_");
					

					
					// compruebo que existen las reads bisulfitadas					
					P.reads=MethylSeqOperations.getListOfReadFiles(P.readsDirectory);
					StringBuilder copy_P_reads=P.reads;// para quedarme una copia donde aun no he cambiado la barra de subdir '/' por '_'
					P.reads=MethylSeqOperations.replaceSlashByHyphen(P.reads);					
					
					if(bisulfitedReads==null) bisulfitedReads=MethylSeqOperations.searchForFiles(P.reads,P.workingDirectory,"bisulfited_CT_");
				
					// una vez que ya tengo los fastqs de reads del workingDirectory con su correspondiente prefijo 'bisulfited_CT_' o GA,
					// se me pueden presentar estas dos situaciones en el array bisulfitedReads:					
					// a) el directorio de READS NO contiene subdirectorios de muestras, todos los fastqs estan juntos, bisulfitedReads contiene
					// algo como 'ES_LIF.txt,ES_RA.txt' separados por comas.
					// b) el directorio de READS contiene subdirectorios de muestras, bisulfitedReads contiene algo como:
					// 	'control/ES_LIF.txt','experiment/ES_RA.txt','experiment/MEF.txt'
					// qué paso previo hay que hacer ahora para poder llamar al método performBowtieAlignmentwithSAMoutput?
					// en el caso a) hay que llamarlo una vez para cada fastq de entrada (antes se le pasaban todos juntos en un array), poniendo como
					// valor para alignmentOutputFile el nombre del fastq de entrada, que a su vez se pasa tambien como valor a String [] bisulfitedReads,
					// de uno en uno cada vez.
					// en el caso b) hay que llamar al metodo de alineamiento pasando como valor para 'String [] bisulfitedReads' todas las reads pertenecientes a la misma
					// muestra (mismo subdir), separadas por comas como siempre, y en el parametro alignmentOutputFile hay que pasar el nombre del fichero de salida, que
					// va a ser unico para todas las reads de la muestra (se alinean todas juntas y su salida va en un fichero comun). El nombre del fichero de salida
					// va a ser el nombre del subdiretorio donde estan las reads
					
					//extrae los nombres de los subirectorios (muestras) en caso de que haya. Si no hay simplemente contiene los fastqs del directorio READS,
					// vamos que es igual al propio bisulfitedReads.
					
					String [] samples=MethylSeqOperations.extractSubdirNames(copy_P_reads.toString().split(","));
					
					// ahora para llamar a bowtie:
					//a) en caso de que todas los fastqs estén en el directorio READs, le paso como argumento 'alignmentOutputFile' uno de esos ficheros fastqs, y como
					// bisulfited reads un array que en este caso contiene un unico elemento cuyo valor es ese mismo fichero.
					//b) en caso de que haya muestras (subdirs), el array 'samples' contiene las muestras, entonces en cada iteración al alinear cada muestra
					// le paso como argumento 'alignmentOutputFile' el nombre de la muestra en cuestion y como bisulfitedReads, el conjunto de ficheros de reads de 
					// esa muestra (subdir) separados por comas, de esta manera bowtie junta todos los ficheros de una misma muestra en uno solo y asi gana profundidad
					// de read
					
					System.out.println("samples: "+Arrays.toString(samples));
					for(int y=0;y<samples.length;y++){
						
						StringBuilder subconjuntoFastqsDeUnaMismaMuestra=new StringBuilder("");
						
						// si es un fastq perteneciente a ese sample en cuestion (deberia llevar el indicativo de sample incluido en el nombre, recordar- > experiment/MEFs_s_8_TtGATT-sequence.txt,experiment/Liver_s_8_TtAGTT-sequence.txt)
						for(int r=0;r<bisulfitedReads.length;r++) if(bisulfitedReads[r].indexOf(samples[y])!=-1) subconjuntoFastqsDeUnaMismaMuestra.append(bisulfitedReads[r]).append(",");
						
						String[] subsetBisulfitedReads=subconjuntoFastqsDeUnaMismaMuestra.toString().split(",");
						
						//Al loro, si hay muestras (subdirs) le adjunto al nombre el prefijo 'bisulfited_CT_', para que aparezca asi en el fichero de salida
						String alignmentOutputFileCT=samples[y];
						String alignmentOutputFileGA=samples[y];
						if(alignmentOutputFileCT.indexOf("bisulfited_CT_")==-1) alignmentOutputFileCT="bisulfited_CT_"+samples[y]+"_against_CTref";
						if(alignmentOutputFileGA.indexOf("bisulfited_CT_")==-1) alignmentOutputFileGA="bisulfited_CT_"+samples[y]+"_against_GAref";
						alignmentOutputFileCT=new StringBuilder(P.outputDirectory).append(alignmentOutputFileCT).append(".aligned").toString();
						alignmentOutputFileGA=new StringBuilder(P.outputDirectory).append(alignmentOutputFileGA).append(".aligned").toString();

						for (int i = 0; i<CT_indexedReferencesGenericName.length; i++){
							String refCT = CT_indexedReferencesGenericName[i];
							String refGA = GA_indexedReferencesGenericName[i];
							try {
								MethylSeqOperations.performBowtieAlignmentwithAmbigousTaggedSAMoutput(alignmentOutputFileCT, alignmentOutputFileGA, refCT, refGA, subsetBisulfitedReads, P.workingDirectory, P.bowtieDirectory, P.outputDirectory, P.pBowtie, P.eBowtie, P.lBowtie, P.nBowtie, P.kBowtie, P.mBowtie, P.chunksBowtie, P.solexaQ);
							} catch (IOException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
							
						}						
						//MethylSeqOperations.replaceBowtieOutputWithOriginalReadsSAMoutput(P.outputDirectory,new String[]{alignmentOutputFileCT});
						//MethylSeqOperations.replaceBowtieOutputWithOriginalReadsSAMoutput(P.outputDirectory,new String[]{alignmentOutputFileGA});
					}
					
					
				}
				// compruebo si he recibido ficheros bed de entrada con anotaciones. Si es asi me los apunto
				String[] bedFiles=null;				
				if(P.bedFile!=null) bedFiles=MethylSeqOperations.searchForFiles(P.bedFile,P.bedFilesDirectory,null);
				
				if(P.j){
					System.out.println("***Required read depth: "+MethylSeqOperations.DEPTH_FILTER);
					// antes de alinear debo comprobar que existen los ficheros indexados de referencias bisulfitadas
					String[] CT_indexedReferencesGenericName=MethylSeqOperations.getIndexedReferences(P.reference,P.workingDirectory,"bisulfited_CT_");
					String[] GA_indexedReferencesGenericName=MethylSeqOperations.getIndexedReferences(P.reference,P.workingDirectory,"bisulfited_GA_");
										
					/* compruebo que no hayan desaparecido los fastqs de reads (que no los haya borrado el usuario en una posterior ejecución de
					este paso)*/
					P.reads=MethylSeqOperations.getListOfReadFiles(P.readsDirectory);
					//reemplazo para el caso de que sean fastqs que están directamente en el directorio de reads,
					// si no hay '/' no interfiere
					String [] samples=MethylSeqOperations.extractSubdirNames(P.reads.toString().split(","));
					
					
					// recupera el nombre de los ficheros con el valor de error rate sobre los barcodes
					//NOTA: en el caso de que haya varios fastqs para una misma muestra, es justo en este paso donde va a calcular
					// el error conjunto de todos los fastqs de esa muestra (como si fuera un único fastq)
					String [] errorRateFiles=MethylSeqOperations.getErrorRateFileNames(P.reads.toString().split(","),P.readsDirectory,P.outputDirectory,P.control_genome);
					
					if(P.control_genome==null){ // en caso de haber genoma de control no realiza este paso						
						System.out.println("**** Considered Error Rates files:");
						for(int q=0;q<errorRateFiles.length;q++)System.out.println((q+1)+": "+errorRateFiles[q]);
					}
					P.reads=MethylSeqOperations.replaceSlashByHyphen(P.reads);
					
					// ahora para llamar a bowtie:
					//a) en caso de que todas los fastqs estén en el directorio READs, le paso como argumento 'alignmentOutputFile' uno de esos ficheros fastqs, y como
					// bisulfited reads un array que en este caso contiene un unico elemento cuyo valor es ese mismo fichero.
					//b) en caso de que haya muestras (subdirs), el array 'samples' contiene las muestras, entonces en cada iteración al alinear cada muestra
					// le paso como argumento 'alignmentOutputFile' el nombre de la muestra en cuestion y como bisulfitedReads, el conjunto de ficheros de reads de 
					// esa muestra (subdir) separados por comas, de esta manera bowtie junta todos los ficheros de una misma muestra en uno solo y asi gana profundidad
					// de read
					
					for(int y=0;y<samples.length;y++){
						
						//StringBuilder subconjuntoFastqsDeUnaMismaMuestra=new StringBuilder("");
						
						// si es un fastq perteneciente a ese sample en cuestion (deberia llevar el indicativo de sample incluido en el nombre, recordar- > experiment/MEFs_s_8_TtGATT-sequence.txt,experiment/Liver_s_8_TtAGTT-sequence.txt)
						//for(int r=0;r<bisulfitedReads.length;r++) if(bisulfitedReads[r].indexOf(samples[y])!=-1) subconjuntoFastqsDeUnaMismaMuestra.append(bisulfitedReads[r]).append(",");
						
						//String[] subsetBisulfitedReads=subconjuntoFastqsDeUnaMismaMuestra.toString().split(",");
						
						//Al loro, si hay muestras (subdirs) le adjunto al nombre el prefijo 'bisulfited_CT_', para que aparezca asi en el fichero de salida
						String alignmentOutputFileCT=samples[y];
						String alignmentOutputFileGA=samples[y];
						if(alignmentOutputFileCT.indexOf("bisulfited_CT_")==-1) alignmentOutputFileCT="bisulfited_CT_"+samples[y]+"_against_CTref";
						if(alignmentOutputFileGA.indexOf("bisulfited_CT_")==-1) alignmentOutputFileGA="bisulfited_CT_"+samples[y]+"_against_GAref";
						alignmentOutputFileCT=new StringBuilder(P.outputDirectory).append(alignmentOutputFileCT).append(".aligned").toString();
						alignmentOutputFileGA=new StringBuilder(P.outputDirectory).append(alignmentOutputFileGA).append(".aligned").toString();
						
						MethylSeqOperations.analyzeMethylation(P.outputDirectory, P.workingDirectory, P.referenceDirectory,P.readsDirectory,MethylSeqOperations.searchForFiles(P.reference, P.referenceDirectory, null),new String[]{alignmentOutputFileCT,alignmentOutputFileGA}, errorRateFiles[y],P.control_genome, P.samtoolsDirectory, P.errorMode, P.fixedErrorRates,P.bedFilesDirectory,bedFiles);
						
					}

				}
				
				if(P.u){
					// antes de alinear debo comprobar que existen los ficheros indexados de referencias bisulfitadas
					String[] CT_indexedReferencesGenericName=MethylSeqOperations.getIndexedReferences(P.reference,P.workingDirectory,"bisulfited_CT_");
					String[] GA_indexedReferencesGenericName=MethylSeqOperations.getIndexedReferences(P.reference,P.workingDirectory,"bisulfited_GA_");
										
					/* compruebo que no hayan desaparecido los fastqs de reads (que no los haya borrado el usuario en una posterior ejecución de
					este paso)*/
					P.reads=MethylSeqOperations.getListOfReadFiles(P.readsDirectory);
					//reemplazo para el caso de que sean fastqs que están directamente en el directorio de reads,
					// si no hay '/' no interfiere
					String [] samples=MethylSeqOperations.extractSubdirNames(P.reads.toString().split(","));
					
					
					// recupera el nombre de los ficheros con el valor de error rate
					//NOTA: en el caso de que haya varios fastqs para una misma muestra, es justo en este paso donde va a calcular
					// el error conjunto de todos los fastqs de esa muestra (como si fuera un único fastq)
					String [] errorRateFiles=MethylSeqOperations.getErrorRateFileNames(P.reads.toString().split(","),P.readsDirectory,P.outputDirectory,P.control_genome);
					P.reads=MethylSeqOperations.replaceSlashByHyphen(P.reads);
					
					// ahora para llamar a bowtie:
					//a) en caso de que todas los fastqs estén en el directorio READs, le paso como argumento 'alignmentOutputFile' uno de esos ficheros fastqs, y como
					// bisulfited reads un array que en este caso contiene un unico elemento cuyo valor es ese mismo fichero.
					//b) en caso de que haya muestras (subdirs), el array 'samples' contiene las muestras, entonces en cada iteración al alinear cada muestra
					// le paso como argumento 'alignmentOutputFile' el nombre de la muestra en cuestion y como bisulfitedReads, el conjunto de ficheros de reads de 
					// esa muestra (subdir) separados por comas, de esta manera bowtie junta todos los ficheros de una misma muestra en uno solo y asi gana profundidad
					// de read
					
					for(int y=0;y<samples.length;y++){
						
						//StringBuilder subconjuntoFastqsDeUnaMismaMuestra=new StringBuilder("");
						
						// si es un fastq perteneciente a ese sample en cuestion (deberia llevar el indicativo de sample incluido en el nombre, recordar- > experiment/MEFs_s_8_TtGATT-sequence.txt,experiment/Liver_s_8_TtAGTT-sequence.txt)
						//for(int r=0;r<bisulfitedReads.length;r++) if(bisulfitedReads[r].indexOf(samples[y])!=-1) subconjuntoFastqsDeUnaMismaMuestra.append(bisulfitedReads[r]).append(",");
						
						//String[] subsetBisulfitedReads=subconjuntoFastqsDeUnaMismaMuestra.toString().split(",");
						
						//Al loro, si hay muestras (subdirs) le adjunto al nombre el prefijo 'bisulfited_CT_', para que aparezca asi en el fichero de salida
						String alignmentOutputFileCT=samples[y];
						String alignmentOutputFileGA=samples[y];
						if(alignmentOutputFileCT.indexOf("bisulfited_CT_")==-1) alignmentOutputFileCT="bisulfited_CT_"+samples[y]+"_against_CTref";
						if(alignmentOutputFileGA.indexOf("bisulfited_CT_")==-1) alignmentOutputFileGA="bisulfited_CT_"+samples[y]+"_against_GAref";
						alignmentOutputFileCT=new StringBuilder(P.outputDirectory).append(alignmentOutputFileCT).append(".aligned").toString();
						alignmentOutputFileGA=new StringBuilder(P.outputDirectory).append(alignmentOutputFileGA).append(".aligned").toString();
						MethylSeqOperations.analyzeMethylationGATK(P.outputDirectory, P.workingDirectory, P.referenceDirectory,P.readsDirectory,MethylSeqOperations.searchForFiles(P.reference, P.referenceDirectory, null), new String[]{alignmentOutputFileCT,alignmentOutputFileGA}, errorRateFiles[y],P.control_genome, P.samtoolsDirectory, P.errorMode, P.fixedErrorRates, Integer.parseInt(P.pBowtie), P.l, P.n, P.c,P.bedFilesDirectory,bedFiles);
					}
					
					
					
				}
				
				
				// hago las cuentas
				if(P.m){
					MethylSeqOperations.count(P.outputDirectory);					
				}
				
				
			}//if(Tools.doesThisFileExist(new StringBuilder(P.projectFolder).append("config.txt").toString()))
			else{
				System.err.println(new StringBuilder("[ERROR]: the file ").append(P.projectFolder).append("config.txt does not exist\n").toString());
				System.exit(1);		
			}
			
		}//if(P.run)
		
		else if(P.delete){// me cargo el proyecto
			System.out.print(new StringBuilder("Deleting ").append(P.projectFolder).append("...... ").toString());
			String [] command=new String[]{"rm","-r",P.projectFolder.toString()};
			int result=Tools.executeProcessWait(command,null,null);
			
			if(result==0) System.out.println("[OK]");
			else{
				System.err.print("\n[ERROR]: Unable to perform:");
				for(int h=0;h<command.length;h++) System.err.print(new StringBuilder(command[h]).append(" ").toString());
				System.err.println("\n");
			}
			
		}
		
		System.out.println("\n[FINISHED]\n");
	}
	
	
	

}


