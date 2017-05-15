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

package es.cnio.bioinfo.bicycle.operations;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.security.Permission;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import es.cnio.bioinfo.bicycle.ErrorRateMode;
import es.cnio.bioinfo.bicycle.MethylationCall;
import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Reference;
import es.cnio.bioinfo.bicycle.RegionMethylation;
import es.cnio.bioinfo.bicycle.Sample;
import es.cnio.bioinfo.bicycle.StandardStreamsToLoggerRedirector;
import es.cnio.bioinfo.bicycle.Tools;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment.Strand;

public class MethylationAnalysis {
	private static final Logger logger = Logger.getLogger(MethylationAnalysis.class.getSimpleName());

	private Project project;

	public MethylationAnalysis(Project p) {
		this.project = p;
	}

	public Project getProject() {
		return project;
	}

	public File getMethylcytosinesFile(Reference reference, Sample sample) {
		return new File(this.project.getOutputDirectory() + File.separator + sample.getName() + "_" + reference
				.getReferenceFile().getName() + ".methylcytosines");
	}

	public File getMethylcytosinesVCFFile(Reference reference, Sample sample) {
		return new File(this.project.getOutputDirectory() + File.separator + sample.getName() + "_" + reference
				.getReferenceFile().getName() + ".methylcytosines.vcf");
	}

	public File getMethylationFile(Strand strand, Reference reference,
								   Sample sample) {
		return new File(this.project.getOutputDirectory() + File.separator + sample.getName() + "_" + reference
				.getReferenceFile().getName() + "_" + strand.name() + ".methylation");
	}

	public File getMethylationVCFFile(Strand strand, Reference reference,
									  Sample sample) {
		return new File(this.project.getOutputDirectory() + File.separator + sample.getName() + "_" + reference
				.getReferenceFile().getName() + "_" + strand.name() + ".methylation.vcf");
	}

	public File getSummaryFile(Reference reference, Sample sample) {
		return new File(this.project.getOutputDirectory() + File.separator + sample.getName() + "_" + reference
				.getReferenceFile().getName() + ".summary");
	}

	public File getMethylatedRegionsFile(Reference reference, Sample sample, File bed) {
		return new File(this.getMethylcytosinesFile(reference, sample).toString().replace("methylcytosines", "") + bed
				.getName().replace("bed", "METHYLATEDregions.txt"));
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
											 List<File> bedFiles) throws IOException, InterruptedException {

		this.analyze(reference, sample, trimreads, trimuntil, removeAmbiguous, removeBad, removeClonal, correctNonCG,
				mindepth, fdr, nThreads, ErrorRateMode.from_barcodes, "", 0d, 0d, bedFiles);
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
												  String controlGenome) throws IOException, InterruptedException {

		this.analyze(reference, sample, trimreads, trimuntil, removeAmbiguous, removeBad, removeClonal, correctNonCG,
				mindepth, fdr, nThreads, ErrorRateMode.from_control_genome, controlGenome, 0d, 0d, bedFiles);
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
										  double watsonError, double crickError) throws IOException,
			InterruptedException {

		this.analyze(reference, sample, trimreads, trimuntil, removeAmbiguous, removeBad, removeClonal, correctNonCG,
				mindepth, fdr, nThreads, ErrorRateMode.FIXED, "", watsonError, crickError, bedFiles);

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
						 List<File> bedFiles) throws IOException, InterruptedException {

		final String command = prepareGATKCommand(reference, sample, trimreads, trimuntil, removeAmbiguous, removeBad,
				removeClonal, correctNonCG, mindepth, fdr, nThreads, errorMode, controlGenome, watsonError,
				crickError, bedFiles);

		logger.info("Starting methylation analysis of sample " + sample.getName());
		//logger.info("GATK command: "+command);
		//final Process p = Runtime.getRuntime().exec(command.split(" "));

		StandardStreamsToLoggerRedirector redirector = null;


		try {
			forbidSystemExitCall();
			redirector = new StandardStreamsToLoggerRedirector(logger, Level.INFO,
					new StandardStreamsToLoggerRedirector.MessageFilter() {

						@Override
						public String filter(String msg) {
							//remove lines containing HelpFormatter
							if (msg.contains("HelpFormatter")) {
								return "";
							}
							if (msg.contains("RestStorageService")) {
								return "";
							}

							if (msg.contains("INFO")) {
								msg = msg.replaceAll("INFO.*- ", "");
							}
							return "GATK: " + msg;
						}
					});
			redirector.redirectStreams();

			org.broadinstitute.sting.gatk.CommandLineGATK.main(command.split(" "));

		} catch (ExitTrappedException e) {
			//SortSam seems to have a System.exit(0) at the end!
		} finally {
			redirector.restoreStreams();
			enableSystemExitCall();
		}

		// patch strange picard behaviour: remove empty directory (named as the user name) that is created in the
		// output directory if sort sams and sequence dictionary (ref_genomes/*.dict) are both created in this
		// process. This occurs the first time the methylationanalsysis is run, where sorted and .dict files have
		// to be created.
		for (File f : this.project.getOutputDirectory().listFiles()) {
			if (f.isDirectory() && f.listFiles().length == 0) {
				f.delete();
			}
		}

		writeRegionsMethylation(reference, sample, bedFiles);

		logger.info("Methylation analysis of sample " + sample.getName() + " OK");

	}

	private String prepareGATKCommand(Reference reference, Sample sample, boolean trimreads, int trimuntil,
									  boolean removeAmbiguous, boolean removeBad, boolean removeClonal, boolean
											  correctNonCG, int mindepth,
									  double fdr, int nThreads, ErrorRateMode errorMode, String controlGenome, double
											  watsonError,
									  double crickError, List<File> bedFiles) throws InterruptedException,
			IOException {

		BowtieAlignment ba = new BowtieAlignment(this.project);
		File samFileCT = ba.getAlignmentOutputFile(Strand.WATSON, sample, reference);
		File samFileGA = ba.getAlignmentOutputFile(Strand.CRICK, sample, reference);
		File fasta = reference.getReferenceFile();


		final PrintStream oldErr = System.err;
		// capture SortSam output and redirect it to our logger

		File sortedCT = new File(samFileCT.getAbsolutePath() + ".sorted.sam");
		File sortedGA = new File(samFileGA.getAbsolutePath() + ".sorted.sam");

		sortSAM(samFileCT, sortedCT);
		sortSAM(samFileGA, sortedGA);

		File outputBamFileCT = buildBAMAndIndex(sortedCT, this.project.getSamtoolsDirectory());
		File outputBamFileGA = buildBAMAndIndex(sortedGA, this.project.getSamtoolsDirectory());

		// RuntimeMXBean runtimemxBean = ManagementFactory.getRuntimeMXBean();
		// String command = "java -Xmx1024M -cp "+runtimemxBean.getClassPath()+"
		// org.broadinstitute.sting.gatk.CommandLineGATK ";

		String command = "-T ListerMethylation -I " + outputBamFileCT + " -I " + outputBamFileGA + " -R " + fasta
				+ " -nt " + nThreads + " --outdir " + project.getOutputDirectory() + " --fdr " + fdr;

		command += " --methylcytosinesfile " + getMethylcytosinesFile(reference, sample).getAbsolutePath();
		command += " --methylcytosinesvcffile " + getMethylcytosinesVCFFile(reference, sample).getAbsolutePath();
		command += " --summaryfile " + getSummaryFile(reference, sample).getAbsolutePath();
		command += " --methylationwatsonfile " + getMethylationFile(Strand.WATSON, reference, sample)
				.getAbsolutePath();
		command += " --methylationcrickfile " + getMethylationFile(Strand.CRICK, reference, sample).getAbsolutePath();
		if (removeClonal) {
			command += " --removeclonal";
		}

		if (bedFiles != null)
			for (File bedfile : bedFiles) {
				command += " -annotation:" + bedfile.getName() + ",bed " + bedfile.getAbsolutePath();
			}

		if (errorMode == ErrorRateMode.from_barcodes) {
			BarcodeErrorComputation bec = new BarcodeErrorComputation(sample);
			double error = bec.computeErrorFromBarcodes();
			watsonError = error;
			crickError = error;
			command += " --errorrate " + watsonError + "," + crickError;

		} else if (errorMode == ErrorRateMode.from_control_genome) {
			command += " --controlgenome " + controlGenome;
		} else {
			command += " --errorrate " + watsonError + "," + crickError;
		}

		if (correctNonCG) {
			command += " --correctnoncg";
		}

		command += " --mindepth " + mindepth;

		if (trimreads) {
			command += " --trim";
		}

		command += " --read_filter Lister";
		if (removeAmbiguous) {
			command += " --removeambiguous";
		}
		if (trimreads) {
			command += " --trimuntil " + trimuntil;
		}
		if (removeBad) {
			command += " --removebad";
		}

		return command;
	}

	private void writeRegionsMethylation(Reference reference, Sample sample, List<File> bedFiles) throws IOException {
		// <---- Cytosine METHYLATION PER ANNOTATED REGION ---->
		// IF BED FILES ARE AVAILABLE, then cytosine methylation per annotated
		// region is also calculated
		// added by osvaldo, 23Jan2016
		if (!bedFiles.isEmpty()) {
			// calculates methylation for each annotated region in each bed file
			// for the methylcytosines file


			// foreach bed file
			for (File bed : bedFiles) {
				logger.info("Checking methylation for regions annotated in " + bed.toString());
				BufferedWriter w = null;
				String annotationSet = bed.getName();
				try {
					// output file that will store methylation values for each
					// annotated region
					File outputFile = this.getMethylatedRegionsFile(reference, sample, bed);

					FileWriter o = new FileWriter(outputFile);

					w = new BufferedWriter(o);
					w.write(RegionMethylation.getMarshallHeader());
					w.newLine();
					w.flush();


					List<RegionMethylation> regionsMethylation =
							computeRegionsMethylation(reference, sample, annotationSet);

					for (RegionMethylation rM : regionsMethylation) {
						w.append(rM.marshall());
						w.newLine();
					}

				} catch (FileNotFoundException e) {
					e.printStackTrace();

				} catch (IOException e) {
					e.printStackTrace();

				} finally {
					if (w != null)
						w.close();
				}

			}
		} // if(!bedFiles.isEmpty())
	}

	public List<RegionMethylation> computeRegionsMethylation(Reference reference, Sample sample, String annotationSet)
			throws FileNotFoundException, IOException {
		File methylcytosinesFile = this.getMethylcytosinesFile(reference, sample);
		logger.info("Calculating methylation per annotated region for: " + methylcytosinesFile.toString().replaceAll
				(project.getOutputDirectory() + File.separator, Project.OUTPUT_DIRECTORY));

		List<RegionMethylation> regionsMethylation = new LinkedList<>();

		// each hash table contains separated values for Watson and
		// Crick
		Map<String, Map<String, Integer>> mCG = new LinkedHashMap<>();
		Map<String, Map<String, Integer>> mCHG = new LinkedHashMap<>();
		Map<String, Map<String, Integer>> mCHH = new LinkedHashMap<>();
		Map<String, Map<String, Integer>> depthCG = new LinkedHashMap<>();
		Map<String, Map<String, Integer>> depthCHG = new LinkedHashMap<>();
		Map<String, Map<String, Integer>> depthCHH = new LinkedHashMap<>();


		boolean firstLine = true;
		int columnOfInterest = -1;

		// opens the methylcytosine file
		FileReader f = new FileReader(methylcytosinesFile);
		BufferedReader b = new BufferedReader(f);
		String line;

		// for each line in the methylcytosines file (starting from
		// line 0, i.e., first line)
		while ((line = b.readLine()) != null) {

			String tokens[] = line.split("\t");


			// if positioned in first (header) line (line 0), it
			// finds out the column number that contains
			// the genomic annotations for the current bed file
			if (firstLine) {
				// checks all first line headers
				for (int pos = 0; pos < tokens.length; pos++) {
					if (tokens[pos].equals(annotationSet)) {
						// the current bed file is in column 'pos'
						columnOfInterest = pos;
						break;
					}
				}

				firstLine = false;

			} // if(firstLine)
			else {// not in first line (not header line)
				MethylationCall call = MethylationCall.unmarshall(line);
				if (!tokens[columnOfInterest].contains("N/A")) {
					List<String> regions = Arrays.asList(tokens[columnOfInterest].split("[|]"));


					for (String region : regions) {
						if (!mCG.containsKey(region)) {
							// initialization
							mCG.put(region, new HashMap<String, Integer>());
							mCG.get(region).put("WATSON", 0);
							mCG.get(region).put("CRICK", 0);

							mCHG.put(region, new HashMap<String, Integer>());
							mCHG.get(region).put("WATSON", 0);
							mCHG.get(region).put("CRICK", 0);

							mCHH.put(region, new HashMap<String, Integer>());
							mCHH.get(region).put("WATSON", 0);
							mCHH.get(region).put("CRICK", 0);

							depthCG.put(region, new HashMap<String, Integer>());
							depthCG.get(region).put("WATSON", 0);
							depthCG.get(region).put("CRICK", 0);

							depthCHG.put(region, new HashMap<String, Integer>());
							depthCHG.get(region).put("WATSON", 0);
							depthCHG.get(region).put("CRICK", 0);

							depthCHH.put(region, new HashMap<String, Integer>());
							depthCHH.get(region).put("WATSON", 0);
							depthCHH.get(region).put("CRICK", 0);
						}

						// value accumulation
						// accumulates methylated cytosines

						String strand = call.getStrand().name();//tokens[2];
						String methylationContext = call.getContext().name();//tokens[3];
						int depth = call.getDepth();//Integer.parseInt(tokens[4]);
						int methylation = call.getCytosines();// Integer.parseInt(tokens[6]);

						// accumulates new values for the current
						// annotated region
						switch (methylationContext) {
							case "CG":
								mCG.get(region).put(strand, mCG.get(region).get(strand) + methylation);
								depthCG.get(region).put(strand, depthCG.get(region).get(strand) + depth);
								break;

							case "CHG":
								mCHG.get(region).put(strand, mCHG.get(region).get(strand) + methylation);
								depthCHG.get(region).put(strand, depthCHG.get(region).get(strand) + depth);
								break;

							case "CHH":
								mCHH.get(region).put(strand, mCHH.get(region).get(strand) + methylation);
								depthCHH.get(region).put(strand, depthCHH.get(region).get(strand) + depth);
								break;

							default:
								System.out.print("[Invalid methylation context]: ");
								System.out.println(line);
						}

					}
				} // if(!tokens[columnOfInterest].equals("N/A")

			} // else{// not in first line (not header line)

		} // while((line=b.readLine())!=null)
		b.close();

		for (String region : mCG.keySet()) {
			regionsMethylation.add(
					new RegionMethylation(
							region,
							mCG.get(region).get("WATSON").intValue(),
							depthCG.get(region).get("WATSON").intValue(),
							mCHG.get(region).get("WATSON").intValue(),
							depthCHG.get(region).get("WATSON").intValue(),
							mCHH.get(region).get("WATSON").intValue(),
							depthCHH.get(region).get("WATSON").intValue(),
							mCG.get(region).get("CRICK").intValue(),
							depthCG.get(region).get("CRICK").intValue(),
							mCHG.get(region).get("CRICK").intValue(),
							depthCHG.get(region).get("CRICK").intValue(),
							mCHH.get(region).get("CRICK").intValue(),
							depthCHH.get(region).get("CRICK").intValue()

					));

		}

		return regionsMethylation;
	}


	private void sortSAM(File sam, File output) throws InterruptedException, IOException {
		//sort the sam

		if (new File(output.getAbsolutePath()).exists() && new File(output.getAbsolutePath()).lastModified() > sam
				.lastModified()) {
			//System.out.println("skipping. found the corresponding sorted file, older than the unsorted input file");
			return;
		}
		logger.info("Sorting " + sam.getAbsolutePath().replaceAll(project.getOutputDirectory() + File.separator,
				Project.OUTPUT_DIRECTORY));
		String outfile = sam.getAbsolutePath() + ".sorted.sam";

		StandardStreamsToLoggerRedirector redirector = null;
		try {
			forbidSystemExitCall();
			redirector = new StandardStreamsToLoggerRedirector(logger, Level.INFO,
					new StandardStreamsToLoggerRedirector.MessageFilter() {
						@Override
						public String filter(String msg) {
							if (!msg.contains("INFO")) return ""; //skip non INFO lines
							if (msg.contains("INPUT=")) return ""; //skip hello message
							return "SortSam (Picard): " + msg.replaceAll("INFO.*SortSam\\s+", "");
						}
					});
			redirector.redirectStreams();
			net.sf.picard.sam.SortSam.main(new String[]{"I=" + sam.getAbsolutePath(), "O=" + outfile, "SO=coordinate",
					"TMP_DIR=" + sam.getAbsoluteFile().getParentFile().getAbsolutePath()});
		} catch (ExitTrappedException e) {
			//SortSam seems to have a System.exit(0) at the end!
		} finally {
			enableSystemExitCall();
			redirector.restoreStreams();
		}

	}

	private static PrintStream err;

	private static void disableSystemErr(PrintStream originalStdErr) {
		if (err == null) err = originalStdErr;
		System.setErr(new PrintStream(new OutputStream() {

			@Override
			public void write(int arg0) throws IOException {

			}

		}));
	}

	private static void enableSystemErr() {
		if (err != null) System.setErr(err);
		err = null;

	}

	private static class ExitTrappedException extends SecurityException {
		private static final long serialVersionUID = 1L;
	}

	private static void forbidSystemExitCall() {
		final PrintStream systemErr = System.err;
		final SecurityManager securityManager = new SecurityManager() {
			private boolean hasExited = false;

			public void checkPermission(Permission perm, Object context) {
				//System.err.println(perm);

			}

			;

			public void checkPermission(Permission permission) {
				if (permission.getName().startsWith("exitVM")) {
					hasExited = true;
					if (Integer.parseInt((permission.getName().split("\\.")[1])) != 0 && !hasExited) {

						throw new RuntimeException("SortSam exited with status: " + Integer.parseInt((permission
								.getName().split("\\.")[1])) + ". Please check if you have sufficient space in file " +
								"system.");
					}
					disableSystemErr(systemErr);
					throw new ExitTrappedException();
				}
			}
		};

		System.setSecurityManager(securityManager);
	}

	private static void enableSystemExitCall() {

		System.setSecurityManager(null);
		enableSystemErr();
	}

	private File buildBAMAndIndex(File samCT, File samtoolsDirectory) {
		File bam = new File(samCT.getAbsolutePath() + ".bam");
		String samtoolsPath = (samtoolsDirectory != null ? samtoolsDirectory.getAbsolutePath() + File
				.separator : "") + "samtools";
		if (bam.exists() && bam.lastModified() > samCT.lastModified()) {
			//	System.out.println("bam exists and is older. Skip");

		} else {
			logger.info("Building BAM for " + samCT.toString().replaceAll(project.getOutputDirectory() + File
					.separator, ""));

			Tools.executeProcessWait(samtoolsPath + " view -S -b -o " + bam.getAbsolutePath() + " " + samCT
					.getAbsolutePath());
			logger.info("BAM built for " + samCT.toString().replaceAll(project.getOutputDirectory() + File.separator,
					""));
		}


		File bai = new File(bam.getAbsolutePath() + ".bai");

		if (bai.exists() && bam.exists() && bai.lastModified() > bam.lastModified()) {
			//	System.out.println("bai exists and is older. Skip");

		} else {
			logger.info("Building index for " + samCT.toString().replaceAll(project.getOutputDirectory() + File
					.separator, Project.OUTPUT_DIRECTORY));

			Tools.executeProcessWait(samtoolsPath + " index " + bam.getAbsolutePath());
			logger.info("Index built for " + samCT.toString().replaceAll(project.getOutputDirectory() + File
					.separator, Project.OUTPUT_DIRECTORY));
		}


		return bam;


	}


}


