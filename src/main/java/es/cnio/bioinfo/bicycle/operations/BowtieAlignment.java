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

import static java.lang.Math.max;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
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
import es.cnio.bioinfo.bicycle.StandardStreamsToLoggerRedirector;
import es.cnio.bioinfo.bicycle.Tools;

public class BowtieAlignment {

	public enum PairedEndOrientation {
		FR("--fr"), FF("--ff"), RF("--rf");
		private String bowtieParameter;

		private PairedEndOrientation(String bowtieParameter) {
			this.bowtieParameter = bowtieParameter;
		}

		public String getBowtieParameter() {
			return bowtieParameter;
		}
	}

	public enum Strand {
		WATSON, CRICK,
	}

	private static final Logger logger = Logger.getLogger(BowtieAlignment.class.getSimpleName());


	private Project project;

	public enum Bowtie1Quals {
		AFTER_1_3("--solexa1.3-quals", "solexa1.3"),
		BEFORE_1_3("--solexa-quals", "solexa"),
		PHRED_33("--phred33-quals", "phred33"),
		PHRED_64("--phred64-quals", "phred64"),
		INTEGER("--integer-quals", "integer");

		private final String bicycleParameterValue;
		private final String parameterValue;

		Bowtie1Quals(String parameterValue, String bicycleParameterValue) {

			this.parameterValue = parameterValue;
			this.bicycleParameterValue = bicycleParameterValue;
		}

		public String getParameterValue() {
			return parameterValue;
		}

		public static Bowtie1Quals parseQuals(String qualsString) {
			for (Bowtie1Quals qualityValue : Bowtie1Quals.values()) {
				if (qualityValue.bicycleParameterValue.equalsIgnoreCase(qualsString)) {
					return qualityValue;
				}
			}
			throw new IllegalArgumentException("invalid qualities parameter: " + qualsString);
		}
	}

	public enum Bowtie2Quals {
		BEFORE_1_3("--solexa-quals", "solexa"),
		PHRED64("--phred64", "phred64"),
		PHRED_33("--phred33", "phred33"),
		INT_QUALS("--int-quals", "int");

		private final String bicycleParameterValue;
		private final String parameterValue;

		Bowtie2Quals(String parameterValue, String bicycleParameterValue) {

			this.parameterValue = parameterValue;
			this.bicycleParameterValue = bicycleParameterValue;
		}

		public String getParameterValue() {
			return parameterValue;
		}

		public static Bowtie2Quals parseQuals(String qualsString) {
			for (Bowtie2Quals qualityValue : Bowtie2Quals.values()) {
				if (qualityValue.bicycleParameterValue.equalsIgnoreCase(qualsString)) {
					return qualityValue;
				}
			}
			throw new IllegalArgumentException("invalid qualities parameter: " + qualsString);
		}
	}

	interface AlignmentScoreFunction {
		int getScore(String SAMLine);
	}

	private static class Bowtie2ScoreFunction implements AlignmentScoreFunction {
		// AS: Alignment score (bowtie 2)
		Pattern scorePattern = Pattern.compile("\\tAS:i:([^\\n\\t]+)");

		@Override
		public int getScore(String SAMLine) {
			Matcher matcher = scorePattern.matcher(SAMLine);
			matcher.find();
			return Integer.parseInt(matcher.group(1));
		}
	}

	private static class Bowtie1ScoreFunction implements AlignmentScoreFunction {

		// NM: edit distance (bowtie 1). The greater the worse, so the score is inverted
		Pattern scorePattern = Pattern.compile("\\tNM:i:([^\\n\\t]+)");

		@Override
		public int getScore(String SAMLine) {
			Matcher matcher = scorePattern.matcher(SAMLine);
			matcher.find();
			return -1 * Integer.parseInt(matcher.group(1));
		}
	}

	public BowtieAlignment(Project p) {
		this.project = p;
	}


	public void buildBowtieIndex(Reference reference) throws IOException {
		buildBowtieIndex(reference, 1);
	}

	public void buildBowtie2Index(Reference reference) {
		buildBowtieIndex(reference, 2);
	}

	private void buildBowtieIndex(Reference reference, int bowtieVersion) {
		ReferenceBisulfitation rb = new ReferenceBisulfitation(this.project);


		File bisulfitedReferenceCT = rb.getBisulfitedReference(ReferenceBisulfitation.Replacement.CT, reference);
		File bisulfitedReferenceGA = rb.getBisulfitedReference(ReferenceBisulfitation.Replacement.GA, reference);

		if (!bisulfitedReferenceCT.exists()) {
			throw new IllegalArgumentException("Reference CtoT bisulfitation " + bisulfitedReferenceCT + " could not" +
					" " +
					"be found. Please, perform reference in-silico bisulfitation first");
		}

		if (!bisulfitedReferenceGA.exists()) {
			throw new IllegalArgumentException("Reference GtoA bisulfitation " + bisulfitedReferenceGA + " could not" +
					" " +
					"be found. Please, perform reference in-silico bisulfitation first");
		}
		buildBowtieIndex(bisulfitedReferenceCT, bowtieVersion);

		buildBowtieIndex(bisulfitedReferenceGA, bowtieVersion);
	}

	private void buildBowtieIndex(File bisulfitedReference, int bowtieVersion) {
		logger.info("Building Bowtie index for " + bisulfitedReference.toString().replaceAll(project
				.getWorkingDirectory() + File.separator, Project.WORKING_DIRECTORY));

		String bowtieBuildPath = "";

		if (bowtieVersion == 1) {
			bowtieBuildPath = (this.project.getBowtieDirectory() != null ? this.project.getBowtieDirectory()
					.getAbsolutePath() + File
					.separator : "") + "bowtie-build";
		} else {
			bowtieBuildPath = (this.project.getBowtie2Directory() != null ? this.project.getBowtie2Directory()
					.getAbsolutePath() + File
					.separator : "") + "bowtie2-build";
		}


		String[] command = new String[]{bowtieBuildPath, bisulfitedReference.getAbsolutePath(), bisulfitedReference
				.getAbsolutePath()};

		StandardStreamsToLoggerRedirector redirector = new StandardStreamsToLoggerRedirector(logger, Level.INFO,
				new StandardStreamsToLoggerRedirector.MessageFilter() {

					@Override
					public String filter(String msg) {
						return "[bowtie-build]: " + msg;
					}
				},
				logger, Level.INFO,
				new StandardStreamsToLoggerRedirector.MessageFilter() {

					@Override
					public String filter(String msg) {
						return "[bowtie-build]: " + msg;
					}
				}
		);
		redirector.redirectStreams();

		try {
			int result = Tools.executeProcessWait(command, System.out, System.err);
			if (result == 0) logger.info("Bowtie index build OK");
			else {
				String commandString = "";

				for (int j = 0; j < command.length; j++)
					commandString += command[j] + " ";

				throw new RuntimeException("Error during bowtie index build. Command was: " + commandString);
			}
		} finally {
			redirector.restoreStreams();
		}


	}

	public File getAlignmentOutputFile(Strand strand, Sample s, Reference r) {
		return new File(this.project.getOutputDirectory() + File.separator + "bisulfited_CT_" + s.getName() +
				"_against_" + r.getReferenceFile().getName() + "_" + strand.name() + ".sam");
	}

	private interface BowtieCommandCreator {
		String[] getCommand(File reference, Sample sample, Strand strand, boolean nohead);
	}

	private class Bowtie1CommandCreator implements BowtieCommandCreator {

		private int e, l, n, chunkmbs, I, X;
		private Bowtie1Quals solexaQ;

		Bowtie1CommandCreator(/* bowtie params */
							  final int e,
							  final int l,
							  final int n,


							  final int chunkmbs,
							  final Bowtie1Quals solexaQ,

				/*bowtie paired-end parameters*/
							  final int I,
							  final int X) {
			this.e = e;
			this.l = l;
			this.n = n;
			this.chunkmbs = chunkmbs;
			this.solexaQ = solexaQ;
			this.I = I;
			this.X = X;
		}

		@Override
		public String[] getCommand(File reference, Sample sample, Strand strand, boolean nohead) {
			String[] command = null;

			String bowtiePath = (sample.getProject().getBowtieDirectory() != null ? sample.getProject()
					.getBowtieDirectory().getAbsolutePath() + File
					.separator : "") + "bowtie";
			final int M = 1; // tag if there are more than one possible alignment with XM:i:>2
			final int k = 1; // report only 1 alignment. This is mandatory in order to postprocessing works

			if (!sample.isPaired()) {
				if (nohead) {
					command = new String[]{
							bowtiePath,
							"-t",
							"--chunkmbs", "" + chunkmbs,
							"--mm",
							solexaQ.getParameterValue(),
							"-e", "" + e,
							"-l", "" + l,
							"-n", "" + n,
							"-k", "" + k,
							"-M", "" + M,
							"--best",
							"-S",
							"--sam-nohead",
							"--sam-RG", "ID:" + strand.name(),
							"--sam-RG", "SM:" + strand.name(),
							"--nomaqround",
							reference.getAbsolutePath(),
							"-"};
				} else {
					command = new String[]{
							bowtiePath,
							"-t",
							"--chunkmbs", "" + chunkmbs,
							solexaQ.getParameterValue(),
							"-e", "" + e,
							"-l", "" + l,
							"-n", "" + n,
							"-k", "" + k,
							"-M", "" + M,
							"--best",
							"-S",
							"--sam-RG", "ID:" + strand.name(),
							"--sam-RG", "SM:" + strand.name(),
							"--nomaqround",
							reference.getAbsolutePath(),
							"-"};
				}
			} else { //paired
				if (nohead) {
					command = new String[]{
							bowtiePath,
							"-t",
							"--chunkmbs", "" + chunkmbs,
							"--mm", solexaQ.getParameterValue(),
							"-e", "" + e,
							"-l", "" + l,
							"-n", "" + n,
							"-k", "" + k,
							"-M", "" + M,
							"--best",
							"-S",
							"--sam-nohead",
							"--sam-RG", "ID:" + strand.name(),
							"--sam-RG", "SM:" + strand.name(),
							"--nomaqround",
							reference.getAbsolutePath(),
							"-I", "" + I,
							"-X", "" + X,
							"--fr",
							"--12", "-"};
				} else {
					command = new String[]{
							bowtiePath,
							"-t",
							"--chunkmbs", "" + chunkmbs,
							solexaQ.getParameterValue(),
							"-e", "" + e,
							"-l", "" + l,
							"-n", "" + n,
							"-k", "" + k,
							"-M", "" + M,
							"--best",
							"-S",
							"--sam-RG", "ID:" + strand.name(),
							"--sam-RG", "SM:" + strand.name(),
							"--nomaqround", reference.getAbsolutePath(),
							"-I", "" + I,
							"-X", "" + X,
							"--fr",
							"--12", "-"};
				}
			}
			return command;
		}
	}

	private class Bowtie2CommandCreator implements BowtieCommandCreator {

		private final boolean local;
		private final String i;
		private final int R;
		private final String sm;
		private int L, N, I, X, D;
		private Bowtie2Quals quals;

		Bowtie2CommandCreator(/* bowtie params */
							  final boolean local,
							  final int D,
							  final int R,
							  final int L,
							  final String i,
							  final String sm,
							  final int N,
							  final Bowtie2Quals quals,

			/*bowtie paired-end parameters*/
							  final int I,
							  final int X) {
			this.local = local;
			this.D = D;
			this.R = R;
			this.L = L;
			this.i = i;
			this.sm = sm;
			this.N = N;
			this.quals = quals;
			this.I = I;
			this.X = X;
		}

		@Override
		public String[] getCommand(File reference, Sample sample, Strand strand, boolean nohead) {
			String[] command = null;

			String bowtiePath = (sample.getProject().getBowtie2Directory() != null ? sample.getProject()
					.getBowtie2Directory().getAbsolutePath() + File
					.separator : "") + "bowtie2";

			// Where is k option? Bowtie 2, by default reports only (the best) alignment for a read. If more than one
			// alignment is possible, the XS tag will be present. If we want to use -k option, it is mandatory that it
			// must be equal to 1. If not, multiple alignments may be reported and this is not supported by our
			// postprocessing step. We discard the use of -k by now (in fact, in bowtie2 this option is not used in
			// any of their presets)

			if (!sample.isPaired()) {
				if (nohead) {
					command = new String[]{
							bowtiePath,
							"-t",
							"--mm",
							local ? "--local" : "",
							"--no-discordant",
							"--no-mixed",
							quals.getParameterValue(),
							"-D", "" + D,
							"-R", "" + R,
							"-L", "" + L,
							"-i", "" + i,
							"-i", "" + i,
							"--score-min", "" + sm,
							"-N", "" + N,
							"--no-hd",
							"--rg-id", strand.name(),
							"--rg", "SM:" + strand.name(),
							"--sam-no-qname-trunc",
							"-x", reference.getAbsolutePath(),
							"-U", "-"};
				} else {
					command = new String[]{
							bowtiePath,
							"-t",
							local ? "--local" : "",
							"--no-discordant",
							"--no-mixed",
							quals.getParameterValue(),
							"-D", "" + D,
							"-R", "" + R,
							"-L", "" + L,
							"-i", "" + i,
							"--score-min", "" + sm,
							"-N", "" + N,
							"--rg-id", strand.name(),
							"--rg", "SM:" + strand.name(),
							"-x", reference.getAbsolutePath(),
							"-U", "-"};
				}
			} else { //paired
				if (nohead) {
					command = new String[]{
							bowtiePath, "-t",
							"--mm",
							local ? "--local" : "",
							"--no-discordant",
							"--no-mixed",
							quals.getParameterValue(),
							"-D", "" + D,
							"-R", "" + R,
							"-L", "" + L,
							"-i", "" + i,
							"--score-min", "" + sm,
							"-N", "" + N,
							"--no-hd",
							"--rg-id", strand.name(),
							"--rg", "SM:" + strand.name(),
							"--sam-no-qname-trunc",
							"-x", reference.getAbsolutePath(),
							"-I", "" + I,
							"-X", "" + X,
							"--fr",
							"--tab5",
							"-"};
				} else {
					command = new String[]{bowtiePath,
							"-t", local ? "--local" : "",
							"--no-discordant",
							"--no-mixed",
							quals.getParameterValue(),
							"-D", "" + D,
							"-R", "" + R,
							"-L", "" + L,
							"-i", "" + i,
							"--score-min", "" + sm,
							"-N", "" + N,
							"--rg-id", strand.name(),
							"--rg", "SM:" + strand.name(),
							"-x", reference.getAbsolutePath(),
							"-I", "" + I,
							"-X", "" + X,
							"--fr",
							"--tab5",
							"-"};
				}
			}
			return command;
		}
	}

	public void performBowtie1Alignment(
			final Sample sample,
			final Reference reference,
			boolean skipUnconverted,
			int threadsNumber,
			
			/* bowtie params */
			final int e,
			final int l,
			final int n,


			final int chunkmbs,
			final Bowtie1Quals solexaQ) throws IOException {
		performBowtie1Alignment(sample, reference, skipUnconverted, threadsNumber, e, l, n, chunkmbs, solexaQ, 0, 250);
	}


	public void performBowtie1Alignment(final Sample sample,
										final Reference reference,
										boolean skipUnconverted,
										int threadsNumber,

			/* bowtie params */
										final int e,
										final int l,
										final int n,


										final int chunkmbs,
										final Bowtie1Quals solexaQ,

			/*bowtie paired-end parameters*/
										final int I,
										final int X) throws IOException {

		Bowtie1CommandCreator commandCreator = new Bowtie1CommandCreator(e, l, n, chunkmbs, solexaQ, I, X);


		performBowtieAlignment(sample, reference, skipUnconverted, threadsNumber, commandCreator, new
				Bowtie1ScoreFunction());
	}

	public void performBowtie2Alignment(final Sample sample,
										final Reference reference,
										boolean skipUnconverted,
										int threadsNumber,

			/* bowtie 2 params */
										final boolean local,
										final int D,
										final int R,
										final int L,
										final String i,
										final String sm,
										final int N,
										final Bowtie2Quals quals) throws IOException {

		Bowtie2CommandCreator commandCreator = new Bowtie2CommandCreator(local, D, R, L, i, sm, N, quals, 0, 250);


		performBowtieAlignment(sample, reference, skipUnconverted, threadsNumber, commandCreator, new
				Bowtie2ScoreFunction());
	}

	public void performBowtie2Alignment(final Sample sample,
										final Reference reference,
										boolean skipUnconverted,
										int threadsNumber,

			/* bowtie 2 params */
										final boolean local,
										final int D,
										final int R,
										final int L,
										final String i,
										final String sm,
										final int N,
										final Bowtie2Quals quals,

		/*bowtie paired-end parameters*/
										final int I,
										final int X) throws IOException {

		Bowtie2CommandCreator commandCreator = new Bowtie2CommandCreator(local, D, R, L, i, sm, N, quals, I, X);

		// AS: Alignment score (bowtie 2)
		performBowtieAlignment(sample, reference, skipUnconverted, threadsNumber, commandCreator, new
				Bowtie2ScoreFunction());
	}

	private void performBowtieAlignment(
			final Sample sample,
			final Reference reference,
			boolean skipUnconverted,
			int threadsNumber,
			BowtieCommandCreator commandCreator, AlignmentScoreFunction scoreFunction) throws IOException {

		logger.info("Peforming alignment of sample " + sample.getName() + " against " + reference.getReferenceFile()
				.toString().replaceAll(project.getReferenceDirectory() + File.separator, ""));


		ReferenceBisulfitation rb = new ReferenceBisulfitation(this.project);

		final File refCT = rb.getBisulfitedReference(ReferenceBisulfitation.Replacement.CT, reference);
		final File refGA = rb.getBisulfitedReference(ReferenceBisulfitation.Replacement.GA, reference);

		if (!refCT.exists()) {
			throw new IllegalArgumentException("cannot find in-silico CtoT bifulfited reference: " + refCT + ". " +
					"Perform reference bisulfitation first.");
		}

		if (!refGA.exists()) {
			throw new IllegalArgumentException("cannot find in-silico GtoA bifulfited reference: " + refGA + ". " +
					"Perform reference bisulfitation first.");
		}

		final File alignmentOutputFileCT = getAlignmentOutputFile(Strand.WATSON, sample, reference);
		final File alignmentOutputFileGA = getAlignmentOutputFile(Strand.CRICK, sample, reference);


		int threads = threadsNumber / 2;
		if (threads == 0) threads = 1;

		List<BufferedReader> streamsWATSON = null;
		List<BufferedReader> streamsCRICK = null;


		if (!sample.isPaired()) {
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
		} else {
			//paired end

			//streamsWATSON
			streamsWATSON = new LinkedList<BufferedReader>();
			List<BufferedReader> mate1Readers = FastqSplitter.splitfastq(sample.getReadsMate1Files(), threads);
			List<BufferedReader> mate2Readers = FastqSplitter.splitfastq(sample.getReadsMate2Files(), threads);

			for (int i = 0; i < mate1Readers.size(); i++) {
				streamsWATSON.add(new PairedEndBowtieReader(sample, mate1Readers.get(i), mate2Readers.get(i), sample
						.isDirectional(), skipUnconverted));
			}

			//streamsCRICK
			streamsCRICK = new LinkedList<BufferedReader>();
			mate1Readers = FastqSplitter.splitfastq(sample.getReadsMate1Files(), threads);
			mate2Readers = FastqSplitter.splitfastq(sample.getReadsMate2Files(), threads);

			for (int i = 0; i < mate1Readers.size(); i++) {
				streamsCRICK.add(new PairedEndBowtieReader(sample, mate1Readers.get(i), mate2Readers.get(i), sample
						.isDirectional(), skipUnconverted));
			}
		}

		abstract class LineProcessor {
			public abstract void processLine(String line);
		}


		class AlignerThread extends Thread {
			File ref;
			LineProcessor out;
			private String logFileName;
			private BufferedReader readsStream;
			private boolean nohead;
			boolean shouldStop = false;
			private Strand strand;

			public AlignerThread(File ref, Strand strand, LineProcessor out, String logfileName, BufferedReader
					readsStream, boolean nohead) {
				this.ref = ref;
				this.out = out;
				this.logFileName = logfileName;
				this.readsStream = readsStream;
				this.nohead = nohead;
				this.strand = strand;
			}

			public void run() {


				logger.info("Aligning " +
						sample.getReadsFiles().toString().replaceAll(project.getReadsDirectory().toString() + File
								.separator, "")
						+ " " +
						"against " +
						"[" + ref.toString().replaceAll(project.getWorkingDirectory().toString() + File.separator, "")
						+ "]...... " +
						"(see .log file)...... ");

				String[] command = commandCreator.getCommand(ref, sample, strand, nohead);


				String outFile = logFileName;


				try {
					FileOutputStream outLog = new FileOutputStream(new File(outFile));

					final Process process = Tools.executeProcess(command, null, outLog);

					BufferedReader stdOut = new BufferedReader(new InputStreamReader(process.getInputStream()));

					String line = null;

					//thread to feed lines
					new Thread() {
						public void run() {
							PrintStream ps = new PrintStream(process.getOutputStream());

							String readsLine = null;
							logger.info("Start read feeding to alignment against " + ref.toString().replaceAll(project
									.getWorkingDirectory() + File.separator, Project.WORKING_DIRECTORY) + ". Log " +
									"file:" +
									" " +
									"" + logFileName
									.replaceAll
											(project.getOutputDirectory() + File.separator, Project.OUTPUT_DIRECTORY));
							try {
								while ((readsLine = readsStream.readLine()) != null && !shouldStop) {
									ps.println(readsLine);

								}
								ps.flush();
								ps.close();

								logger.info("Finished read feeding to alignment against " + ref.toString().replaceAll
										(project
												.getWorkingDirectory() + File.separator, Project.WORKING_DIRECTORY));
							} catch (Exception e) {
								throw new RuntimeException(e);
							}
						}

						;
					}.start();

					try {
						while ((line = stdOut.readLine()) != null) {
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

		class AlignerPostprocessor {
			private final AlignmentScoreFunction scoreFunction;
			private int id;
			private String CTLine = null;
			private String GALine = null;

			private int tagCount = 0;
			private int mergeCount = 0;

			StringBuilder outputBufferCT = new StringBuilder(100000);
			StringBuilder outputBufferGA = new StringBuilder(100000);

			private Pattern scorePattern;

			public AlignerPostprocessor(int id, AlignmentScoreFunction scoreFunction) {
				this.id = id;
				this.scoreFunction = scoreFunction;
			}

			public void close() {
				flushBuffer();
				logger.info("Both alignments have finished. Ambigous reads: " + tagCount);

			}

			// for single-end
			private String directionalPreviousLineCT;
			private boolean directionalPreviousLineCTWasAligned;
			private String directionalPreviousLineGA;
			private boolean directionalPreviousLineGAWasAligned;
			private int directionalScore = Integer.MIN_VALUE;


			// in paired-end
			private int previousScoreMate1 = Integer.MIN_VALUE;
			private int previousScoreMate2 = Integer.MIN_VALUE;
			private int nonDirectionalPreviousScoreMate1 = Integer.MIN_VALUE;
			private String directionalPreviousLineCTMate1;
			private boolean directionalPreviousLineCTMate1WasAligned;
			private String directionalPreviousLineGAMate1;
			private boolean directionalPreviousLineGAMate1WasAligned;
			private String directionalPreviousLineCTMate2;
			private boolean directionalPreviousLineCTMate2WasAligned;
			private String directionalPreviousLineGAMate2;
			private boolean directionalPreviousLineGAMate2WasAligned;
			private String nonDirectionalPreviousLineCTMate1;
			private boolean nonDirectionalPreviousLineCTMate1WasAligned;
			private String nonDirectionalPreviousLineGAMate1;
			private boolean nonDirectionalPreviousLineGAMate1WasAligned;
			// in both PE and SE
			private boolean directionalWasAmbiguous = false;
			private boolean nonDirectionalIsAmbiguous = false;


			public void merge() {

				CTLine = replaceOriginalRead(CTLine).trim();
				GALine = replaceOriginalRead(GALine).trim();

				boolean ambiguous = false;
				if (!CTLine.startsWith("@")) {
					String[] tokensCT = CTLine.split("\t");
					String[] tokensGA = GALine.split("\t");

					//System.out.println(tokensCT[0]+" = "+tokensGA[0]);
					if (!sample.isPaired() && !tokensCT[0].equals(tokensGA[0])) {
						// Note: this does not happen when bowtie says "Exhausted best-first chunk memory for read"

						logger.severe("BUG: reading two samrecords from CT and GA alignments with are a " +
								"different read	CT:" + CTLine + "\nGA:" + GALine);
						System.exit(1);
					} else if (sample.isPaired() && !tokensCT[0].substring(0, tokensCT[0].length() - 1).equals
							(tokensGA[0].substring(0, tokensGA[0].length() - 1))) {
						logger.severe("BUG: reading two samrecords from CT and GA alignments with are a different " +
								"read (ignoring last character)\nCT:" + CTLine + "\nGA:" + GALine);
						System.exit(1);
					}

					if (!tokensCT[5].equals("*") && !tokensGA[5].equals("*")) {
						ambiguous = true;
						tagCount++;
					}

					if (sample.isDirectional()) {
						mergeDirectionalSAMLines(ambiguous);
					} else {
						mergeNonDirectionalSAMRecords(ambiguous, tokensCT, tokensGA);
					}


				} else {
					// header line (starting with @)
					outputBufferCT.append(CTLine + "\n");
					outputBufferGA.append(GALine + "\n");
				}


				mergeCount++;

				if (mergeCount % 10000 == 0) {
					flushBuffer();
					logger.info("Aligner thread " + this.id + ": " + mergeCount + " reads processed");
				}

				CTLine = null;
				GALine = null;
				//System.out.println("merged!");

			}

			private void mergeNonDirectionalSAMRecords(boolean ambiguous, String[] tokensCT, String[] tokensGA) {
				//get score
				int currentScore = Integer.MIN_VALUE;
				boolean alignedInCT = false;
				boolean alignedInGA = false;
				if (!tokensCT[5].equals("*")) {
					alignedInCT = true;
					currentScore = scoreFunction.getScore(CTLine);
				}

				if (!tokensGA[5].equals("*")) {
					alignedInGA = true;
					int GAScore = scoreFunction.getScore(GALine);
					if (GAScore > currentScore) {
						currentScore = GAScore;
					}
				}

				if (sample.isPaired()) {
					mergeNonDirectionalPairedEndSAMLines(ambiguous, currentScore, alignedInCT, alignedInGA);

				} else {
					mergeNonDirectionalSingleEndSAMLines(ambiguous, currentScore, alignedInCT, alignedInGA);
				}
			}

			private void mergeNonDirectionalSingleEndSAMLines(boolean ambiguous, int currentScore, boolean alignedInCT, boolean alignedInGA) {
				if (directionalPreviousLineCT == null) {
					// reading the directional attempt
					directionalPreviousLineCT = CTLine + "\tRG:Z:" + Strand.WATSON.name();
					directionalPreviousLineCTWasAligned = alignedInCT;
					directionalPreviousLineGA = GALine + "\tRG:Z:" + Strand.CRICK.name();
					directionalPreviousLineGAWasAligned = alignedInGA;
					if (ambiguous) {
						directionalWasAmbiguous = true;
					}
					directionalScore = currentScore;
				} else {
					int nonDirectionalScore = currentScore;
					// reading the non-directional attempt

					// which of the two pairs should be written, the one with the best alignment. If none
					// has alignment put the directional unaligned sam records in the output
					if ((directionalScore >= nonDirectionalScore)) {
						// no better score or no-alignment, so print only the previous (the directional)

						// the direction tag is present only if the read was aligned
						String directionTagCT = directionalPreviousLineCTWasAligned?"\tZD:A:F":"";
						String directionTagGA = directionalPreviousLineGAWasAligned?"\tZD:A:F":"";

						if (directionalWasAmbiguous || nonDirectionalScore > Integer.MIN_VALUE) {
							// the directional alignment is ambiguous when it was aligned in WATSON and
							// CRICK, but also if the non-directional attempt (currentScore) was also
							// aligned (nonDirectionalScore > MIN_VALUE)
							String ambiguousTypeTag;

							if (directionalWasAmbiguous && nonDirectionalScore > Integer.MIN_VALUE) {
								ambiguousTypeTag = "ZT:A:B";
							} else if (directionalWasAmbiguous) {
								ambiguousTypeTag = "ZT:A:S";
							} else {
								ambiguousTypeTag = "ZT:A:D";
							}

							outputBufferCT.append(directionalPreviousLineCT + directionTagCT + "\tZA:A:Y\t" +
									ambiguousTypeTag
									+ "\n");
							outputBufferGA.append(directionalPreviousLineGA + directionTagGA + "\tZA:A:Y\t" +
									ambiguousTypeTag + "\n");
						} else {
							outputBufferCT.append(directionalPreviousLineCT + directionTagCT + "\n");
							outputBufferGA.append(directionalPreviousLineGA + directionTagGA + "\n");
						}

					} else {
						// we have better score, ignore previous, print the non-directional instead

						// the direction tag is present only if the read was aligned
						String directionTagCT = alignedInCT?"\tZD:A:R":"";
						String directionTagGA = alignedInGA?"\tZD:A:R":"";

						if (ambiguous || directionalScore > Integer.MIN_VALUE) {
							// if the non-directional is ambiguous (CRICK/WATSON) or the directional also
							// was aligned, mark this alignment as ambiguous
							// ambiguity means both WASTON/CRICK or DIRECTIONAL/NON-DIRECTIONAL confusion
							nonDirectionalIsAmbiguous = true;

							String ambiguousTypeTag;
							if (ambiguous && directionalScore > Integer.MIN_VALUE) {
								ambiguousTypeTag = "ZT:A:B";
							} else if (ambiguous) {
								ambiguousTypeTag = "ZT:A:S";
							} else {
								ambiguousTypeTag = "ZT:A:D";
							}
							outputBufferCT.append(CTLine + directionTagCT + "\tZA:A:Y" + "\t" + ambiguousTypeTag +
									"\tRG:Z:" + Strand.WATSON
									.name() +
									"\n");
							outputBufferGA.append(GALine + directionTagGA + "\tZA:A:Y" + "\t" + ambiguousTypeTag +
									"\tRG:Z:" + Strand.CRICK
									.name() +
									"\n");
						} else {
							outputBufferCT.append(CTLine + directionTagCT + "\tRG:Z:" + Strand.WATSON.name() + "\n");
							outputBufferGA.append(GALine + directionTagGA + "\tRG:Z:" + Strand.CRICK.name() + "\n");
						}
					}

					directionalPreviousLineCT = null;
					directionalPreviousLineCTWasAligned = false;
					directionalPreviousLineGA = null;
					directionalPreviousLineGAWasAligned = false;
					directionalScore = Integer.MIN_VALUE;
					directionalWasAmbiguous = false;
					nonDirectionalIsAmbiguous = false;

				}
			}

			private void mergeNonDirectionalPairedEndSAMLines(boolean ambiguous, int currentScore, boolean alignedInCT, boolean alignedInGA) {
				if (directionalPreviousLineCTMate1 == null) {
					// reading the first mate alignments of the directional attempt

					directionalPreviousLineCTMate1 = CTLine + "\tRG:Z:" + Strand.WATSON.name();
					directionalPreviousLineCTMate1WasAligned = alignedInCT;
					directionalPreviousLineGAMate1 = GALine + "\tRG:Z:" + Strand.CRICK.name();
					directionalPreviousLineGAMate1WasAligned = alignedInGA;
					previousScoreMate1 = currentScore;
				} else if (directionalPreviousLineCTMate2 == null) {
					// reading the second mate alignments of the directional attempt
					directionalPreviousLineCTMate2 = CTLine + "\tRG:Z:" + Strand.WATSON.name();
					directionalPreviousLineCTMate2WasAligned = alignedInCT;
					directionalPreviousLineGAMate2 = GALine + "\tRG:Z:" + Strand.CRICK.name();
					directionalPreviousLineGAMate2WasAligned = alignedInGA;
					previousScoreMate2 = currentScore;
					if (ambiguous) {
						directionalWasAmbiguous = true;
					}
				} else if (nonDirectionalPreviousLineCTMate1 == null) {
					// reading the first mate alignments of the non-directional attempt
					nonDirectionalPreviousLineCTMate1 = CTLine + "\tRG:Z:" + Strand.WATSON.name();
					nonDirectionalPreviousLineCTMate1WasAligned = alignedInCT;
					nonDirectionalPreviousLineGAMate1 = GALine + "\tRG:Z:" + Strand.CRICK.name();
					nonDirectionalPreviousLineGAMate1WasAligned = alignedInGA;
					nonDirectionalPreviousScoreMate1 = currentScore;
				} else {
					// reading the second mate alignments of the non-directional attempt

					//which mate to print?? those with better score
					int directionalScore = max(previousScoreMate1, previousScoreMate2);
					int nonDirectionalScore = max(nonDirectionalPreviousScoreMate1, currentScore);

					if (directionalScore >= nonDirectionalScore) {
						// the directional is better

						// the direction tag is present only if the read was aligned
						String directionTagCT = directionalPreviousLineCTMate1WasAligned?"\tZD:A:F":"";
						String directionTagGA = directionalPreviousLineGAMate1WasAligned?"\tZD:A:F":"";

						if (directionalWasAmbiguous || nonDirectionalScore > Integer.MIN_VALUE) {
							// the directional alignment is ambiguous when it was aligned in WATSON and
							// CRICK, but also if the non-directional attempt (currentScore) was also
							// aligned (nonDirectionalScore > MIN_VALUE)

							String ambiguousTypeTag;
							if (directionalWasAmbiguous && nonDirectionalScore > Integer.MIN_VALUE) {
								ambiguousTypeTag = "ZT:A:B";
							} else if (directionalWasAmbiguous) {
								ambiguousTypeTag = "ZT:A:S";
							} else {
								ambiguousTypeTag = "ZT:A:D";
							}
							outputBufferCT.append(directionalPreviousLineCTMate1 + directionTagCT + "\tZA:A:Y\t" +
									ambiguousTypeTag + "\n");
							outputBufferGA.append(directionalPreviousLineGAMate1 + directionTagGA + "\tZA:A:Y\t" +
									ambiguousTypeTag + "\n");
							outputBufferCT.append(directionalPreviousLineCTMate2 + directionTagCT + "\tZA:A:Y\t" +
									ambiguousTypeTag + "\n");
							outputBufferGA.append(directionalPreviousLineGAMate2 + directionTagGA + "\tZA:A:Y\t" +
									ambiguousTypeTag + "\n");
						} else {
							outputBufferCT.append(directionalPreviousLineCTMate1 + directionTagCT + "\n");
							outputBufferGA.append(directionalPreviousLineGAMate1 + directionTagGA + "\n");
							outputBufferCT.append(directionalPreviousLineCTMate2 + directionTagCT + "\n");
							outputBufferGA.append(directionalPreviousLineGAMate2 + directionTagGA + "\n");
						}
					} else {
						//we have better score, the non-directional is better

						// the direction tag is present only if the read was aligned
						String directionTagCT = nonDirectionalPreviousLineCTMate1WasAligned?"\tZD:A:R":"";
						String directionTagGA = nonDirectionalPreviousLineGAMate1WasAligned?"\tZD:A:R":"";

						if (ambiguous || directionalScore > Integer.MIN_VALUE) {
							// if the non-directional is ambiguous (CRICK/WATSON) or the directional also
							// was aligned, mark this alignment as ambiguous
							// ambiguity means both WASTON/CRICK or DIRECTIONAL/NON-DIRECTIONAL confusion
							String ambiguousTypeTag;
							if (ambiguous && directionalScore > Integer.MIN_VALUE) {
								ambiguousTypeTag = "ZT:A:B";
							} else if (ambiguous) {
								ambiguousTypeTag = "ZT:A:S";
							} else {
								ambiguousTypeTag = "ZT:A:D";
							}
							//print first mate
							outputBufferCT.append(nonDirectionalPreviousLineCTMate1 + directionTagCT + "\tZA:A:Y" +
									"\t" + ambiguousTypeTag + "\n");
							outputBufferGA.append(nonDirectionalPreviousLineGAMate1 + directionTagGA + "\tZA:A:Y" +
									"\t" +
									ambiguousTypeTag + "\n");

							// print second mate (current lines)
							outputBufferCT.append(CTLine + "\tRG:Z:" + Strand.WATSON.name() + directionTagCT + "\tZA:A:Y"
									+ "\t" + ambiguousTypeTag + "\n");
							outputBufferGA.append(GALine + "\tRG:Z:" + Strand.CRICK.name() + directionTagGA +
									"\tZA:A:Y" +
									"\t" + ambiguousTypeTag + "\n");
						} else {
							//print first mate
							outputBufferCT.append(nonDirectionalPreviousLineCTMate1 + directionTagCT + "\n");
							outputBufferGA.append(nonDirectionalPreviousLineGAMate1 + directionTagGA + "\n");

							// print second mate (current lines)
							outputBufferCT.append(CTLine + directionTagCT + "\tRG:Z:" + Strand.WATSON.name() + "\n");
							outputBufferGA.append(GALine + directionTagGA + "\tRG:Z:" + Strand.CRICK.name() + "\n");
						}

					}
					previousScoreMate1 = Integer.MIN_VALUE;
					previousScoreMate2 = Integer.MIN_VALUE;
					nonDirectionalPreviousScoreMate1 = Integer.MIN_VALUE;
					directionalPreviousLineCTMate1 = null;
					directionalPreviousLineCTMate1WasAligned = false;
					directionalPreviousLineGAMate1 = null;
					directionalPreviousLineGAMate1WasAligned = false;
					directionalPreviousLineCTMate2 = null;
					directionalPreviousLineCTMate2WasAligned = false;
					directionalPreviousLineGAMate2 = null;
					directionalPreviousLineGAMate2WasAligned = false;
					nonDirectionalPreviousLineCTMate1 = null;
					nonDirectionalPreviousLineCTMate1WasAligned = false;
					nonDirectionalPreviousLineGAMate1 = null;
					nonDirectionalPreviousLineGAMate1WasAligned = false;
					directionalWasAmbiguous = false;
				}
			}

			private void mergeDirectionalSAMLines(boolean ambiguous) {
				if (ambiguous) {
					outputBufferCT.append(CTLine + "\tZA:A:Y\tRG:Z:" + Strand.WATSON.name() + "\n");
					outputBufferGA.append(GALine + "\tZA:A:Y\tRG:Z:" + Strand.CRICK.name() + "\n");
				} else {
					outputBufferCT.append(CTLine + "\tRG:Z:" + Strand.WATSON.name() + "\n");
					outputBufferGA.append(GALine + "\tRG:Z:" + Strand.CRICK.name() + "\n");
				}
			}

			private void flushBuffer() {
				synchronized (outCT) {
					outCT.print(outputBufferCT.toString());
					outCT.flush();
					outputBufferCT.setLength(0);

					outGA.print(outputBufferGA.toString());
					outGA.flush();
					outputBufferGA.setLength(0);
				}

			}

			public LineProcessor CTProcessor = new LineProcessor() {

				@Override
				public void processLine(String line) {
					synchronized (AlignerPostprocessor.this) {
						while (CTLine != null) {
							try {
								AlignerPostprocessor.this.wait(10000);
								if (CTLine != null) {
									logger.info("awaking, but CTLine is still not null (if you see this message " +
											"continously, bowtie may be not responding), it is: " + CTLine +
											"\nProcessing new CT line: " + line);
								}
							} catch (InterruptedException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
						}
						CTLine = line;

						if (GALine != null) {
							merge();
							AlignerPostprocessor.this.notify();
						}
					}

				}

			};

			public LineProcessor GAProcessor = new LineProcessor() {

				@Override
				public void processLine(String line) {
					synchronized (AlignerPostprocessor.this) {
						while (GALine != null) {
							try {
								AlignerPostprocessor.this.wait(10000);
								if (GALine != null) {
									logger.info("awaking, but GALine is still not null, it is: " + GALine +
											"\nProcessing new GA line: " + line);
								}
							} catch (InterruptedException e) {

								e.printStackTrace();
							}
						}
						GALine = line;

						if (CTLine != null) {
							merge();
							AlignerPostprocessor.this.notify();
						}
					}

				}

			};
		}

		//write the header of the sam doing a "dummy alignment"
		AlignerPostprocessor dummyposprocessor = new AlignerPostprocessor(0, scoreFunction);
		LineProcessor dummyctProcessor = dummyposprocessor.CTProcessor;
		LineProcessor dummygaProcessor = dummyposprocessor.GAProcessor;
		AlignerThread dummyThreadCT = new AlignerThread(refCT, Strand.WATSON, dummyctProcessor, alignmentOutputFileCT
				+ "_p_head.log", new BufferedReader(new InputStreamReader(new ByteArrayInputStream(new byte[0]))),
				false);
		AlignerThread dummyThreadGA = new AlignerThread(refGA, Strand.CRICK, dummygaProcessor, alignmentOutputFileGA +
				"_p_head.log", new BufferedReader(new InputStreamReader(new ByteArrayInputStream(new byte[0]))),
				false);
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

		for (int i = 0; i < streamsWATSON.size(); i++) {
			AlignerPostprocessor postprocessor = new AlignerPostprocessor(i + 1, scoreFunction);
			postprocessors.add(postprocessor);
			LineProcessor ctProcessor = postprocessor.CTProcessor;
			LineProcessor gaProcessor = postprocessor.GAProcessor;
			BufferedReader streamCT = streamsWATSON.get(i);
			BufferedReader streamGA = streamsCRICK.get(i);
			AlignerThread threadCT = new AlignerThread(refCT, Strand.WATSON, ctProcessor, alignmentOutputFileCT +
					"_p_" + i + ".log", streamCT, true);
			AlignerThread threadGA = new AlignerThread(refGA, Strand.CRICK, gaProcessor, alignmentOutputFileGA + "_p_"
					+ i + ".log", streamGA, true);

			threadCT.start();
			threadGA.start();

			alignerThreads.add(threadCT);
			alignerThreads.add(threadGA);

		}

		for (Thread t : alignerThreads) {
			try {
				t.join();

			} catch (InterruptedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		}
		for (AlignerPostprocessor postprocessor : postprocessors) {
			postprocessor.close();
		}
		outCT.flush();
		outGA.flush();
		outCT.close();
		outGA.close();

		logger.info("Alignment of sample " + sample.getName() + " OK");

	}

	private String getReverseComplementary(String sequence) {
		StringBuilder complementaria = new StringBuilder("");

		// 1ero construyo la complementaria
		for (int i = 0; i < sequence.length(); i++) {
			switch (sequence.charAt(i)) {
				case 'A':
					complementaria.append("T");
					break;
				case 'T':
					complementaria.append("A");
					break;
				case 'C':
					complementaria.append("G");
					break;
				case 'G':
					complementaria.append("C");
					break;
				default:
					complementaria.append(sequence.charAt(i));
			}
		}

		//ahora la reversa
		StringBuilder reversa = new StringBuilder("");
		for (int i = complementaria.length() - 1; i > -1; i--) {
			reversa.append(complementaria.charAt(i));
		}

		return (reversa.toString());
	}

	private String replaceOriginalRead(String samline) {
		if (!samline.startsWith("@")) {
			//System.out.println(samline);
			final String[] tokens = samline.split("[\t]");
			// recupero la secuencia inicial
			final String[] firstColumn = tokens[0].split("[|][|]");
			String originalRead = null;

			int flag = Integer.parseInt(tokens[1]);
			boolean mate1 = false;
			boolean paired = false;
			if ((flag & 0x0001) == 0x0001) {
				//paired!
				paired = true;
				if ((flag & 0x0040) == 0x0040) {
					//mate1
					mate1 = true;
					originalRead = firstColumn[1];
				} else if ((flag & 0x0080) == 0x0080) {
					//mate2
					mate1 = false;
					originalRead = firstColumn[2];
				} else {
					throw new RuntimeException("Malformed FLAG in SAM. It says that is a paired read, but it is not " +
							"the first nor the second pair");
				}
			} else {
				originalRead = firstColumn[1];
			}
			final StringBuilder lineModified = new StringBuilder();

			if ((flag & 0x0010) == 0x0010) {
				originalRead = getReverseComplementary(originalRead);
			}


			for (int j = 0; j < tokens.length; j++) {
				// le adjunto la cabecera original de la read
				if (j == 0) {
					lineModified.append(firstColumn[0]);
					if (paired) {
						if (mate1) {
							lineModified.append("/1");
						} else {
							lineModified.append("/2");
						}
					}
					lineModified.append("\t");
				}
				// le adjunto la read original en lugar de la que tenia
				else if (j == 9) lineModified.append(originalRead).append("\t");

					//else if(j==tokens.length-1) lineModified=new StringBuilder(lineModified).append(tokens[j]);
				else lineModified.append(tokens[j]).append("\t");

			}
			return lineModified.toString();
		}//if(!thisLine.startsWith("@"))
		return samline;
	}


}


