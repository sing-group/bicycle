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

package es.cnio.bioinfo.bicycle.cli;

import java.io.File;
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
		boolean onlyWithOneAlignment = parameters.containsKey(this.findOption("o"));
		boolean removeClonal = parameters.containsKey(this.findOption("c"));
		boolean removeSAM = parameters.containsKey(this.findOption("m"));
		boolean correctNonCG = parameters.containsKey(this.findOption("g"));
		int trimuntil = Integer.parseInt(parameters.get(this.findOption("t")));
		boolean trimreads = trimuntil != 0;
		int mindepth = Integer.parseInt(parameters.get(this.findOption("d")));

		double fdr = Double.parseDouble(parameters.get(this.findOption("f")));

		List<File> bedFiles = new LinkedList<File>();

		if (parameters.containsKey(this.findOption("b"))) {
			String bedFilesString = parameters.get(this.findOption("b"));
			if (bedFilesString != null) {
				String[] tokens = bedFilesString.split(",");

				for (String token : tokens) {
					File bedFile = new File(token);
					if (!bedFile.exists()) {
						throw new IllegalArgumentException("BED file not found: " + bedFile);
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
		if (errorMode == ErrorRateMode.from_control_genome) {

			if (errorTokens.length < 2) {
				throw new IllegalArgumentException("control genome must be set");
			}

			String controlGenome = errorTokens[1];

			for (Sample sample : project.getSamples()) {
				for (Reference reference : project.getReferences()) {
					ma.analyzeWithErrorFromControlGenome(reference, sample, trimreads, trimuntil, removeAmbiguous,
							onlyWithOneAlignment, removeBad, removeClonal, correctNonCG, mindepth, fdr, nThreads,
							bedFiles, controlGenome, removeSAM);
				}
			}

		} else if (errorMode == ErrorRateMode.FIXED) {
			if (errorTokens.length < 2) {
				throw new IllegalArgumentException("error rates must be set");
			}
			String[] errorRates = errorTokens[1].split(",");
			if (errorRates.length != 2) {
				throw new IllegalArgumentException("bad error rates. It must be <watson_rate>,<crick_rate>. E.g.: " +
						"0.01,0.02");
			}

			double watsonError = Double.parseDouble(errorRates[0]);
			double crickError = Double.parseDouble(errorRates[1]);

			for (Sample sample : project.getSamples()) {
				for (Reference reference : project.getReferences()) {
					ma.analyzeWithFixedErrorRate(reference, sample, trimreads, trimuntil, removeAmbiguous,
							onlyWithOneAlignment, removeBad,
							removeClonal, correctNonCG, mindepth, fdr, nThreads, bedFiles, watsonError, crickError, removeSAM);
				}
			}
		} else if (errorMode == ErrorRateMode.from_barcodes) {
			for (Sample sample : project.getSamples()) {
				for (Reference reference : project.getReferences()) {
					ma.analyzeWithErrorFromBarcodes(reference, sample, trimreads, trimuntil, removeAmbiguous,
							onlyWithOneAlignment, removeBad, removeClonal, correctNonCG, mindepth, fdr, nThreads,
							bedFiles, removeSAM);
				}
			}
		}
	}

	@Override
	protected List<Option> createOptions() {
		List<Option> toret = super.createOptions();

		toret.add(new DefaultValuedOption("threads", "n",
				"number of threads to analyze", "4"));

		toret.add(new Option("remove-uncorrectly-converted", "r",
				"ignore non-correctly bisulfite-converted reads", true, false));

		toret.add(new Option("remove-ambiguous", "a",
				"ignore reads aligned to both Watson and Crick strands", true, false));
		toret.add(new Option("only-with-one-alignment", "o",
				"ignore reads with more than one possible alignment", true, false));
		toret.add(new DefaultValuedOption("trim-reads", "t",
				"Trim reads to the <t> mismatch. 0 means no trim", "4"));

		toret.add(new DefaultValuedOption("min-depth", "d",
				"Ignore positions with less than <d> reads", "1"));

		toret.add(new DefaultValuedOption("fdr", "f",
				"FDR threshold", "0.01"));

		toret.add(new DefaultValuedOption("error-mode", "e",
				"Error rate computation mode. Valid options are: " + ErrorRateMode.from_control_genome +
						"=<control_genome_name>, " + ErrorRateMode.from_barcodes + ", " + ErrorRateMode.FIXED.name() +
						"=<watson_error_rate,crick_error_rate>", ErrorRateMode.FIXED.name() + "=0.01,0.01"));

		toret.add(new Option("annotate-beds", "b",
				"Comma-separated (with no spaces) list of BED files to annotate cytosines", true, true));

		toret.add(new Option("remove-clonal", "c",
				"Remove clonal reads", true, false));

		toret.add(new Option("correct non-CG to CG", "g",
				"Correct non-CG", true, false));
		
		toret.add(new Option("remove-sam", "m",
      "Remove SAM files. Only keep generated bams", true, false));

		return toret;
	}

}
