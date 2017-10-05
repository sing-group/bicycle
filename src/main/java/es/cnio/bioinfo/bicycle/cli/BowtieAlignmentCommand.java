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

import java.util.List;
import java.util.Map;

import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Reference;
import es.cnio.bioinfo.bicycle.Sample;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment.Bowtie1Quals;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment.Bowtie2Quals;

public class BowtieAlignmentCommand extends ProjectCommand {

	@Override
	public String getName() {
		return "align";
	}

	@Override
	public String getDescription() {
		return "Aligns with Bowtie against both references using multiple bowties (CtoT and GtoA)";
	}

	@Override
	public void executeImpl(CLIApplication app, Project project, Map<Option, String> parameters) throws Exception {


		int t = Integer.parseInt(parameters.get(this.findOption("t")));
		int v = Integer.parseInt(parameters.get(this.findOption("v")));
		if (v != 1 && v != 2) {
			throw new IllegalArgumentException("bowtie version must be 1 or 2");
		}

		boolean skipUnconverted = parameters.containsKey(this.findOption("u"));

		final BowtieAlignment ba = new BowtieAlignment(project);
		for (Sample sample : project.getSamples()) {
			for (Reference reference : project.getReferences()) {
				if (v == 1) { /* bowtie 1 */
					int e = Integer.parseInt(parameters.get(this.findOption("e")));
					int l = Integer.parseInt(parameters.get(this.findOption("l")));
					int n = Integer.parseInt(parameters.get(this.findOption("n")));
					int c = Integer.parseInt(parameters.get(this.findOption("c")));
					int I = Integer.parseInt(parameters.get(this.findOption("I")));
					int X = Integer.parseInt(parameters.get(this.findOption("X")));
					Bowtie1Quals quals = Bowtie1Quals.parseQuals(parameters.get(this.findOption("q")));
					ba.performBowtie1Alignment(sample, reference, skipUnconverted, t, e, l, n, c, quals, I, X);
				} else { /* bowtie 2 */
					boolean local = parameters.containsKey(this.findOption("o"));
					int D = Integer.parseInt(parameters.get(this.findOption("D")));
					int R = Integer.parseInt(parameters.get(this.findOption("R")));

					int L = Integer.parseInt(parameters.get(this.findOption("L2")));
					int N = Integer.parseInt(parameters.get(this.findOption("N2")));
					int I = Integer.parseInt(parameters.get(this.findOption("I2")));
					int X = Integer.parseInt(parameters.get(this.findOption("X2")));
					String scoreMinFunction = (!local) ? "L,-0.6,-0.6" : "G,20,8";
					if (parameters.containsKey(findOption("sm"))) {
						scoreMinFunction = parameters.get(findOption("sm"));
					}

					String i = (!local) ? "S,1,1.15" : "S,1,0.75";
					if (parameters.containsKey(findOption("f"))) {
						i = parameters.get(findOption("f"));
					}

					Bowtie2Quals quals = Bowtie2Quals.parseQuals(parameters.get(this.findOption("q2")));
					ba.performBowtie2Alignment(sample, reference, skipUnconverted, t, local, D, R, L, i,
							scoreMinFunction,
							N, quals, I,
							X);
				}
			}
		}
	}

	@Override
	protected List<Option> createOptions() {
		List<Option> toret = super.createOptions();

		toret.add(new DefaultValuedOption("threads", "t",
				"number of threads per sample and ref alignment", "4"));

		toret.add(new Option("skip-unconverted-barcodes", "b",
				"skip reads with unconverted barcodes. The barcode should be on the name of the read files and " +
						"delimited by '_' and '-'. For example, a valid reads file name would be: " +
						"AML_s_8_TGtATT-reads" +
						".fastq, so the barcode is TGtATT", true, false));

		toret.add(new DefaultValuedOption("bowtie-version", "v",
				"bowtie version to use (valid options are 1 or 2)"
				, "2"));
		toret.add(new DefaultValuedOption("bowtie-maqerr", "e",
				"Maximum permitted total of quality values at all mismatched read positions throughout the entire " +
						"alignment, not just in the \"seed\"", "140"));

		toret.add(new DefaultValuedOption("bowtie-seedlen", "l",
				"The \"seed length\"; i.e., the number of bases on the high-quality end of the read to which the -n " +
						"ceiling applies. The lowest permitted setting is 5 and the default is 28. bowtie is faster " +
						"for larger values of -l.", "20"));

		toret.add(new DefaultValuedOption("bowtie-seedmms", "n",
				"Maximum number of mismatches permitted in the \"seed\", i.e. the first L base pairs of the read " +
						"(where L is set with -l/--bowtie-seedlen). This may be 0, 1, 2 or 3", "0"));

		toret.add(new DefaultValuedOption("bowtie-I", "I",
				"The minimum insert size for valid paired-end alignments (paired-end projects only)", "0"));
		toret.add(new DefaultValuedOption("bowtie-X", "X",
				"The maximum insert size for valid paired-end alignments (paired-end projects only)", "250"));		
		toret.add(new DefaultValuedOption("bowtie-chunkmbs", "c",
				"The number of megabytes of memory a given thread is given to store path descriptors", "64"));

		toret.add(new DefaultValuedOption("bowtie-quals", "q",
				"How qualities will be treated. Valid values are: solexa1.3, solexa, phred33, phred64, integer",
				"solexa1.3"));

		/* bowtie 2 options */
		toret.add(new Option("bowtie2-local", "o",
				"Enables --local mode (by default the --end-to-end mode is used). In this mode, Bowtie 2 does not " +
						"require that the entire " +
						"read align from one end to the other. Rather," +
						" some characters may be omitted (\"soft clipped\") from the ends in order to achieve the " +
						"greatest possible alignment score", true, false));
		toret.add(new DefaultValuedOption("bowtie2-D", "D",
				"How many consecutive seed extension attempts can \"fail\" before Bowtie 2 moves on, using the " +
						"alignments found so far. A seed extension \"fails\" if it does not yield a new best or a new" +
						" " +
						"second-best alignment", "15"));
		toret.add(new DefaultValuedOption("bowtie2-R", "R",
				"Maximum number of times Bowtie 2 will \"re-seed\" reads with repetitive seeds. When " +
						"\"re-seeding,\" Bowtie 2 simply chooses a new set of reads (same length, same number of " +
						"mismatches allowed) at different offsets and searches for more alignments. A read is " +
						"considered to have repetitive seeds if the total number of seed hits divided by the number " +
						"of" +
						" seeds that aligned at least once is greater than 300", "2"));
		toret.add(new DefaultValuedOption("bowtie2-L", "L2",
				"Sets the length of the seed substrings to align during multiseed alignment. Smaller values make " +
						"alignment slower but more sensitive", "20"));
		toret.add(new Option("bowtie2-i-func", "f",
				"Sets a function governing the interval between seed substrings to use during multiseed alignment. " +
						"See bowtie2 manual for details. The default in --end-to-end mode is S,1,1.15 " +
						"and S,1,0.75 in --local mode", true, true));
		toret.add(new DefaultValuedOption("bowtie2-N", "N2",
				"Sets the number of mismatches to allowed in a seed alignment during multiseed alignment." +
						" Can be set to 0 or 1. Setting this higher makes alignment slower (often much slower) " +
						"but increases sensitivity. Default: 0", "0"));
		toret.add(new DefaultValuedOption("bowtie2-I", "I2",
				"The minimum fragment length for valid paired-end alignments (paired-end projects only)", "0"));
		toret.add(new DefaultValuedOption("bowtie2-X", "X2",
				"The maximum fragment length for valid paired-end alignments (paired-end projects only)", "500"));
		toret.add(new DefaultValuedOption("bowtie2-quals", "q2",
				"How qualities will be treated. Valid values are: solexa, phred33, phred64, int", "phred64"));
		toret.add(new Option("bowtie2-score-min", "sm",
				"Sets a function governing the minimum alignment score needed for an alignment to be considered " +
						"\"valid\" (i.e. good enough to report). This is a function of read length. For instance, " +
						"specifying L,0,-0.6 sets the minimum-score function f to f(x) = 0 + -0.6 * x, where x is the" +
						" read length. See bowtie2 manual for details. The default in --end-to-end mode is L," +
						"-0.6,-0.6 and in --local mode is G,20,8.", true, true));
		return toret;
	}

}
