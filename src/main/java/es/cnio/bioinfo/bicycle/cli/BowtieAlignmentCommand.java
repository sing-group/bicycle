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
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment.Quals;

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
		int e = Integer.parseInt(parameters.get(this.findOption("e")));
		int l = Integer.parseInt(parameters.get(this.findOption("l")));
		int n = Integer.parseInt(parameters.get(this.findOption("n")));
		//int k = Integer.parseInt(parameters.get(this.findOption("k")));
		//int m = Integer.parseInt(parameters.get(this.findOption("m")));		
		int c = Integer.parseInt(parameters.get(this.findOption("c")));
		int I = Integer.parseInt(parameters.get(this.findOption("I")));
		int X = Integer.parseInt(parameters.get(this.findOption("X")));

		boolean skipUnconverted = parameters.containsKey(this.findOption("u"));

		BowtieAlignment.PairedEndOrientation orientation = BowtieAlignment.PairedEndOrientation.FR;

		/*if (parameters.containsKey(this.findOption("rf"))){
			orientation = BowtieAlignment.PairedEndOrientation.RF;
		}else if (parameters.containsKey(this.findOption("ff"))){
			orientation = BowtieAlignment.PairedEndOrientation.FF;
		}
		*/
		String qualsString = parameters.get(this.findOption("q"));
		Quals quals = null;
		if (qualsString.equalsIgnoreCase("solexa1.3")) {
			quals = Quals.AFTER_1_3;
		} else if (qualsString.equalsIgnoreCase("solexa")) {
			quals = Quals.BEFORE_1_3;
		} else if (qualsString.equalsIgnoreCase("phred33")) {
			quals = Quals.PHRED_33;
		}
		if (quals == null) {
			throw new IllegalArgumentException("invalid qualities parameter: " + qualsString);
		}

		BowtieAlignment ba = new BowtieAlignment(project);
		for (Sample sample : project.getSamples()) {
			for (Reference reference : project.getReferences()) {
				ba.performBowtieAlignment(sample, reference, skipUnconverted, t, e, l, n, c, quals, I, X);
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
		/*toret.add(new Option("bowtie-fr", "fr", 
				"The upstream/downstream mate orientations for a valid paired-end alignment against the forward
				reference strand: forward-reverse (this is the default)", true, false));
		toret.add(new Option("bowtie-fr", "rf", 
				"The upstream/downstream mate orientations for a valid paired-end alignment against the forward
				reference strand: reverse-forward", true, false));
		toret.add(new Option("bowtie-fr", "ff", 
				"The upstream/downstream mate orientations for a valid paired-end alignment against the forward
				reference strand: forward-forward", true, false));
		*/
		
		/*toret.add(new DefaultValuedOption("bowtie-k", "k", 
				"Report up to <k> valid alignments per read or pair", "10"));
		*/
		/*toret.add(new DefaultValuedOption("bowtie-m", "m", 
				"Suppress all alignments for a particular read or pair if more than <m> reportable alignments exist
				for it", "1"));
		*/
		toret.add(new DefaultValuedOption("bowtie-chunkmbs", "c",
				"The number of megabytes of memory a given thread is given to store path descriptors", "64"));

		toret.add(new DefaultValuedOption("bowtie-quals", "q",
				"How qualities will be treated. Valid values are: solexa1.3, solexa, phred33", "solexa1.3"));


		return toret;
	}

}
