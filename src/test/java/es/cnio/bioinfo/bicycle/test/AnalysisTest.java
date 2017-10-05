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

package es.cnio.bioinfo.bicycle.test;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Reference;
import es.cnio.bioinfo.bicycle.Sample;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment.Bowtie1Quals;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment.Strand;
import es.cnio.bioinfo.bicycle.operations.MethylationAnalysis;
import es.cnio.bioinfo.bicycle.operations.ReferenceBisulfitation;
import es.cnio.bioinfo.bicycle.operations.ReferenceBisulfitation.Replacement;

@RunWith(Parameterized.class)
public class AnalysisTest {
	private final int bowtieVersion;
	private final boolean bowtie2Local;

	@Parameters
	public static Collection<Object[]> data() {
		return Arrays.asList(new Object[][]{{1, false}, {2,false}, {2, true}});
	}

	public AnalysisTest(int bowtieVersion, boolean bowtie2Local)
	{
		this.bowtieVersion = bowtieVersion;
		this.bowtie2Local = bowtie2Local;
	}

	private Project prepareProject() throws IOException {

		File tempDir = Utils.generateTempDirName("newproject");
		System.err.println("CREATED PROJECT IN " + tempDir);
		Project p = Project.buildNewProject(
				tempDir,
				new File(Utils.getReferenceDirectory()),
				new File(Utils.getReadsDirectory()),
				new File(Utils.getBowtiePath()),
				new File(Utils.getBowtie2Path()),
				new File(Utils.getSamtoolsPath()),
				true);

		ReferenceBisulfitation rb = new ReferenceBisulfitation(p);
		BowtieAlignment ba = new BowtieAlignment(p);

		for (Reference ref : p.getReferences()) {
			rb.computeReferenceBisulfitation(Replacement.CT, ref, true);
			rb.computeReferenceBisulfitation(Replacement.GA, ref, true);
			if (this.bowtieVersion == 1) {
				ba.buildBowtieIndex(ref);
			} else {
				ba.buildBowtie2Index(ref);
			}
		}
		for (Sample sample : p.getSamples()) {
			/*SampleBisulfitation sb = new SampleBisulfitation(sample);
			sb.computeSampleBisulfitation(true);*/
			for (Reference reference : p.getReferences()) {
				if (bowtieVersion == 1) {
					ba.performBowtie1Alignment(sample, reference, false, 4, 140, 20, 0, 64, Bowtie1Quals.BEFORE_1_3);
				} else {
					ba.performBowtie2Alignment(sample, reference, false, 4, this.bowtie2Local, 15, 2, 20,
							!this.bowtie2Local?"S,1,1.15":"S,1,0.75",
							!this.bowtie2Local?"L,-0.6,-0.6":"G,20,8",
							0, BowtieAlignment.Bowtie2Quals.BEFORE_1_3);
				}
			}
		}

		return p;
	}

	@Test
	public void analysisMultithread() throws IOException, InterruptedException {
		Project project = prepareProject();
		try {

			MethylationAnalysis ma = new MethylationAnalysis(project);

			List<File> bedFiles = Arrays.asList(new File(Utils.getBedsDirectory()).listFiles(new FilenameFilter() {

				@Override
				public boolean accept(File arg0, String arg1) {
					return arg1.endsWith(".bed");
				}
			}));

			for (Sample sample : project.getSamples()) {
				for (Reference reference : project.getReferences()) {

					ma.analyzeWithErrorFromControlGenome(
							reference,
							sample,
							true,
							4,
							true,
							false,
							true,
							false,
							true,
							1,
							0.01,
							4,
							bedFiles,
							"control");
					assertTrue(ma.getSummaryFile(reference, sample).exists());

					System.err.println("====METHYLATION-WATSON=====");
					System.err.println(Utils.readFile(ma.getMethylationFile(Strand.WATSON, reference, sample)));
					assertTrue(Utils.readFile(ma.getMethylationFile(Strand.WATSON, reference, sample)).indexOf
								("chr10\t6\tWATSON\tCG") != -1);
					assertTrue(Utils.readFile(ma.getMethylationFile(Strand.WATSON, reference, sample)).indexOf
								("chr10\t14\tWATSON\tCHG") != -1);
					System.err.println("===========================");

					System.err.println("====METHYLATION-CRICK=====");
					System.err.println(Utils.readFile(ma.getMethylationFile(Strand.CRICK, reference, sample)));
					System.err.println("===========================");


					System.err.println("=====METHYLCYTOSINES=======");
					System.err.println(Utils.readFile(ma.getMethylcytosinesFile(reference, sample)));
					System.err.println(Utils.readFile(ma.getMethylcytosinesVCFFile(reference, sample)));
					System.err.println("===========================");

					System.err.println("==========SUMMARY==========");
					System.err.println(Utils.readFile(ma.getSummaryFile(reference, sample)));
					System.err.println("===========================");

				}
			}
		} finally {
			Utils.deleteDirOnJVMExit(project.getProjectDirectory());
		}

	}

	@Test
	public void analysisTrimAndBad() throws IOException, InterruptedException {
		Project project = prepareProject();
		try {

			MethylationAnalysis ma = new MethylationAnalysis(project);

			List<File> bedFiles = Arrays.asList(new File(Utils.getBedsDirectory()).listFiles(new FilenameFilter() {

				@Override
				public boolean accept(File arg0, String arg1) {
					return arg0.getName().endsWith(".bed");
				}
			}));


			for (Sample sample : project.getSamples()) {
				for (Reference reference : project.getReferences()) {

					ma.analyzeWithErrorFromControlGenome(
							reference,
							sample,
							true, //trim
							4,
							true, //ambiguous
							false, //only with one aligment
							true, //bad
							false,
							true,
							1,
							0.01,
							1,
							bedFiles,
							"control");
					assertTrue(ma.getSummaryFile(reference, sample).exists());

					System.err.println("====METHYLATION-WATSON=====");
					System.err.println(Utils.readFile(ma.getMethylationFile(Strand.WATSON, reference, sample)));
					//trim, this line must have depth 4, not 5 due to trimming...

						assertTrue(Utils.readFile(ma.getMethylationFile(Strand.WATSON, reference, sample)).indexOf
								("chr10\t59\tWATSON\tCG\t4\t4\t4\t1.0\tCCCC") != -1);
					System.err.println("===========================");

					System.err.println("====METHYLATION-CRICK=====");
					System.err.println(Utils.readFile(ma.getMethylationFile(Strand.CRICK, reference, sample)));
					//trim test, the following line must be missing
					assertTrue(Utils.readFile(ma.getMethylationFile(Strand.CRICK, reference, sample)).indexOf
							("chr10\t16\tCRICK\tCHG\t1\t1\t1\t1.0\tG\t0.0\tfalse\tfalse") == -1);
					System.err.println("===========================");


					System.err.println("=====METHYLCYTOSINES=======");
					System.err.println(Utils.readFile(ma.getMethylcytosinesFile(reference, sample)));
					System.err.println("===========================");

					System.err.println("==========SUMMARY==========");
					System.err.println(Utils.readFile(ma.getSummaryFile(reference, sample)));
					System.err.println("===========================");

				}
			}
		} finally {
			Utils.deleteDirOnJVMExit(project.getProjectDirectory());
		}
	}


	@Test
	public void analysisSingleThread() throws IOException, InterruptedException {
		Project project = prepareProject();
		try {

			MethylationAnalysis ma = new MethylationAnalysis(project);

			List<File> bedFiles = Arrays.asList(new File(Utils.getBedsDirectory()).listFiles(new FilenameFilter() {

				@Override
				public boolean accept(File dir, String name) {
					return name.endsWith(".bed");
				}
			}));

			for (Sample sample : project.getSamples()) {
				for (Reference reference : project.getReferences()) {

					ma.analyzeWithErrorFromControlGenome(
							reference,
							sample,
							false, //notrim
							4,
							true,
							false,
							true,
							false,
							true,
							1,
							0.01,
							1,
							bedFiles,
							"control");
					assertTrue(ma.getSummaryFile(reference, sample).exists());

					System.err.println("====METHYLATION-WATSON=====");
					System.err.println(Utils.readFile(ma.getMethylationFile(Strand.WATSON, reference, sample)));
						assertTrue(Utils.readFile(ma.getMethylationFile(Strand.WATSON, reference, sample)).indexOf
								("chr10\t6\tWATSON\tCG") != -1);
						assertTrue(Utils.readFile(ma.getMethylationFile(Strand.WATSON, reference, sample)).indexOf
								("chr10\t14\tWATSON\tCHG") != -1);
					System.err.println("===========================");

					System.err.println("====METHYLATION-CRICK=====");
					System.err.println(Utils.readFile(ma.getMethylationFile(Strand.CRICK, reference, sample)));
					if (!this.bowtie2Local) {
						assertTrue(Utils.readFile(ma.getMethylationFile(Strand.CRICK, reference, sample)).indexOf
								("chr10\t36\tCRICK\tCHH\t1\t1\t0\t0.0\tA\t1.0\tfalse\tfalse") != -1);
					}
					System.err.println("===========================");


					System.err.println("=====METHYLCYTOSINES=======");
					System.err.println(Utils.readFile(ma.getMethylcytosinesFile(reference, sample)));
					System.err.println("===========================");

					System.err.println("==========SUMMARY==========");
					System.err.println(Utils.readFile(ma.getSummaryFile(reference, sample)));
					System.err.println("===========================");

					System.err.println("==========REGIONS==========");
					for (File bed : bedFiles) {
						System.err.println(Utils.readFile(ma.getMethylatedRegionsFile(reference, sample, bed)));
						if (bed.getName().contains("track2")) {
							assertTrue(Utils.readFile(ma.getMethylatedRegionsFile(reference, sample, bed)).indexOf
										("From26to70\t5") != -1);
							assertTrue(Utils.readFile(ma.getMethylatedRegionsFile(reference, sample, bed)).indexOf
										("0.1296") != -1);
							assertTrue(Utils.readFile(ma.getMethylatedRegionsFile(reference, sample, bed)).indexOf
										("From51to80\t5") != -1);
						}
					}
					System.err.println("===========================");


				}
			}
		} finally {
			Utils.deleteDirOnJVMExit(project.getProjectDirectory());
		}
	}

	@Test
	public void testDepth() throws IOException, InterruptedException {
		Project project = prepareProject();
		try {

			MethylationAnalysis ma = new MethylationAnalysis(project);

			List<File> bedFiles = Arrays.asList(new File(Utils.getBedsDirectory()).listFiles(new FilenameFilter() {

				@Override
				public boolean accept(File arg0, String arg1) {
					return arg0.getName().endsWith(".bed");
				}
			}));

			int mindepth = 2;
			for (Sample sample : project.getSamples()) {
				for (Reference reference : project.getReferences()) {

					ma.analyzeWithErrorFromControlGenome(
							reference,
							sample,
							true,
							4,
							true,
							true,
							true,
							false,
							true,
							mindepth,
							0.01,
							4,
							bedFiles,
							"control");
					assertTrue(ma.getSummaryFile(reference, sample).exists());

					System.err.println("====METHYLATION-WATSON=====");
					System.err.println(Utils.readFile(ma.getMethylationFile(Strand.WATSON, reference, sample)));

					System.err.println("===========================");

					System.err.println("====METHYLATION-CRICK=====");
					System.err.println(Utils.readFile(ma.getMethylationFile(Strand.CRICK, reference, sample)));
					System.err.println("===========================");


					System.err.println("=====METHYLCYTOSINES=======");
					System.err.println(Utils.readFile(ma.getMethylcytosinesFile(reference, sample)));
					String[] lines = Utils.readFile(ma.getMethylcytosinesFile(reference, sample)).split("\n");
					for (String line : lines) {
						if (!line.startsWith("#") && line.length() > 0)
							assertTrue(Integer.parseInt(line.split("\t")[5]) >= mindepth);
					}

					System.err.println("===========================");

					System.err.println("==========SUMMARY==========");
					System.err.println(Utils.readFile(ma.getSummaryFile(reference, sample)));
					System.err.println("===========================");

				}
			}
		} finally {
			Utils.deleteDirOnJVMExit(project.getProjectDirectory());
		}

	}

	@Test
	public void testRemoveClonal() throws IOException, InterruptedException {
		Project project = prepareProject();
		try {

			MethylationAnalysis ma = new MethylationAnalysis(project);

			List<File> bedFiles = Arrays.asList(new File(Utils.getBedsDirectory()).listFiles(new FilenameFilter() {

				@Override
				public boolean accept(File arg0, String arg1) {
					return arg0.getName().endsWith(".bed");
				}
			}));

			int mindepth = 2;
			for (Sample sample : project.getSamples()) {
				for (Reference reference : project.getReferences()) {

					ma.analyzeWithErrorFromControlGenome(
							reference,
							sample,
							false,
							4,
							false,
							true,
							false,
							true, //remove clonal
							true,
							mindepth,
							0.01,
							4,
							bedFiles,
							"control");
					assertTrue(ma.getSummaryFile(reference, sample).exists());

					System.err.println("====METHYLATION-WATSON=====");
					System.err.println(Utils.readFile(ma.getMethylationFile(Strand.WATSON, reference, sample)));
					if (!this.bowtie2Local) {
						assertTrue(Utils.readFile(ma.getMethylationFile(Strand.WATSON, reference, sample)).indexOf
								("chr10\t59\tWATSON\tCG\t5\t5\t5\t1.0\tCCCCC") != -1);
					} else {
						assertTrue(Utils.readFile(ma.getMethylationFile(Strand.WATSON, reference, sample)).indexOf
								("chr10\t59\tWATSON\tCG\t2\t2\t2\t1.0\tCC") != -1);
					}


					System.err.println("===========================");

					System.err.println("====METHYLATION-CRICK=====");
					System.err.println(Utils.readFile(ma.getMethylationFile(Strand.CRICK, reference, sample)));
					System.err.println("===========================");


					System.err.println("=====METHYLCYTOSINES=======");
					System.err.println(Utils.readFile(ma.getMethylcytosinesFile(reference, sample)));
					String[] lines = Utils.readFile(ma.getMethylcytosinesFile(reference, sample)).split("\n");
					for (String line : lines) {
						if (!line.startsWith("#") && line.length() > 0)
							assertTrue(Integer.parseInt(line.split("\t")[5]) >= mindepth);
					}

					System.err.println("===========================");

					System.err.println("==========SUMMARY==========");
					System.err.println(Utils.readFile(ma.getSummaryFile(reference, sample)));
					System.err.println("===========================");

				}
			}
		} finally {
			Utils.deleteDirOnJVMExit(project.getProjectDirectory());
		}



	}

	@Test
	public void testWithUniqueAlignments() throws IOException, InterruptedException {
		Project project = prepareProject();
		try {

			MethylationAnalysis ma = new MethylationAnalysis(project);

			List<File> bedFiles = Arrays.asList(new File(Utils.getBedsDirectory()).listFiles(new FilenameFilter() {

				@Override
				public boolean accept(File arg0, String arg1) {
					return arg0.getName().endsWith(".bed");
				}
			}));

			int mindepth = 2;
			for (Sample sample : project.getSamples()) {
				for (Reference reference : project.getReferences()) {

					ma.analyzeWithErrorFromControlGenome(
							reference,
							sample,
							false,
							4,
							false,
							true, //with unique alignments
							false,
							true,
							true,
							mindepth,
							0.01,
							4,
							bedFiles,
							"control");
					assertTrue(ma.getSummaryFile(reference, sample).exists());

					System.err.println("====METHYLATION-WATSON=====");
					System.err.println(Utils.readFile(ma.getMethylationFile(Strand.WATSON, reference, sample)));

					// with the uniquealignments option, the --local mode generates in this dataset many reads with
					// multiple alignments, with are then filtered
					if (!this.bowtie2Local) {
						assertTrue(Utils.readFile(ma.getMethylationFile(Strand.WATSON, reference, sample)).indexOf
								("chr10\t59\tWATSON\tCG\t5\t5\t5\t1.0\tCCCCC") != -1);
					} else {
						assertTrue(Utils.readFile(ma.getMethylationFile(Strand.WATSON, reference, sample)).indexOf
								("chr10\t59\tWATSON\tCG\t2\t2\t2\t1.0\tCC") != -1);
					}
					System.err.println("===========================");

					System.err.println("====METHYLATION-CRICK=====");
					System.err.println(Utils.readFile(ma.getMethylationFile(Strand.CRICK, reference, sample)));
					System.err.println("===========================");


					System.err.println("=====METHYLCYTOSINES=======");
					System.err.println(Utils.readFile(ma.getMethylcytosinesFile(reference, sample)));
					System.err.println("===========================");

					System.err.println("==========SUMMARY==========");
					System.err.println(Utils.readFile(ma.getSummaryFile(reference, sample)));
					System.err.println("===========================");

				}
			}
		} finally {
			Utils.deleteDirOnJVMExit(project.getProjectDirectory());
		}
	}

}
