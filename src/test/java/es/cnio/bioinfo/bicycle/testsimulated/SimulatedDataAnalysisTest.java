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

package es.cnio.bioinfo.bicycle.testsimulated;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Reference;
import es.cnio.bioinfo.bicycle.Sample;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment.Bowtie1Quals;
import es.cnio.bioinfo.bicycle.operations.MethylationAnalysis;
import es.cnio.bioinfo.bicycle.operations.ReferenceBisulfitation;
import es.cnio.bioinfo.bicycle.operations.ReferenceBisulfitation.Replacement;
import es.cnio.bioinfo.bicycle.test.Utils;

@RunWith(Parameterized.class)
public class SimulatedDataAnalysisTest {

	private final int bowtieVersion;
	private final boolean bowtie2Local;

	@Parameters
	public static Collection<Object[]> data() {
		return Arrays.asList(new Object[][]{{1, false}, {2,false}, {2, true}});
	}

	public SimulatedDataAnalysisTest(int bowtieVersion, boolean bowtie2Local)
	{
		this.bowtieVersion = bowtieVersion;
		this.bowtie2Local = bowtie2Local;
	}

	private Project prepareProject() throws IOException {

		File tempDir = Utils.generateTempDirName("newproject-simulated-data");

		Project p = Project.buildNewProject(
				tempDir,
				new File(Utils.getSimulatedDataReferenceDirectory()),
				new File(Utils.getSimulatedDataReadsDirectory()),
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
		//	Project project = Project.readFromDirectory(new File
		// ("/tmp/newproject-simulated-data181c07eb-a533-4eec-91ab-47177f88f7ff"));
		try {

			MethylationAnalysis ma = new MethylationAnalysis(project);

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
							1,
							0.01,
							1,
							new ArrayList<File>(),
							"Ecoli",
							false);
					assertTrue(ma.getSummaryFile(reference, sample).exists());


					System.err.println("==========SUMMARY==========");
					System.err.println(Utils.readFile(ma.getSummaryFile(reference, sample)));
					System.err.println("===========================");

					// the 30% methylation level on CG
					Matcher matcher = Pattern.compile("CG: [0-9]+/[0-9]+ \\((0\\.[0-9]+)\\)").matcher(Utils.readFile(ma.getSummaryFile(reference, sample)));
					assertTrue(matcher.find());
					assertTrue(Math.abs(Double.parseDouble(matcher.group(1)) - 0.30d) < 0.01);
					System.err.println("===========================");
				}


			}
		} finally {
			Utils.deleteDirOnJVMExit(project.getProjectDirectory());
		}
	}
}
