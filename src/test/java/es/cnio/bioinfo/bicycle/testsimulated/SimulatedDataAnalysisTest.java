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

import org.junit.Test;

import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Reference;
import es.cnio.bioinfo.bicycle.Sample;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment.Quals;
import es.cnio.bioinfo.bicycle.operations.MethylationAnalysis;
import es.cnio.bioinfo.bicycle.operations.ReferenceBisulfitation;
import es.cnio.bioinfo.bicycle.operations.ReferenceBisulfitation.Replacement;
import es.cnio.bioinfo.bicycle.operations.SampleBisulfitation;
import es.cnio.bioinfo.bicycle.test.Utils;

public class SimulatedDataAnalysisTest {

	
	
	private Project prepareProject() throws IOException {
		
		File tempDir = Utils.generateTempDirName("newproject-simulated-data");	
		
		Project p = Project.buildNewProject(
				tempDir,
				new File(Utils.getSimulatedDataReferenceDirectory()), 
				new File(Utils.getSimulatedDataReadsDirectory()), 
				new File(Utils.getBowtiePath()),
				new File(Utils.getSamtoolsPath()));
		
		ReferenceBisulfitation rb = new ReferenceBisulfitation(p);
		BowtieAlignment ba = new BowtieAlignment(p);
		
		for (Reference ref : p.getReferences()){				
			rb.computeReferenceBisulfitation(Replacement.CT, ref, true);
			rb.computeReferenceBisulfitation(Replacement.GA, ref, true);
			ba.buildBowtieIndex(ref);
		}
		for (Sample sample : p.getSamples()){
			SampleBisulfitation sb = new SampleBisulfitation(sample);
			sb.computeSampleBisulfitation(true);
			for (Reference reference : p.getReferences()){
				ba.performBowtieAlignment(sample, reference, 4, 140, 20, 0, 64, Quals.BEFORE_1_3);
			}
		}
			
		return p;		
	}
	
	@Test
	public void analysisMultithread() throws IOException, InterruptedException{
		Project project = prepareProject();
	//	Project project = Project.readFromDirectory(new File("/tmp/newproject-simulated-data181c07eb-a533-4eec-91ab-47177f88f7ff"));
		try{
			
			MethylationAnalysis ma = new MethylationAnalysis(project);

			
			for (Sample sample: project.getSamples()){
				for (Reference reference : project.getReferences()){
					
					ma.analyzeWithErrorFromControlGenome(
							reference, 
							sample, 
							true, 
							4, 
							true, 
							true,
							false,
							true,
							1,
							0.01, 
							4,
							new ArrayList<File>(), 
							"Ecoli");
					assertTrue(ma.getSummaryFile(reference, sample).exists());
					

					
					
					System.err.println("==========SUMMARY==========");
					System.err.println(Utils.readFile(ma.getSummaryFile(reference, sample)));
					assertTrue(Utils.readFile(ma.getSummaryFile(reference, sample)).indexOf("0.29")!=-1); //assert the 30% methylation level on CG
					
					System.err.println("===========================");
					
				}
			}
		}finally{
			//TestUtils.deleteDir(project.getProjectDirectory());
		}
		
	}
	
	
}
