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

import java.io.File;
import java.io.IOException;

import org.broadinstitute.sting.utils.baq.BAQ.QualityMode;
import org.junit.Test;

import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Reference;
import es.cnio.bioinfo.bicycle.Sample;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment;
import es.cnio.bioinfo.bicycle.operations.ReferenceBisulfitation;
import es.cnio.bioinfo.bicycle.operations.SampleBisulfitation;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment.Quals;
import es.cnio.bioinfo.bicycle.operations.ReferenceBisulfitation.Replacement;

public class AlignmentTest {

	@Test(expected=IllegalArgumentException.class)
	public void buildIndexWithNoBisulfitedRef() throws IOException {
		File tempDir = Utils.generateTempDirName("newproject");		
		try{
			Project p = Project.buildNewProject(
					tempDir,
					new File(Utils.getReferenceDirectory()), 
					new File(Utils.getReadsDirectory()), 
					new File(Utils.getBowtiePath()),
					new File(Utils.getSamtoolsPath()));
			
			boolean bisulfitedReferenceExists = false;
			ReferenceBisulfitation rb = new ReferenceBisulfitation(p);
			for (Reference ref : p.getReferences()){	
				if (rb.getBisulfitedReference(Replacement.CT, ref).exists()){
					rb.getBisulfitedReference(Replacement.CT, ref).delete();
					bisulfitedReferenceExists = true;
				}
				if (rb.getBisulfitedReference(Replacement.GA, ref).exists()){
					rb.getBisulfitedReference(Replacement.GA, ref).delete();
					bisulfitedReferenceExists = true;
				}
			}
			
			if (bisulfitedReferenceExists){
				Utils.deleteDir(tempDir);
				tempDir = Utils.generateTempDirName("newproject");		
				p = Project.buildNewProject(
						tempDir,
						new File(Utils.getReferenceDirectory()), 
						new File(Utils.getReadsDirectory()), 
						new File(Utils.getBowtiePath()),
						new File(Utils.getSamtoolsPath()));
			}
			BowtieAlignment ba = new BowtieAlignment(p);
			
			for (Reference ref : p.getReferences()){
				
				ba.buildBowtieIndex(ref);
			}
		}finally{
			Utils.deleteDir(tempDir);
		}		
	}
	
	@Test()
	public void buildIndex() throws IOException {
		File tempDir = Utils.generateTempDirName("newproject");		
		try{
			Project p = Project.buildNewProject(
					tempDir,
					new File(Utils.getReferenceDirectory()), 
					new File(Utils.getReadsDirectory()), 
					new File(Utils.getBowtiePath()),
					new File(Utils.getSamtoolsPath()));
			
			ReferenceBisulfitation rb = new ReferenceBisulfitation(p);
			BowtieAlignment ba = new BowtieAlignment(p);
			
			for (Reference ref : p.getReferences()){				
				rb.computeReferenceBisulfitation(Replacement.CT, ref, false);
				rb.computeReferenceBisulfitation(Replacement.GA, ref, false);
				ba.buildBowtieIndex(ref);
			}
		}finally{
			Utils.deleteDir(tempDir);
		}		
	}
	
	@Test()
	public void align() throws IOException {
		File tempDir = Utils.generateTempDirName("newproject");		
		try{
			Project p = Project.buildNewProject(
					tempDir,
					new File(Utils.getReferenceDirectory()), 
					new File(Utils.getReadsDirectory()), 
					new File(Utils.getBowtiePath()),
					new File(Utils.getSamtoolsPath()));
			
			ReferenceBisulfitation rb = new ReferenceBisulfitation(p);
			BowtieAlignment ba = new BowtieAlignment(p);
			
			for (Reference ref : p.getReferences()){				
				rb.computeReferenceBisulfitation(Replacement.CT, ref, false);
				rb.computeReferenceBisulfitation(Replacement.GA, ref, false);
				ba.buildBowtieIndex(ref);
			}
			for (Sample sample : p.getSamples()){
				SampleBisulfitation sb = new SampleBisulfitation(sample);
				sb.computeSampleBisulfitation(true);
				for (Reference reference : p.getReferences()){
					ba.performBowtieAlignment(sample, reference, 4, 140, 20, 0, 64, Quals.BEFORE_1_3);
				}
			}
			
		}finally{
			Utils.deleteDir(tempDir);
		}		
	}
	

}
