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
import java.util.Arrays;
import java.util.Collection;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Reference;
import es.cnio.bioinfo.bicycle.Sample;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment.Bowtie1Quals;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment.Bowtie2Quals;
import es.cnio.bioinfo.bicycle.operations.ReferenceBisulfitation;
import es.cnio.bioinfo.bicycle.operations.ReferenceBisulfitation.Replacement;

@RunWith(Parameterized.class)
public class AlignmentTest {

	private final int bowtieVersion;
	private final boolean bowtie2Local;

	@Parameters
	public static Collection<Object[]> data() {
		return Arrays.asList(new Object[][]{{1, false}, {2,false}, {2, true}});
	}
	
	public AlignmentTest(int bowtieVersion, boolean bowtie2Local) {
		this.bowtieVersion = bowtieVersion;
		this.bowtie2Local = bowtie2Local;
	}
	
	@Test(expected = IllegalArgumentException.class)
	public void buildIndexWithNoBisulfitedRef() throws IOException {
		File tempDir = Utils.generateTempDirName("newproject");
		try {
			Project p = Project.buildNewProject(
					tempDir,
					new File(Utils.getReferenceDirectory()),
					new File(Utils.getReadsDirectory()),
					new File(Utils.getBowtiePath()),
					new File(Utils.getBowtie2Path()),
					new File(Utils.getSamtoolsPath()),
					true);

			boolean bisulfitedReferenceExists = false;
			ReferenceBisulfitation rb = new ReferenceBisulfitation(p);
			for (Reference ref : p.getReferences()) {
				if (rb.getBisulfitedReference(Replacement.CT, ref).exists()) {
					rb.getBisulfitedReference(Replacement.CT, ref).delete();
					bisulfitedReferenceExists = true;
				}
				if (rb.getBisulfitedReference(Replacement.GA, ref).exists()) {
					rb.getBisulfitedReference(Replacement.GA, ref).delete();
					bisulfitedReferenceExists = true;
				}
			}

			if (bisulfitedReferenceExists) {
				Utils.deleteDir(tempDir);
				tempDir = Utils.generateTempDirName("newproject");
				p = Project.buildNewProject(
						tempDir,
						new File(Utils.getReferenceDirectory()),
						new File(Utils.getReadsDirectory()),
						new File(Utils.getBowtiePath()),
						new File(Utils.getBowtie2Path()),
						new File(Utils.getSamtoolsPath()),
						true);
			}
			BowtieAlignment ba = new BowtieAlignment(p);

			for (Reference ref : p.getReferences()) {
				if (this.bowtieVersion == 1) {
					ba.buildBowtieIndex(ref);
				} else {
					ba.buildBowtie2Index(ref);
				}
			}
		} finally {
			Utils.deleteDir(tempDir);
		}
	}

	@Test()
	public void buildIndex() throws IOException {
		File tempDir = Utils.generateTempDirName("newproject");
		try {
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
				rb.computeReferenceBisulfitation(Replacement.CT, ref, false);
				rb.computeReferenceBisulfitation(Replacement.GA, ref, false);

				if (this.bowtieVersion == 1) {
					ba.buildBowtieIndex(ref);
				} else {
					ba.buildBowtie2Index(ref);
				}
			}
		} finally {
			Utils.deleteDir(tempDir);
		}
	}

	@Test()
	public void align() throws IOException {
		File tempDir = Utils.generateTempDirName("newproject");
		try {
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
				rb.computeReferenceBisulfitation(Replacement.CT, ref, false);
				rb.computeReferenceBisulfitation(Replacement.GA, ref, false);
				if (this.bowtieVersion == 1) {
					ba.buildBowtieIndex(ref);
				} else {
					ba.buildBowtie2Index(ref);
				}
			}
			for (Sample sample : p.getSamples()) {

				for (Reference reference : p.getReferences()) {
					if (bowtieVersion == 1) {
						ba.performBowtie1Alignment(sample, reference, false, 4, 140, 20, 0, 64, Bowtie1Quals.BEFORE_1_3);
					} else {
						ba.performBowtie2Alignment(sample, reference, false, 4, this.bowtie2Local, 15, 2, 20,
										!this.bowtie2Local?"S,1,1.15":"S,1,0.75",
										!this.bowtie2Local?"L,-0.6,-0.6":"G,20,8",
										0, Bowtie2Quals.BEFORE_1_3);
					}
				}
			}
		} finally {
			Utils.deleteDir(tempDir);
		}
	}
}
