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
import java.io.IOException;

import org.junit.Test;

import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Reference;
import es.cnio.bioinfo.bicycle.operations.ReferenceBisulfitation;
import es.cnio.bioinfo.bicycle.operations.ReferenceBisulfitation.Replacement;

public class ReferenceBisulfitationTest {

	@Test
	public void test() throws IOException {
		String aSequence = "ACCCCGGGTTT";
		String aBisulfitedCTSequence = "ATTTTGGGTTT";
		String aBisulfitedGASequence = "ACCCCAAATTT";

		File tempDir = Utils.generateTempDirName("newproject");
		File refsDir = Utils.generateTempDirName("refs");
		refsDir.mkdir();
		File genome1 = Utils.touchFile(refsDir, "genome01.fasta");
		Utils.append(genome1, ">seq\n" + aSequence + "\n");
		File genome2 = Utils.touchFile(refsDir, "genome02.fasta");
		Utils.append(genome2, ">seq\n" + aSequence + "\n");
		try {
			Project p = Project.buildNewProject(
					tempDir,
					refsDir,
					new File(Utils.getReadsDirectory()),
					new File(Utils.getBowtiePath()),
					new File(Utils.getBowtie2Path()),
					new File(Utils.getSamtoolsPath()),
					true);

			ReferenceBisulfitation rb = new ReferenceBisulfitation(p);


			for (Reference ref : p.getReferences()) {
				rb.computeReferenceBisulfitation(Replacement.CT, ref, false);
				rb.computeReferenceBisulfitation(Replacement.GA, ref, false);

				assertTrue(rb.getBisulfitedReference(Replacement.CT, ref).exists());
				assertTrue(rb.getBisulfitedReference(Replacement.GA, ref).exists());
				assertTrue(Utils.readFile(rb.getBisulfitedReference(Replacement.CT, ref)).indexOf
						(aBisulfitedCTSequence) != -1);
				assertTrue(Utils.readFile(rb.getBisulfitedReference(Replacement.GA, ref)).indexOf
						(aBisulfitedGASequence) != -1);
			}


		} finally {
			Utils.deleteDir(tempDir);
			Utils.deleteDir(refsDir);
		}
	}

}
