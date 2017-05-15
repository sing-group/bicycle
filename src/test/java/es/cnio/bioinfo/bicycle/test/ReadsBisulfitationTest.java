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
import es.cnio.bioinfo.bicycle.Sample;
import es.cnio.bioinfo.bicycle.operations.SampleBisulfitation;

public class ReadsBisulfitationTest {

	@Test
	public void test() throws IOException {
		File tempDir = Utils.generateTempDirName("newproject");
		File readsDir = Utils.generateTempDirName("reads");
		readsDir.mkdir();
		File sample1 = new File(readsDir.getAbsolutePath() + File.separator + "sample01");
		sample1.mkdir();
		File sample2 = new File(readsDir.getAbsolutePath() + File.separator + "sample02");
		sample2.mkdir();
		Utils.touchFile(sample1, "reads01.fastq");
		Utils.touchFile(sample1, "reads02.fastq");
		Utils.touchFile(sample1, "reads03.fastq");
		Utils.touchFile(sample1, "reads04.fastq");
		try {


			Project p = Project.buildNewProject(
					tempDir,
					new File(Utils.getReferenceDirectory()),
					readsDir,
					new File(Utils.getBowtiePath()),
					new File(Utils.getSamtoolsPath()),
					true);
			for (Sample sample : p.getSamples()) {
				SampleBisulfitation sb = new SampleBisulfitation(sample);
				sb.computeSampleBisulfitation(false);
				for (File readsFile : sample.getReadsFiles()) {
					assertTrue(sb.getBisulfitedFile(readsFile).exists());
				}

			}
		} finally {
			Utils.deleteDir(tempDir);
			Utils.deleteDir(readsDir);
		}
	}

}
