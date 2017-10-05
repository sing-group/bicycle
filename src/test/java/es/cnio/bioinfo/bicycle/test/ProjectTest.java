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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.junit.Test;

import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Sample;
import es.cnio.bioinfo.bicycle.cli.BicycleApplication;
import es.cnio.bioinfo.bicycle.cli.CLIApplication;

public class ProjectTest {

	@Test(expected = IllegalArgumentException.class)
	public void testCreateProjectOnExistentPath() throws IOException {
		File tempDir = Utils.generateTempDirName("newproject");
		try {
			tempDir.mkdir(); //ouch!

			Project.buildNewProject(
					tempDir,
					new File(Utils.getReferenceDirectory()),
					new File(Utils.getReadsDirectory()),
					new File(Utils.getBowtiePath()),
					new File(Utils.getBowtie2Path()),
					new File(Utils.getSamtoolsPath()),
					true);
		} finally {
			Utils.deleteDir(tempDir);
		}
	}

	@Test
	public void testCreateProject() throws IOException {
		File tempDir = Utils.generateTempDirName("newproject");
		try {
			Project.buildNewProject(
					tempDir,
					new File(Utils.getReferenceDirectory()),
					new File(Utils.getReadsDirectory()),
					new File(Utils.getBowtiePath()),
					new File(Utils.getBowtie2Path()),
					new File(Utils.getSamtoolsPath()),
					true);

			assertTrue(new File(tempDir.getAbsolutePath() + File.separator + Project.OUTPUT_DIRECTORY).exists());
			assertTrue(new File(tempDir.getAbsolutePath() + File.separator + Project.WORKING_DIRECTORY).exists());
			assertTrue(new File(tempDir.getAbsolutePath() + File.separator + Project.CONFIG_FILE).exists());
		} finally {
			Utils.deleteDir(tempDir);
		}
	}

	@Test
	public void testLoadProject() throws IOException {
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

			Project p2 = Project.readFromDirectory(tempDir);

			assertEquals(p.getReadsDirectory(), p2.getReadsDirectory());
			assertEquals(p.getBowtieDirectory(), p2.getBowtieDirectory());
			assertEquals(p.getSamtoolsDirectory(), p2.getSamtoolsDirectory());
			assertEquals(p.getSamples(), p2.getSamples());
		} finally {
			Utils.deleteDir(tempDir);
		}
	}


	@Test
	public void testSamples() throws IOException {
		File tempDir = Utils.generateTempDirName("newproject");
		File samplesDir = Utils.generateTempDirName("samples");
		samplesDir.mkdir();
		try {
			//prepare 3 samples (one in the root directory, and two in subdirectories)
			File sample0 = new File(samplesDir + File.separator + "reads00.fastq");
			List<File> sample0files = new LinkedList<File>();
			sample0files.add(sample0);
			Utils.touchFile(sample0);


			File sample1 = new File(samplesDir + File.separator + "SAMPLE_1");
			sample1.mkdir();
			List<File> sample1files = new LinkedList<File>();
			File fastq = new File(sample1 + File.separator + "reads01.fastq");
			Utils.touchFile(fastq);
			sample1files.add(fastq);
			fastq = new File(sample1 + File.separator + "reads02.fastq");
			Utils.touchFile(fastq);
			sample1files.add(fastq);


			List<File> sample2files = new LinkedList<File>();
			File sample2 = new File(samplesDir + File.separator + "SAMPLE_2");
			sample2.mkdir();
			fastq = new File(sample2 + File.separator + "reads03.fastq");
			Utils.touchFile(fastq);
			sample2files.add(fastq);
			fastq = new File(sample2 + File.separator + "reads04.fastq");
			Utils.touchFile(fastq);
			sample2files.add(fastq);

			Project p = Project.buildNewProject(
					tempDir,
					new File(Utils.getReferenceDirectory()),
					samplesDir,
					new File(Utils.getBowtiePath()),
					new File(Utils.getBowtie2Path()),
					new File(Utils.getSamtoolsPath()),
					true);

			List<Sample> samples = p.getSamples();
			assertEquals(3, samples.size());

			List<String> sampleNames = new LinkedList<String>();
			HashMap<String, List<File>> sampleFileNames = new HashMap<String, List<File>>();
			for (Sample s : p.getSamples()) {
				sampleNames.add(s.getName());
				List<File> fileNames = sampleFileNames.get(s.getName());
				if (fileNames == null) {
					fileNames = new LinkedList<File>();
					sampleFileNames.put(s.getName(), fileNames);
				}
				fileNames.addAll(s.getReadsFiles());

			}
			assertTrue(sampleNames.contains("SAMPLE_1"));
			assertTrue(sampleNames.contains("SAMPLE_2"));
			assertTrue(sampleNames.contains("reads00.fastq"));

			assertTrue(sampleFileNames.get("reads00.fastq").containsAll(sample0files));
			assertTrue(sampleFileNames.get("SAMPLE_1").containsAll(sample1files));
			assertTrue(sampleFileNames.get("SAMPLE_2").containsAll(sample2files));

		} finally {
			Utils.deleteDir(tempDir);
			Utils.deleteDir(samplesDir);
		}
	}

	@Test
	public void testCLI() throws IOException {
		File tempDir = Utils.generateTempDirName("newproject");
		try {
			CLIApplication methylapp = new BicycleApplication();
			methylapp.run(new String[]{
					"create-project",
					"-p", tempDir.getAbsolutePath(),
					"-r", Utils.getReferenceDirectory(),
					"-f", Utils.getReadsDirectory(),
					"-b", Utils.getBowtiePath(),
					"-samtools-directory", Utils.getSamtoolsPath()
			});
			assertTrue(new File(tempDir.getAbsolutePath() + File.separator + Project.OUTPUT_DIRECTORY).exists());
			assertTrue(new File(tempDir.getAbsolutePath() + File.separator + Project.WORKING_DIRECTORY).exists());
			assertTrue(new File(tempDir.getAbsolutePath() + File.separator + Project.CONFIG_FILE).exists());
		} finally {
			Utils.deleteDir(tempDir);
		}
	}

}
