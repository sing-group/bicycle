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

package es.cnio.bioinfo.bicycle;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

public class Project {

	private static final Logger logger = Logger.getLogger(Project.class.getSimpleName());

	public static final String WORKING_DIRECTORY = "workingDirectory" + File.separator;
	public static final String OUTPUT_DIRECTORY = "output" + File.separator;
	public static final String CONFIG_FILE = "config.txt";

	private File bowtieDirectory, bowtie2Directory, samtoolsDirectory, projectDirectory, referenceDirectory,
			readsDirectory;
	private boolean paired = false;
	private boolean directional = true;
	private String mate1regexp = null;

	private List<Sample> samples = new LinkedList<Sample>();

	private Map<String, String> properties = new HashMap<String, String>();

	protected Project() {
	}

	public File getBowtieDirectory() {
		return bowtieDirectory;
	}

	public File getBowtie2Directory() {
		return bowtie2Directory;
	}

	public File getSamtoolsDirectory() {
		return samtoolsDirectory;
	}

	public File getProjectDirectory() {
		return projectDirectory;
	}

	public File getReferenceDirectory() {
		return referenceDirectory;
	}

	public File getReadsDirectory() {
		return readsDirectory;
	}

	public File getWorkingDirectory() {
		return new File(projectDirectory.getAbsolutePath() + File.separator + WORKING_DIRECTORY);
	}

	public File getOutputDirectory() {
		return new File(projectDirectory.getAbsolutePath() + File.separator + OUTPUT_DIRECTORY);
	}

	public File getConfigFile() {
		return new File(projectDirectory.getAbsolutePath() + File.separator + CONFIG_FILE);
	}

	public List<Sample> getSamples() {
		return Collections.unmodifiableList(samples);
	}

	public boolean isPaired() {
		return paired;
	}

	public boolean isDirectional() {
		return directional;
	}

	String getMate1regexp() {
		return mate1regexp;
	}

	public List<Reference> getReferences() {
		List<File> files = Arrays.asList(this.getReferenceDirectory().listFiles(new FileFilter() {
			@Override
			public boolean accept(File file) {
				return file.isFile() && file.getName().endsWith(".fa");
			}
		}));

		ArrayList<Reference> toret = new ArrayList<Reference>();
		for (File refFile : files) {
			toret.add(new Reference(this, refFile));
		}

		return toret;
	}

	public void saveProject() throws IOException {
		BufferedWriter wr = new BufferedWriter(new FileWriter(this.getConfigFile()));
		wr.write("project_directory:" + this.getProjectDirectory().getAbsolutePath());
		wr.newLine();

		if (bowtieDirectory != null) {
			wr.write("bowtie_directory:" + bowtieDirectory.getAbsolutePath());
			wr.newLine();
		}

		if (bowtie2Directory != null) {
			wr.write("bowtie2_directory:" + bowtie2Directory.getAbsolutePath());
			wr.newLine();
		}

		if (samtoolsDirectory != null) {
			wr.write("samtools_directory:" + samtoolsDirectory.getAbsolutePath());
			wr.newLine();
		}

		wr.write("reference_directory:" + referenceDirectory.getAbsolutePath());
		wr.newLine();
		wr.write("reads_directory:" + readsDirectory.getAbsolutePath());
		wr.newLine();
		wr.write("paired:" + paired);
		wr.newLine();
		wr.write("directional:" + directional);
		wr.newLine();
		if (paired) {
			wr.write("mate1regexp:" + mate1regexp);
			wr.newLine();
		}

		for (String prop : this.properties.keySet()) {
			wr.write(prop + ":" + this.properties.get(prop));
			wr.newLine();
		}
		wr.close();

	}

	public void addProperty(String name, String value) {
		this.properties.put(name, value);
	}

	public String getProperty(String name) {
		return this.properties.get(name);
	}

	public static Project buildNewProject(File projectDirectory,
										  File referenceDirectory, File readsDirectory, File bowtieDirectory,
										  File bowtie2Directory, File samtoolsDirectory, boolean directional) throws IOException {
		return buildNewProject(projectDirectory, referenceDirectory,
				readsDirectory, bowtieDirectory, bowtie2Directory, samtoolsDirectory, directional, false, null);
	}

	public static Project buildNewProject(File projectDirectory,
										  File referenceDirectory, File readsDirectory, File bowtieDirectory,
										  File bowtie2Directory,
										  File samtoolsDirectory, boolean isDirectional, boolean isPaired, String
												  mate1regexp)
			throws IOException {
		if (projectDirectory.exists()) {
			throw new IllegalArgumentException("The project folder already exists, please remove it first\n");
		}

		logger.info("Creating new project");
		Project p = new Project();
		p.projectDirectory = projectDirectory;

		p.getProjectDirectory().mkdirs();
		p.getWorkingDirectory().mkdirs();
		p.getOutputDirectory().mkdirs();


		p.bowtieDirectory = bowtieDirectory;
		p.bowtie2Directory = bowtie2Directory;

		p.samtoolsDirectory = samtoolsDirectory;
		p.paired = isPaired;
		p.directional = isDirectional;
		p.mate1regexp = mate1regexp;
		//check the reference directory
		logger.info("Looking for the reference(s) directory " + referenceDirectory.getAbsolutePath() + "...");
		if (referenceDirectory.exists()) {
			p.referenceDirectory = referenceDirectory;
		} else {
			String error = "unable to find " + referenceDirectory.getAbsolutePath();
			logger.severe(error);
			throw new IllegalArgumentException(error);
		}


		//check reads directory
		logger.info("Looking for the read(s) directory " + readsDirectory.getAbsolutePath() + "...");
		if (readsDirectory.exists()) {
			p.readsDirectory = readsDirectory;
			p.samples = Sample.buildSamples(p, p.paired, p.directional, p.mate1regexp);
			logger.info("Num of samples: " + p.samples.size());
		} else {
			String error = "unable to find " + readsDirectory.getAbsolutePath();
			logger.severe(error);
			throw new IllegalArgumentException(error);
		}
		p.saveProject();
		logger.info("Project creation OK");
		return p;
	}

	public static Project readFromDirectory(File projectDirectory) {
		if (projectDirectory.exists() && projectDirectory.isDirectory()) {
			Project p = new Project();
			String[] lines = Tools.readFile(projectDirectory + File.separator + "config.txt").split("[\n]");
			for (String line : lines) {
				if (line.startsWith("#")) continue;
				String[] tokens = line.split("[:]");
				if (tokens.length < 2) continue; //not key:value line, ignore

				String key = tokens[0];
				String value = tokens[1];
				if (key.equals("project_directory"))
					p.projectDirectory = new File(value);
				else if (key.equals("reference_directory"))
					p.referenceDirectory = new File(value);
				else if (key.equals("reads_directory"))
					p.readsDirectory = new File(value);
				else if (key.equals("bowtie_directory"))
					p.bowtieDirectory = new File(value);
				else if (key.equals("bowtie2_directory"))
					p.bowtie2Directory = new File(value);
				else if (key.equals("samtools_directory"))
					p.samtoolsDirectory = new File(value);
				else if (key.equals("paired")) {
					p.paired = Boolean.parseBoolean(value);
				} else if (key.equals("directional")) {
					p.directional = Boolean.parseBoolean(value);
				} else if (key.equals("mate1regexp"))
					p.mate1regexp = value;
				else
					p.addProperty(key, value);
			}
			p.samples.addAll(Sample.buildSamples(p, p.paired, p.directional, p.mate1regexp));
			return p;
		} else {
			throw new RuntimeException("Path " + projectDirectory + " was not found or it is not a directory or it cannot be accessed");
		}
	}
}
