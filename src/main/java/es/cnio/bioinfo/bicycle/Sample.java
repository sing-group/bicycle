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

import java.io.File;
import java.io.FilenameFilter;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;

public class Sample {
	
	private String name="unnamed-sample";
	private List<File> readsFiles = new LinkedList<File>();
	private boolean paired = false;
	private String mate1Regexp="_1.fastq";
	//for paired-end
	private List<File> readsFilesM1 = new LinkedList<File>();
	private List<File> readsFilesM2 = new LinkedList<File>();
	
	private Project project;
	
	Sample(Project project, String name, boolean paired, String mate1Regexp) {
		super();
		this.name = name;
		this.project = project;
		this.paired = paired;
		this.mate1Regexp = mate1Regexp;
	}
	Sample(Project project, String name) {
		super();
		this.name = name;
		this.project = project;
	}
	
	public String getName() {
		return name;
	}
	
	public boolean isPaired() {
		return paired;
	}
	public List<File> getReadsFiles(){
		return Collections.unmodifiableList(readsFiles);
	}
		
	public List<File> getReadsMate1Files(){
		if (!this.project.isPaired()){
			throw new IllegalArgumentException("The project is not a paired-end project");
		}
		return Collections.unmodifiableList(readsFilesM1);
	}
	public List<File> getReadsMate2Files(){
		if (!this.project.isPaired()){
			throw new IllegalArgumentException("The project is not a paired-end project");
		}
		return Collections.unmodifiableList(readsFilesM2);
		
	}
	public Project getProject() {
		return project;
	}
	
	private void addReadsFile(File f){
		this.readsFiles.add(f);
		if (this.paired){
			if (f.getName().matches(".*"+this.mate1Regexp+".*")){
				this.readsFilesM1.add(f);
			}else{
				this.readsFilesM2.add(f);
			}
		}
	}
	
	public static List<Sample> buildSamples(Project project){
		return buildSamples(project, false, "");
	}
	
	/**
	 * Builds samples from a directory which should contain reads files. Each reads file in the root of
	 * the path is a independent sample (with 1 reads file). If there are subdirectories, each subdirectory is a
	 * sample with all reads files inside it belonging to this sample.
	 * 
	 * @param readsPath
	 * @return
	 */
	public static List<Sample> buildSamples(Project project, boolean paired, String mate1Regexp){
		File readsPath = project.getReadsDirectory();
		if (!readsPath.exists() || !readsPath.isDirectory()){
			throw new IllegalArgumentException("Cannot find "+readsPath.getAbsolutePath()+" or it is not a directory, or it is not accessible");
		}
		List<Sample> toret = new LinkedList<Sample>();
		for (File file : readsPath.listFiles(new FilenameFilter() {
				@Override
				public boolean accept(File dir, String name) {
					return !name.equals("CVS");
				}
				})){
			
			if (file.isDirectory()){
				if (file.isHidden()) continue;
				Sample s = new Sample(project, file.getName(), paired, mate1Regexp);
				
				File[] files = file.listFiles();
				Arrays.sort(files, new Comparator<File>(){
					@Override
					public int compare(File arg0, File arg1) {
						return arg0.getName().compareTo(arg1.getName());
					}
					
				});
				for (File subdirFile : files){
					if (!subdirFile.isDirectory() && !file.isHidden()){
						s.addReadsFile(subdirFile);
					}
				}
				
				s.validate();
				toret.add(s);
			}else if (!file.isHidden()){
				if (paired) throw new IllegalArgumentException("Paired-end samples files must be inside a directory per sample");
				Sample s = new Sample(project, file.getName(), paired, mate1Regexp);
				s.addReadsFile(file);
				toret.add(s);
				s.validate();
			}
		}
		return toret;
	}
	
	private void validate() {
		if (this.paired && this.readsFilesM1.size() != this.readsFilesM2.size()){
			throw new IllegalArgumentException("This sample is paired end, but has a different number of read files for each mate. mate 1: "+this.readsFilesM1+". mate 2: "+this.readsFilesM2);
		}
		
	}
	@Override
	public boolean equals(Object obj) {
		if (!( obj instanceof Sample)){
			return false;
		}
		Sample sample = (Sample) obj;
		return sample.getName().equals(this.getName()) && sample.readsFiles.equals(this.readsFiles);
	}
	
}
