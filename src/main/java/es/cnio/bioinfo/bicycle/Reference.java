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
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;

public class Reference {

	private File referenceFile;
	private Project project;

	Reference(Project project, File referenceFile) {
		this.project = project;
		this.referenceFile = referenceFile;
	}
	
	public File getReferenceFile() {
		return referenceFile;
	}
	
	public Project getProject() {
		return project;
	}
	
	@Override
	public boolean equals(Object obj) {
		return (obj instanceof Reference && ((Reference) obj).getReferenceFile().equals(this.referenceFile));
	}
	
	@Override
	public int hashCode() {
		return this.referenceFile.hashCode();
	}
	
	public List<String> getSequenceNames() {
		try {
			File index = new File(this.project.getWorkingDirectory()+File.separator+this.referenceFile+".index"); 
			if (!index.exists()) {
					createReferenceIndex(index);
			}
			
			Scanner sc = new Scanner(index);
			
			List<String> sequenceNames = new LinkedList<>();
			while (sc.hasNextLine()) {
				sequenceNames.add(sc.nextLine().split("\t")[0]);
			}
			sc.close();
			
			return sequenceNames;
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}
	}

	private void createReferenceIndex(File index) throws FileNotFoundException {
		PrintStream indexOut = new PrintStream(new FileOutputStream(index));
		Scanner sc = new Scanner(this.referenceFile);
		while (sc.hasNextLine()) {
			String line = sc.nextLine();
			if (line.startsWith(">")) {
				indexOut.println(line.substring(1));
			}
		}
		sc.close();
		indexOut.close();
			
	}
}
