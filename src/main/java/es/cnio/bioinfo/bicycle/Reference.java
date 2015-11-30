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
}
