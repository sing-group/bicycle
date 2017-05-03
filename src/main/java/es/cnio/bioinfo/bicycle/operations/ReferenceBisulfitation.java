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

package es.cnio.bioinfo.bicycle.operations;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Logger;

import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Reference;

public class ReferenceBisulfitation {
	private static final Logger logger = Logger.getLogger(ReferenceBisulfitation.class.getSimpleName());
	
	private static final String BISULFITED_DIR_PROPERTY="bisulfited_reference_dir";
	private Project project;

	public ReferenceBisulfitation(Project p) {	
		this.project = p;
	}
	public File getBisulfitedReference(Replacement replacement, Reference ref){
		if (!project.getReferences().contains(ref)){
			throw new IllegalArgumentException("This project doesn't has this reference: "+ref.getReferenceFile());
		}
		String bisulfitedReferenceDirName = this.project.getProperty(BISULFITED_DIR_PROPERTY+replacement.name());
		if (bisulfitedReferenceDirName != null){
			return new File(bisulfitedReferenceDirName+File.separator+ref.getReferenceFile().getName()+"_bisulfited_"+replacement.name());
		}else{
			File tryOnReference = new File(this.project.getReferenceDirectory()+File.separator+ref.getReferenceFile().getName()+"_bisulfited_"+replacement.name());
			File tryOnWorking = new File(this.project.getWorkingDirectory()+File.separator+ref.getReferenceFile().getName()+"_bisulfited_"+replacement.name());
			
			if (tryOnReference.exists()){
				return tryOnReference;
			}else if (tryOnWorking.exists()){
				return tryOnWorking;
			}
			throw new IllegalArgumentException("property not found in project "+BISULFITED_DIR_PROPERTY+replacement.name()+". Did you perform reference bisulfitation?");
			
		}
		
	}
	public void computeReferenceBisulfitation(Replacement replacement, Reference reference, boolean onWorkingDir) throws IOException{
		
		
		File inputFile = reference.getReferenceFile();
		logger.info("Starting "+replacement.name()+" in-silico bisulfitation for reference file: "+inputFile
				.toString().replaceAll
				(project.getReferenceDirectory()+File.separator, ""));
		this.project.addProperty(BISULFITED_DIR_PROPERTY+replacement.name(), onWorkingDir?project.getWorkingDirectory().getAbsolutePath():project.getReferenceDirectory().getAbsolutePath());
		this.project.saveProject();
		
		
		BufferedReader br = new BufferedReader(new FileReader(inputFile));
		
		File outputFile=getBisulfitedReference(replacement, reference);
		if (outputFile.exists()) outputFile.delete();
		BufferedWriter wr = new BufferedWriter(new FileWriter(outputFile));			
		
		
		for (String line = br.readLine(); line!=null; line = br.readLine()){				
	        if(!line.startsWith(">")){
        	 	
        	 	final String modifiedLine=replacement.replace(line);
        	 	wr.write(modifiedLine);
        	 	wr.newLine();
			
	        }else{
	         	String newFastaHeader=new StringBuilder(line.replace(' ','_')).toString();
	         	wr.write(newFastaHeader);
	         	wr.newLine();
	        }
	    } // end while
		
		br.close();			
		wr.close();			
		logger.info("In-silico bisulfitation OK");
		
	}
	
	public enum Replacement{
		CT('c', 't'), GA('g','a');
		
		private char from;
		private char to;
		private char FROM, TO;
		private Replacement(char from, char to){
			this.from=from;
			this.to = to;
			this.FROM = Character.toUpperCase(from);
			this.TO = Character.toUpperCase(to);
		}
		public String replace(String line){
			return line.replace(from, to).replace(FROM, TO);
		}
	}
}
