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

package es.cnio.bioinfo.bicycle.gatk;

import java.io.BufferedReader;
import java.io.IOException;

import net.sf.samtools.SAMSequenceDictionary;

/**
 * Reads several sorted GPFiles, merging them and giving the lines sorted by position
 * @author lipido
 *
 */
public class GPFilesReader{

	private BufferedReader[] readers;
	private String[] current;
	private SAMSequenceDictionary sequenceDictionary;
	public GPFilesReader(SAMSequenceDictionary masterSequenceDictionary, BufferedReader ... bufferedReaders) throws IOException {
	
		readers = bufferedReaders;
		current = new String[readers.length];
		this.sequenceDictionary = masterSequenceDictionary;
		int i = 0;
		for (BufferedReader reader : readers){
			current[i++] = reader.readLine();
		}
	}
	
	
	public String readLine() throws IOException{
		String minimum = null;
		String minimumContig = null;
		long minimumPos = Long.MAX_VALUE;
		
		int i = 0;
		int whoIsMinimum=0;
		for (String currentLine : current){
			if (currentLine != null){
				if (minimumContig == null){
					String[] tokens = currentLine.split("\t",3);					
					minimum = currentLine;
					minimumContig = tokens[0];
					minimumPos = Long.parseLong(tokens[1]);
					whoIsMinimum = i;
				}else{
					String[] tokens = currentLine.split("\t",3);
					if (minimumContig.equals(tokens[0])){
						long currentPos = Long.parseLong(tokens[1]);
						if (currentPos < minimumPos){
							minimum = currentLine;
							minimumContig = tokens[0];
							minimumPos = currentPos;
							whoIsMinimum = i;
						}
					}else if (this.sequenceDictionary.getSequenceIndex(minimumContig)>this.sequenceDictionary.getSequenceIndex(tokens[0])){
						minimumContig = tokens[0];
						minimum = currentLine;
						minimumPos = Long.parseLong(tokens[1]);
						whoIsMinimum=i;
					}
				}
				
			}
			i++;
		}
		if (minimum != null){
			this.current[whoIsMinimum] = this.readers[whoIsMinimum].readLine();
		}
		return minimum;
	}
}
