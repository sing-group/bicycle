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
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Reads several sorted GPFiles, merging them and giving the lines sorted by position
 * @author lipido
 *
 */
public class GPFilesReader{

	private BufferedReader[] readers;
	private String[] current;

	private int lastLineReaderIndex = -1;
	
	private Map<String, Integer> sequenceIndexes = new HashMap<>();
	public GPFilesReader(List<String> sequenceNames, BufferedReader ... bufferedReaders) throws IOException {
		
		int i = 0;
		for (String sequenceName: sequenceNames) {
			sequenceIndexes.put(sequenceName, i++);
		}
		
		readers = bufferedReaders;
		current = new String[readers.length];
		
		i = 0;
		for (BufferedReader reader : readers){
			current[i++] = readLine(reader);
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
					}else if (this.sequenceIndexes.get(minimumContig)>this.sequenceIndexes.get(tokens[0])){
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
			this.current[whoIsMinimum] = readLine(this.readers[whoIsMinimum]);
			this.lastLineReaderIndex = whoIsMinimum;
		} else {
			this.lastLineReaderIndex = -1;
		}
		
		return minimum;
	}

	private String readLine(BufferedReader reader) throws IOException {
		String line = null;
		do {
			line = reader.readLine();
		} while (line != null && line.startsWith("#"));
		
		return line;
	}
	
	/**
	 * Returns the index of the buffer where the last returned line was extracted from
	 *  If no line was currently extracted of null was returned, this method returns -1.
	 *  
	 * @return The index of the buffer where the last line was extracted from.
	 */
	public int getLastLineReaderIndex() {
		return lastLineReaderIndex;
	}
}
