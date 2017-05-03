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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class FastqSplitter {


	public static List<BufferedReader> splitfastq(List<File> fastqs, int maxChunks) throws IOException {
		
		//get total length
		long totalLength = 0;
		
		// which is the global starting byte (inclusive) of each file
		List<Long> fileStarts = new LinkedList<Long>();
		
		for (File f: fastqs){
			fileStarts.add(totalLength);
			totalLength += f.length();
		}
	
	
		
		long chunkStep = totalLength / maxChunks;
		//System.out.println("chunk step: "+chunkStep);
		
		long currentStart = 0;
		List<Long> chunks= new LinkedList<Long>();
		
		for (int chunk = 0; chunk<maxChunks; chunk++){
			long currentEnd = currentStart + chunkStep;
			
			final int fileIndex = getFileForPosition(currentEnd, fileStarts);
			File f = fastqs.get(fileIndex);
		
		
			final long previous = currentEnd - fileStarts.get(fileIndex);
			currentEnd = adjustPos(f, previous) + fileStarts.get(fileIndex);
			//System.out.println("adjusted "+f+" from "+previous+" to "+currentEnd);


			chunks.add(currentEnd);			
			currentStart = currentEnd + 1;
			
			if (currentEnd == totalLength - 1){
				break;
			}
		}
		
		currentStart = 0;		
		final List<BufferedReader> toret = new LinkedList<BufferedReader>();
		for (Long end: chunks ){
			final InputStream is = new FileSequenceInputStream(fastqs, end, currentStart);
			
			toret.add(new BufferedReader(new InputStreamReader(is)));
			currentStart = end;
		}
		
		return toret;
		
	}

	
	private static long adjustPos(File f, long skip) throws IOException {
		//long previous = skip;
		
		// find the next @ after the skip-th byte in f
		InputStream is = new FileInputStream(f);
		is.skip(skip);
		byte bytes[] = new byte[1024]; //chunk big enough 
		int readed;
		while ((readed= is.read(bytes))!=-1){
			//find @
			String s = new String(bytes, 0, readed);
			
			int index =-1;
			while((index = s.indexOf("\n@"))!=-1){
				//are we in a sequence name line, or in the middle of the quality string?
				//if we are in the sequence name line, after two eols we must find a '+' char
				//if not, we must skip this quality line
				
				String[] tokens = s.substring(index+1).split("\n",3);
				
				if (tokens.length<3){
					throw new RuntimeException("cant find 3 \\n chars in fastq, which is expected in the readed buffer. Buffer is: "+s);
				}
				if (tokens[2].startsWith("+")){
					return skip + index + 1;					
				}else{
					int index2 =  s.substring(index+1).indexOf('\n');
					if (index2 ==-1){
						throw new RuntimeException("cant find \\n char in fastq, which is expected in the readed buffer. Buffer is: "+s);
					}
					return skip + index + 1 + index2;
				}
			}
		}
		return f.length();
	}

	private static int getFileForPosition(long currentEnd, List<Long> fileStarts) {
		//System.out.println(fileStarts);
		int i = 0;
		while(i<fileStarts.size() && fileStarts.get(i)<=currentEnd){
			i++;
		}
		return i-1;
	}
	
	public static void main(String[] args) throws IOException{
		File f1 = File.createTempFile("fastqexample", ".fastq"); f1.deleteOnExit();
		System.out.println("f1: "+f1);
		PrintStream ps1 = new PrintStream(new FileOutputStream(f1));
		
		File f2 = File.createTempFile("fastqexample", ".fastq"); f2.deleteOnExit();
		PrintStream ps2 = new PrintStream(new FileOutputStream(f2));
		System.out.println("f2: "+f2);
		
		File f3 = File.createTempFile("fastqexample", ".fastq"); f3.deleteOnExit();
		PrintStream ps3 = new PrintStream(new FileOutputStream(f3));
		System.out.println("f3: "+f3);
		
		ps1.println("@sequence1name\nsequence1\n+a\nquality1\n@sequence2name\nsequence2\n+\nquality2");
		ps2.println("@sequence3name\nsequence3\n+\nquality3");
		ps3.println("@sequence4name\nsequence4\n+\nquality4\n@sequence5name\nsequence5\n+\nquality5");
		
		ps1.close();
		ps2.close();
		ps3.close();
		
		List<BufferedReader> iss = splitfastq(Arrays.asList(new File[]{f1,f2,f3}), 2);
		
		for (BufferedReader is : iss){
			System.out.println("New input stream\n---------------");
		
			String line = null;
			while ((line=is.readLine())!=null){
				System.out.println(line);
			}
		}
		System.out.println("****************************");
		iss = splitfastq(Arrays.asList(new File[]{f1,f2,f3}), 20);
		
		for (BufferedReader is : iss){
			System.out.println("New input stream\n---------------");
			String line = null;
			while ((line=is.readLine())!=null){
				System.out.println(line);
			}
		}
		
		System.out.println("****************************");
		iss = splitfastq(Arrays.asList(new File[]{f1,f2,f3}), 1);
		
		for (BufferedReader is : iss){
			System.out.println("New input stream\n---------------");
			String line = null;
			while ((line=is.readLine())!=null){
				System.out.println(line);
			}
		}
		
		System.out.println("****************************");
		
		InputStream is = new FileSequenceInputStream(Arrays.asList(new File[]{f1,f2,f3}),125,0);
		System.out.println("Last input stream\n---------------");
		BufferedReader reader = new BufferedReader(new InputStreamReader(is));
		String line = null;
		while ((line=reader.readLine())!=null){
			
			System.out.println("line: "+line);
		}
		
	}
	
	
	
}
