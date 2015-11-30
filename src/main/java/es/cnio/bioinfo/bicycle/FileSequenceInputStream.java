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
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class FileSequenceInputStream extends InputStream{

	private List<FileInputStream> inputStreams = new LinkedList<FileInputStream>(); 
	private List<File> files = new LinkedList<File>();
	private long totalReaded = 0;
	private long currentReaded = 0;
	private FileInputStream current;
	private File currentFile;
	private int readedFiles = 0;
	private long maxRead;
	public FileSequenceInputStream(List<File> arg0, long maxRead, long initialSkip) throws IOException {
		//System.out.println("created seq input stream for "+arg0+". maxread: "+maxRead+" skip: "+initialSkip);
		this.files = arg0;
		for (File f : files){
			this.inputStreams.add(new FileInputStream(f));
		}
		current = inputStreams.get(0);
		currentFile = files.get(0);
		this.maxRead = maxRead;

		long totalLength = 0;
		for (File f: this.files){
			totalLength+=f.length();
		}
		if (totalLength<maxRead){
			maxRead = totalLength;
		}
		
		this.skip(initialSkip);
	}
	
	
	private long currentLength(){
		return this.currentFile.length();
	}
	
	private boolean inLastFile(){
		
		return readedFiles == this.files.size()-1;
	}
	@Override
	public long skip(long n) throws IOException {		
		if (maxRead - totalReaded < n){
			return skip(maxRead - totalReaded);
		}
		
		long remaining = currentLength() - currentReaded;
		if (remaining >= n){
			long skipped = this.current.skip(n);
			totalReaded += skipped;
			currentReaded += skipped;
			return skipped;
		}else{ //we need to jump to the next file, if it exists
			if (inLastFile()){
				long skipped = this.current.skip(n);
				totalReaded += skipped;
				currentReaded += skipped;
				return skipped;
			}else{
				totalReaded += remaining;
				nextFile();
				return this.skip(n-remaining);
				/*long skipped = this.current.skip(n - remaining);
				totalReaded += skipped;
				currentReaded = skipped;
				return remaining + skipped;*/
				
			}
		}
	}
	
	private void nextFile(){
		if (!inLastFile()){
			readedFiles ++;
			this.currentFile = this.files.get(readedFiles);
			this.current = this.inputStreams.get(readedFiles);
			currentReaded = 0;
			
		}else{
			throw new RuntimeException("you should not see this!");
		}
	}
	@Override
	public int read() throws IOException {
		if (this.totalReaded == this.maxRead){
			return -1;
		}
		byte[] b = new byte[1];
		int readed = this.read(b, 0, 1);
		
		if (readed != -1){
			return b[0];
		}else{
			return -1;
		}
	}
	
	@Override
	public int read(byte[] b, int off, int len) throws IOException {
		
		int toret = -1;
		
		if (this.totalReaded == this.maxRead){
			toret = -1;
		}
		else if (maxRead - totalReaded < len){
			toret= read(b, off, (int)(maxRead-totalReaded));
		}
		else{
			long remaining = currentLength() - currentReaded;
			if (remaining >= len){
				int readed = this.current.read(b, off,len);
				totalReaded += readed;
				currentReaded += readed;
				toret = readed;
			}else{
				if (inLastFile()){
					int readed = this.current.read(b, off,len);
					totalReaded += readed;
					currentReaded += readed;
					toret = readed;
				}else{							
					
					int readed = this.current.read(b, off,(int)remaining);
					
					totalReaded += readed;
					currentReaded += readed;
					if (readed == remaining){
						nextFile();
						
						readed += this.read(b, off+readed, len-readed);
					}				
					toret = (int)(readed);
					
				}
			}
		}
		
		return toret;
	}
	
	
	public static void main(String[] args) throws IOException{
		//test 1, try several files
		File file1 = File.createTempFile("sequenceinputstream", ".dat"); file1.deleteOnExit();
		FileOutputStream os1 = new FileOutputStream(file1);
		
		File file2 = File.createTempFile("sequenceinputstream", ".dat"); file2.deleteOnExit();
		FileOutputStream os2 = new FileOutputStream(file2);
		
		File file3 = File.createTempFile("sequenceinputstream", ".dat"); file3.deleteOnExit();
		FileOutputStream os3 = new FileOutputStream(file3);
		
		//write three bytes on file1
		os1.write(1);
		os1.write(2);
		os1.write(3);
		os1.close();
		
		os2.write(4);
		os2.write(5);
		os2.close();
		
		os3.write(6);
		os3.write(7);
		os3.write(8);
		os3.close();
		

		//test 1
		System.out.println("test 1");
		
		FileSequenceInputStream is = new FileSequenceInputStream(Arrays.asList(new File[]{file1,file2,file3}), 7, 4);
		
		int readed;
		
		
		while((readed=is.read())!=-1){
			
			System.out.println(readed);
		}
		
		//test 2
		System.out.println("test 2");
		is = new FileSequenceInputStream(Arrays.asList(new File[]{file1,file2,file3}), 1000, 0);
		
		byte[] bytes = new byte[3];
		
	
		while((readed=is.read(bytes))!=-1){			
			System.out.println(readed+": "+Arrays.toString(bytes));
		}
	
		
		
	}
	

}
