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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Properties;
import java.util.UUID;
import java.util.logging.Logger;

public class Utils {

	private static final Logger logger = Logger.getLogger(Utils.class.getName());
	private static File tmpDir = new File(System.getProperty("java.io.tmpdir"));
	public static File generateTempDirName(String prefix){
		//File tmpDir = new File(System.getProperty("java.io.tmpdir"));
		logger.info("creating temp SO temp dir is "+tmpDir);
		
		return new File(tmpDir.getAbsoluteFile()+File.separator+prefix+UUID.randomUUID().toString());
		
	}
	
	public static void deleteDirOnJVMExit(final File dir) {
		Runtime.getRuntime().addShutdownHook(new Thread() {
			@Override
			public void run() {
				deleteDir(dir);
			}
		});
	}
	public static boolean deleteDir(File dir) {
	    if (dir.isDirectory()) {
	        String[] children = dir.list();
	        for (int i=0; i<children.length; i++) {
	            boolean success = deleteDir(new File(dir, children[i]));
	            if (!success) {
	                return false;
	            }
	        }
	    }

	    // The directory is now empty so delete it
	    logger.info("Deleting: "+dir);
	    return dir.delete();
	}
	
	public static File touchFile(File directory, String name) throws FileNotFoundException{
		File newfile = new File(directory.getAbsolutePath()+File.separator+name);		
		touchFile(newfile);
		return newfile;
	}
	public static void touchFile(File f) throws FileNotFoundException{
		if (!f.exists()){
			new FileOutputStream(f);
		}		
	}
	
	public static void append(File f, String s) throws IOException{
		FileWriter writer = new FileWriter(f);
		writer.write(s);
		writer.close();
	}
	public static String readFile(File f) throws IOException{
		StringBuilder builder = new StringBuilder();
		BufferedReader reader = new BufferedReader(new FileReader(f));
		boolean first = true;
		for (String line = reader.readLine(); line!=null; line=reader.readLine(), first=false){
			if (!first) builder.append("\n");
			builder.append(line);			
		}
		return builder.toString();
	}
	
	
	public static String getBowtiePath() throws IOException{
		return getProperty("bowtie.path");
	}
	public static String getSamtoolsPath() throws IOException{
		return getProperty("samtools.path");
	}
	public static String getReadsDirectory(){
		return getClassesPath()+File.separator+"testdata"+File.separator+"READS";
	}
	
	public static String getPairedEndReadsDirectory(){
		return getClassesPath()+File.separator+"testdata"+File.separator+"READSPAIRED";
	}
	
	public static String getNonDirectionalSingleEndDirectory(){
		return getClassesPath()+File.separator+"testdata"+File.separator+"READSCOKUS";
	}
	
	public static String getNonDirectionalPairedEndDirectory(){
		return getClassesPath()+File.separator+"testdata"+File.separator+"READSPAIREDCOKUS";
	}
	
	public static String getReferenceDirectory(){
		return getClassesPath()+File.separator+"testdata"+File.separator+"REFERENCES";
	}
	
	public static String getSimulatedDataReferenceDirectory(){
		return getClassesPath()+File.separator+File.separator+"sampledata"+File.separator+"ref_genomes";
	}
	
	public static String getSimulatedDataReadsDirectory(){
		return getClassesPath()+File.separator+File.separator+"sampledata"+File.separator+"reads";
	}
	
	public static String getSimulatedNonDirectionalDataReadsDirectory(){
		return getClassesPath()+File.separator+File.separator+"sampledata"+File.separator+"reads-non-directional";
	}
	
	public static String getSimulatedDataPairedReadsDirectory(){
		return getClassesPath()+File.separator+"sampledata"+File.separator+"reads-paired";
	}
	
	public static String getSimulatedNonDirectionalPairedDataReadsDirectory(){
		return getClassesPath()+File.separator+"sampledata"+File.separator+"reads-non-directional-paired";
	}
	
	public static String getBedsDirectory(){
		return getClassesPath()+"testdata"+File.separator+"BED";
	}
	
	private static String getProperty(String prop) throws IOException{
		//first, look in system properties, else look in configuration file
		String value = System.getProperty(prop);
		if (value != null)
			return value;
		
		FileInputStream myInputStream = new FileInputStream(getClassesPath()+"/"+"nativedirectories.properties");
		Properties myProps = new Properties();		
		myProps.load(myInputStream);
		return myProps.getProperty(prop);
	}
	
	public static String getClassesPath() {
		return Utils.class.getProtectionDomain().getCodeSource().getLocation().getPath().toString();
	}
	
	public static void main(String[] args) throws IOException{
		System.out.println(getBowtiePath());
	}
	
	
}
