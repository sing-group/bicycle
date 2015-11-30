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
import java.io.IOException;

import es.cnio.bioinfo.bicycle.operations.BowtieAlignment.Strand;


/**
 * Takes two BufferedReaders from two fastq files (mate1 mate2) and provides
 * a reader compatible to feed a bowtie process in paired-end mode with the --12 - 
 * input.
 * 
 * This input must be one line per mate following the format:
 * read name |tab| sequence mate 1 |tab| quality mate 1 |tab| sequence mate 2 |tab| quality mate 2 
 * 
 * 
 * @author lipido
 *
 */
class PairedEndBowtieReader extends BufferedReader{

	private BufferedReader mate1;
	private BufferedReader mate2;
	private Strand strand;
	public PairedEndBowtieReader(BufferedReader mate1fastq, BufferedReader mate2fastq, Strand strand) {
		super(mate1fastq); //we must call super constructor, but this is only to avoid npe
		this.mate1 = mate1fastq;
		this.mate2 = mate2fastq;
		this.strand = strand;
	}
	@Override
	public String readLine() throws IOException {
		//try to read four lines from each stream
		String m1ReadName = mate1.readLine();
		if (m1ReadName == null) return null;
		
		String m1Sequence = mate1.readLine();
		if (m1Sequence == null) return null;
		
		//throw away
		if (mate1.readLine()==null) return null;
		
		String m1Qual = mate1.readLine();
		if (m1Qual == null) return null;
		
		String m2ReadName = mate2.readLine();
		if (m2ReadName == null) return null;
		
		String m2Sequence = mate2.readLine();
		if (m2Sequence == null) return null;
		
		//throw away
		if (mate2.readLine()==null) return null;
		
		String m2Qual = mate2.readLine();
		if (m2Qual == null) return null;
		
		String[] m1NameTokens = m1ReadName.split("[|][|]");
		String[] m2NameTokens = m2ReadName.split("[|][|]");
		
		
		StringBuilder toretb = new StringBuilder();
		
		if (this.strand == Strand.WATSON){
			toretb.append(m1NameTokens[0].substring(1))./*append("WATSON").*/
			append("||").
			append(m1NameTokens[1]). 
			append("||").
			//append(getReverseComplementary(m2NameTokens[1])).
			append(m2NameTokens[1]).
			append("\t").
			append(m1Sequence). 
			append("\t").
			append(m1Qual). 
			append("\t").
			//append(getReverseComplementary(m2NameTokens[1]).replaceAll("C", "T")). //we are performing here the in-silico bisulfitation :-(
			append(m2NameTokens[1].replaceAll("G", "A")). //we are performing here the in-silico bisulfitation :-(
			
			append("\t").
//			append(reverse(m2Qual));
			append(m2Qual);
			
			System.err.println("a WATSON input: "+toretb.toString());
		}else{	
			
			toretb.append(m1NameTokens[0].substring(1))./*append("CRICK").*/
			append("||").
			append(m1NameTokens[1]). 
			append("||").
//			append(getReverseComplementary(m2NameTokens[1])).
			append(m2NameTokens[1]).
			append("\t").
			append(m1Sequence). 
			append("\t").
			append(m1Qual). 
			append("\t").
//			append(getReverseComplementary(m2NameTokens[1]).replaceAll("C", "T")). //we are performing here the in-silico bisulfitation :-(
			append(m2NameTokens[1].replaceAll("G", "A")). //we are performing here the in-silico bisulfitation :-(
			append("\t").
//			append(reverse(m2Qual));
			append(m2Qual);
			System.err.println("a CRICK input: "+toretb.toString());
		}
		return toretb.toString();
	}
	
	private String reverse(String s){
		StringBuilder builder = new StringBuilder();
		for (int i = s.length()-1; i>=0; i--){
			builder.append(s.charAt(i));
		}
		return builder.toString();
	}
	/*
	private int lineno=0;
	@Override
	public String readLine() throws IOException {
		String toret = null;
		
		if (lineno < 4){
			//mate 1
			toret = mate1.readLine();
						
		}else{
			//mate2
			toret = mate2.readLine();
			
		}
		if (lineno==8){
			lineno=0;
		}

		if (toret !=null){
			if (lineno==0){
				toret=toret.replaceAll("[|][|]", "//1||");
			}else if (lineno==4){
				toret=toret.replaceAll("[|][|]", "//2||");
			}

			lineno++;
		}
		return toret;
	}
	
	*/
	private String getReverseComplementary(String sequence){
		StringBuilder complementaria=new StringBuilder("");
		
		// 1ero construyo la complementaria
		for(int i=0;i<sequence.length();i++){
				switch(sequence.charAt(i)){
				case 'A': complementaria.append("T"); break;
				case 'T': complementaria.append("A"); break;
				case 'C': complementaria.append("G"); break;
				case 'G': complementaria.append("C"); break;
				default: complementaria.append(sequence.charAt(i));
			}
		}
		
		//ahora la reversa
		StringBuilder reversa=new StringBuilder("");
		for(int i=complementaria.length()-1;i>-1;i--){
			reversa.append(complementaria.charAt(i));				
		}

		return(reversa.toString());
	}
}