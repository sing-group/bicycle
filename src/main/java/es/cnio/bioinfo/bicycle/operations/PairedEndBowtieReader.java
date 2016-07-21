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
import java.util.LinkedList;
import java.util.Queue;

import es.cnio.bioinfo.bicycle.Sample;
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
class PairedEndBowtieReader extends ReadsReader{

	private BufferedReader mate1;
	private BufferedReader mate2;
	private boolean isDirectional;
	
	private boolean skipUnconverted = false;
	public PairedEndBowtieReader(Sample sample, BufferedReader mate1fastq, BufferedReader mate2fastq, boolean isDirectional) {
		this(sample, mate1fastq, mate2fastq, isDirectional, false);
	}
	
	public PairedEndBowtieReader(Sample sample, BufferedReader mate1fastq, BufferedReader mate2fastq, boolean isDirectional, boolean skipUnconverted) {
		super(sample, mate1fastq); //we must call super constructor, but this is only to avoid npe
		this.mate1 = mate1fastq;
		this.mate2 = mate2fastq;
		this.isDirectional = isDirectional;
		this.skipUnconverted = skipUnconverted;
	}
	private Queue<String> nextString = new LinkedList<String>();
	
	@Override
	public String readLine() throws IOException {
		if (!nextString.isEmpty()) {
			return nextString.poll();
		}
		String m1ReadName = null;
		String m1Sequence = null;
		String m1Qual = null;
		String m2ReadName = null;
		String m2Sequence = null;
		String m2Qual = null;

		do {
			//try to read four lines from each stream
			m1ReadName = mate1.readLine();
			if (m1ReadName == null) return null;
			
			m1Sequence = mate1.readLine();
			if (m1Sequence == null) return null;
			
			//throw away
			if (mate1.readLine()==null) return null;
			
			m1Qual = mate1.readLine();
			if (m1Qual == null) return null;
			
			m2ReadName = mate2.readLine();
			if (m2ReadName == null) return null;
			
			m2Sequence = mate2.readLine();
			if (m2Sequence == null) return null;
			
			//throw away
			if (mate2.readLine()==null) return null;
			
			m2Qual = mate2.readLine();
			if (m2Qual == null) return null;
		
		} while (shouldSkip(m1ReadName) || shouldSkip(m2ReadName));
		
		/*String[] m1NameTokens = m1ReadName.split("[|][|]");
		String[] m2NameTokens = m2ReadName.split("[|][|]");
		*/
		
		StringBuilder toretb = new StringBuilder();
		
		toretb.append(m1ReadName.substring(1).replace(' ', '_')).
		append("||").
		append(m1Sequence). 
		append("||").
		//append(getReverseComplementary(m2NameTokens[1])).
		append(m2Sequence).
		append("\t").
		append(m1Sequence.replaceAll("C", "T")). 
		append("\t").
		append(m1Qual). 
		append("\t").
		append(m2Sequence.replaceAll("G", "A")). 
		
		append("\t").
		append(m2Qual);
		
		if (!this.isDirectional) {
			StringBuilder nonDirectionalSB = new StringBuilder();
			nonDirectionalSB.append(m1ReadName.substring(1).replace(' ', '_')).
			append("||").
			append(m1Sequence). 
			append("||").
			//append(getReverseComplementary(m2NameTokens[1])).
			append(m2Sequence).
			append("\t").
			append(m1Sequence.replaceAll("G", "A")).
			append("\t").
			append(m1Qual). 
			append("\t").
			append(m2Sequence.replaceAll("C", "T")).
			
			append("\t").
			append(m2Qual);
			
			nextString.offer(nonDirectionalSB.toString());
		}
		
		return toretb.toString();
	}

	private boolean shouldSkip(String readName) {
		return this.skipUnconverted && this.hasUnconvertedBarcode(readName);
	}
	
}