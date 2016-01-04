/*

Copyright 2012 Daniel Gonzalez Pe��a, Osvaldo Gra��a


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

import java.io.PrintStream;
import java.util.LinkedList;
import java.util.List;
public class MethylationCall {
	private Context context;
	private double pval;
	private int depth;
	private int CTdepth;
	private String pileup;
	private int cytosines;
	private String contig;
	private long position;
	private Strand strand;
	boolean correctedFromNonCG;
	private boolean addedByCorrection;
	private double cutOff=-1;
	private List<String> annotations;
	private double betaScore;
	
	
	public MethylationCall(String contig, long position, Strand strand, Context context, double pval,
			int depth, int CTdepth, int cytosines, String pileup, boolean correctedFromNonCG, boolean addedByCorrection, List<String> annotations, double betaScore) {
		super();
		this.contig = contig;
		this.position = position;
		this.strand = strand;
		
		this.context = context;
		this.pval = pval;
		this.depth = depth;
		this.CTdepth = CTdepth;
		this.cytosines = cytosines;
		this.pileup = pileup;
		this.correctedFromNonCG = correctedFromNonCG;
		
		this.addedByCorrection = addedByCorrection;
		this.annotations = annotations;
		this.betaScore=betaScore;
	}

	

	public Context getContext() {
		return context;
	}

	public double getPval() {
		return pval;
	}

	public int getDepth() {
		return depth;
	}

	public int getCytosines() {
		return cytosines;
	}

	public String getContig() {
		return contig;
	}

	public long getPosition() {
		return position;
	}
	public Strand getStrand() {
		return strand;
	}
	
	public void correctToCG(){
		this.context = Context.CG;
		this.correctedFromNonCG = true;
	}
		
	public boolean isCorrectedFromNonCG() {
		return correctedFromNonCG;
	}

	public boolean isAddedByCorrection() {
		return addedByCorrection;
	}

	public void setCutOff(double cutOff){
		this.cutOff = cutOff;
	}
	
	public boolean isMethylated(){
		return this.pval < this.cutOff;
	}
	
	public double getCutOff() {
		return cutOff;
	}

	public String getPileup() {
		return pileup;
	}
	
	public int getCTdepth() {
		return CTdepth;
	}

	//added (osvaldo, 3jan2016)	
	public double getBetaScore() {
		return betaScore;
	}
	
	public static String getMarshallHeader(){
		String sep="\t";
        //return "#SEQUENCE"+sep+"POS"+sep+"STRAND"+sep+"CONTEXT"+sep+"DEPTH"+sep+"CT.DEPTH"+sep+"CYTOSINE.COUNT"+sep+"PILEUP"+sep+"PVAL"+sep+"CORRECTED_FROM_NON_CG"+sep+"ADDED_BY_CORRECTION";
		        
		//modified (osvaldo, 31dec2015)
        return "#SEQUENCE"+sep+"POS"+sep+"STRAND"+sep+"CONTEXT"+sep+"DEPTH"+sep+"CT.DEPTH"+sep+"CYTOSINE.COUNT"+sep+"BETA.SCORE"+sep+"PILEUP"+sep+"PVAL"+sep+"CORRECTED_FROM_NON_CG"+sep+"ADDED_BY_CORRECTION";
        
	}

	@Override
	public String toString() {
		return this.getContig()+"\t"+this.getPosition()+"\t"+this.getStrand()+"\t"+this.getContext()+"\t"+this.getDepth()+"\t"+this.getCTdepth()+"\t"+this.getCytosines()+"\t"+this.getPileup()+"\t"+this.getPval()+"\t"+this.correctedFromNonCG+"\t"+this.addedByCorrection;	               
	}
	
	public static MethylationCall unmarshall(String line){
		
		String[] tokens = line.split("\t");
		List<String> annotations=new LinkedList<String>();
		//if (tokens.length>=11)
		
		//modified (osvaldo, 3jan2016): required after adding betaScore in between (increments one position)
		int tokenPos=12;
		if (tokens.length>=tokenPos)	
			for (int i = tokenPos;i<tokens.length;i++){
				annotations.add(tokens[i]);
			}
			
		
		//return new MethylationCall(tokens[0], Long.parseLong(tokens[1]), Strand.valueOf(tokens[2]), Context.valueOf(tokens[3]), Double.parseDouble(tokens[8]), Integer.parseInt(tokens[4]), Integer.parseInt(tokens[5]), Integer.parseInt(tokens[6]), tokens[7], Boolean.valueOf(tokens[9]), Boolean.valueOf(tokens[10]), annotations );
		
		//modified (osvaldo, 31dec2015)
		return new MethylationCall(tokens[0], Long.parseLong(tokens[1]), Strand.valueOf(tokens[2]), Context.valueOf(tokens[3]), Double.parseDouble(tokens[9]), Integer.parseInt(tokens[4]), Integer.parseInt(tokens[5]), Integer.parseInt(tokens[6]), tokens[8], Boolean.valueOf(tokens[10]), Boolean.valueOf(tokens[11]), annotations, Double.parseDouble(tokens[7]) );
	}
	
	public List<String> getAnnotations() {
		return annotations;
	}
	
	public String marshall(){
/*		String annotations="";
        for (String annotation: this.annotations){           
            annotations+="\t"+annotation;
        }
       
        return this.getContig()+"\t"+this.getPosition()+"\t"+this.getStrand()+"\t"+this.getContext()+"\t"+this.getDepth()+"\t"+this.getCTdepth()+"\t"+this.getCytosines()+"\t"+this.getPileup()+"\t"+this.getPval()+"\t"+this.correctedFromNonCG+"\t"+this.addedByCorrection+annotations;
*/
       
        //modified (osvaldo, 31dec2015)
		//this method writes in '*methylation' files
        StringBuffer annotations=new StringBuffer();
        for (String annotation: this.annotations){           
            annotations.append("\t").append(annotation);
        }       

        StringBuffer sb=new StringBuffer();
        sb.append(this.getContig());
        sb.append("\t");
        sb.append(this.getPosition());
        sb.append("\t");
        sb.append(this.getStrand());
        sb.append("\t");
        sb.append(this.getContext());
        sb.append("\t");
        sb.append(this.getDepth());
        sb.append("\t");
        sb.append(this.getCTdepth());
        sb.append("\t");
        sb.append(this.getCytosines());
        sb.append("\t");       
        sb.append(this.getBetaScore());
        sb.append("\t");
        sb.append(this.getPileup());
        sb.append("\t");
        sb.append(this.getPval());
        sb.append("\t");
        sb.append(this.correctedFromNonCG);
        sb.append("\t");
        sb.append(this.addedByCorrection);
        sb.append(annotations);
               
        return sb.toString();
	}
	
	public void marshall(PrintStream out){
		out.print(marshall());
	}
	
	
}
