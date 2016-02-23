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

import java.util.HashMap;

import es.cnio.bioinfo.bicycle.MethylationCall;

class MethylationStatistics{
	
	private HashMap<Context, Integer> cCounts = new HashMap<Context, Integer>();
	private HashMap<Context, Integer> mCCounts = new HashMap<Context, Integer>();
	
	private HashMap<Context, Integer> ctReadsCounts = new HashMap<Context, Integer>();
	private HashMap<Context, Integer> mcReadsCounts = new HashMap<Context, Integer>();
	
	private int cCount = 0;
	private int mCCount = 0;
	
	private int ctReadsCount = 0;
	private int mcReadsCount = 0;
	
	private int corrected = 0;
	public MethylationStatistics() {
		for (Context c : Context.values()){
			cCounts.put(c, 0);
			mCCounts.put(c,0);
			
			ctReadsCounts.put(c,  0);
			mcReadsCounts.put(c, 0);
		}
	}
	public void add(MethylationCall call){
		cCount ++;
		cCounts.put(call.getContext(), cCounts.get(call.getContext())+1);
		
		ctReadsCount += call.getCTdepth();
		ctReadsCounts.put(call.getContext(), ctReadsCounts.get(call.getContext())+call.getCTdepth());
		
		mcReadsCount += call.getCytosines();
		mcReadsCounts.put(call.getContext(), mcReadsCounts.get(call.getContext())+call.getCytosines());
		
		if (call.isMethylated()){
			if (call.isCorrectedFromNonCG()){
				corrected ++;
			}
			mCCount ++;
			mCCounts.put(call.getContext(), mCCounts.get(call.getContext())+1);
		}
		
	}
	
	public String toString(){
		StringBuffer toret = new StringBuffer();
		
		
		toret.append("Called methylcytosines (pval<cutoff)\n");
		toret.append(" total: "+mCCount+"/"+cCount+" ("+(double)mCCount/(double)cCount+")\n");
		
		toret.append(" per context called methylcytosines: ");
		for (Context c : Context.values()){
			toret.append(" "+c+":"+(double)mCCounts.get(c)/(double) mCCount);
		}
		toret.append("\n");
		
		for (Context c : Context.values()){
			toret.append(" "+c+" called methylcytosines: "+mCCounts.get(c)+"/"+cCounts.get(c)+" ("+(double)mCCounts.get(c)/(double) cCounts.get(c)+")\n");
		}
		toret.append("Methylation Levels:\n");
		for (Context c : Context.values()){
			toret.append(" "+c+": "+mcReadsCounts.get(c)+"/"+ctReadsCounts.get(c)+" ("+(double)mcReadsCounts.get(c)/(double) ctReadsCounts.get(c)+")\n");
		}
		toret.append("non-CG corrections: "+this.corrected+"\n");
		return toret.toString();
	}
}
public class GlobalMethylationStatistics{

	private MethylationStatistics globalStatistics = new MethylationStatistics();
	
	private HashMap<Strand, MethylationStatistics> perStrand = new HashMap<Strand, MethylationStatistics>();

	public GlobalMethylationStatistics() {
		for (Strand s : Strand.values()){
			perStrand.put(s, new MethylationStatistics());
		}
	}
	
	public void add(MethylationCall call){
		globalStatistics.add(call);
		
		perStrand.get(call.getStrand()).add(call);
		
	}
	@Override
	public String toString() {
		StringBuffer toret = new StringBuffer();
		
		toret.append("---- GLOBAL --------\n");
		toret.append(globalStatistics+"\n");
		for (Strand s: Strand.values()){
			toret.append("---- "+s+" --------\n");
			toret.append(perStrand.get(s)+"\n");
		}
		return toret.toString();
	}
	
}

