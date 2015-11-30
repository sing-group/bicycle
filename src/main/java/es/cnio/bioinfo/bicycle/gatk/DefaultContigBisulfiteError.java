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

public class DefaultContigBisulfiteError implements ContigBisulfiteError{
	

	public HashMap<Strand, HashMap<Context, DefaultBisulfiteError>> perStrandErrors = new HashMap<Strand, HashMap<Context, DefaultBisulfiteError>>();
	
	public DefaultContigBisulfiteError() {
		
		for (Strand strand: Strand.values()){
			HashMap<Context, DefaultBisulfiteError> errorMap = new HashMap<Context, DefaultBisulfiteError>();
			perStrandErrors.put(strand, errorMap);
			for (Context context: Context.values()){
				errorMap.put(context, new DefaultBisulfiteError(context, strand));
			}
		}
	
	}
	
	public void addError(Strand strand, Context context, int reads, int errorReads){
		HashMap<Context, DefaultBisulfiteError> errorMap = perStrandErrors.get(strand);
		
		errorMap.get(context).add(reads, errorReads);
		
	}
	
	public DefaultBisulfiteError getError(Strand strand, Context context){
		HashMap<Context, DefaultBisulfiteError> errorMap = perStrandErrors.get(strand);
		return errorMap.get(context);
	}
	
	public String toString(){
		StringBuffer toretb = new StringBuffer();
		for (Strand strand: new Strand[]{Strand.WATSON, Strand.CRICK}){
			toretb.append(strand+"={");
			boolean first=true;
			for (Context context: Context.values()){
				if (!first){
					toretb.append(", ");
				}else first = false;
				toretb.append(context+"="+this.perStrandErrors.get(strand).get(context).getError()+" ("+this.perStrandErrors.get(strand).get(context).getErrorReads()+"/"+this.perStrandErrors.get(strand).get(context).getTotalReads()+")");
			}
			toretb.append("}\n");
			
		}
		return toretb.toString();
	}
}
