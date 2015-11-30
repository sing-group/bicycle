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

public class DefaultBisulfiteError implements BisulfiteError{

	
	public int totalReads=0;
	public int errorReads=0;
	
	public Context context;
	public Strand strand;
	
	public DefaultBisulfiteError(Context context, Strand strand) {
		super();
		this.context = context;
		this.strand = strand;
	}

	public void add(int totalReads, int errorReads){
		if (errorReads > totalReads || totalReads<0 || errorReads <0){
			throw new IllegalArgumentException("total reads cannot be greater than errorReads and both must be positive");
		}
		this.totalReads += totalReads;
		this.errorReads += errorReads;
	}
	
	@Override
	public double getError(){
		if (totalReads == 0){
			return 0d;
		}
		return (double)errorReads/(double)totalReads;
	}
	
	
	public int getTotalReads() {
		return totalReads;
	}
	public int getErrorReads() {
		return errorReads;
	}

	public Context getContext() {
		return context;
	}

	public Strand getStrand() {
		return strand;
	}
	
	
}
