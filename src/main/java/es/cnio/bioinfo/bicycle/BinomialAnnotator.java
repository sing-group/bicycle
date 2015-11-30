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
import java.io.IOException;
import java.io.InputStreamReader;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;

public class BinomialAnnotator {

	 public static void main(String[] args) throws IOException, MathException{
		 BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		 
		 double p = Double.parseDouble(args[0]);
		 String line = null;
		 
		 
			
		 while ((line = in.readLine())!=null){
		
				 String[] tokens = line.split("\t");
				 int reads = Integer.parseInt(tokens[5]);
				 int cCount = Integer.parseInt(tokens[4]);
				 BinomialDistribution binomial = new BinomialDistributionImpl(reads, p);
				 double pval = (reads==0)? 1.0d: (1.0d-binomial.cumulativeProbability(cCount-1));
				 if(System.out.checkError()){
					 System.exit(1);
				 }
				 System.out.println(line+"\t"+pval);
			 
		 }
	 }
}
