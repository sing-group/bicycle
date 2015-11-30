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

package es.cnio.bioinfo.bicycle.cli;

import java.util.List;
import java.util.Map;

import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Sample;
import es.cnio.bioinfo.bicycle.operations.SampleBisulfitation;

public class ReadsBisulfitationCommand extends ProjectCommand {

	@Override
	public String getName() {
		return "reads-bisulfitation";
	}

	@Override
	public String getDescription() {
		return "Performs reads in-silico bisulfitation (CtoT)";
	}

	@Override
	public void executeImpl(CLIApplication app, Project project, Map<Option, String> parameters) throws Exception {
		
		
		for (Sample s: project.getSamples()){
			SampleBisulfitation sb = new SampleBisulfitation(s);
			sb.computeSampleBisulfitation(parameters.containsKey(this.findOption("b")));
		}
	}

	@Override
	protected List<Option> createOptions() {
		List<Option> toret = super.createOptions();
		
		/*toret.add(new Option("remove-bad-bisulfited", "r", 
				"remove unconverted reads (those were the bisulfite " +
				"treatment failed, following the rule applied in Lister et al., Nature 2009)", true, false));
		*/
		toret.add(new Option("remove-unconverted-barcodes", "b", "remove reads with unconverted barcodes", true, false));
		
		return toret;
	}

}
