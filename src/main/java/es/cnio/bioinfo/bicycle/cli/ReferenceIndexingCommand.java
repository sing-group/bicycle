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
import es.cnio.bioinfo.bicycle.Reference;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment;

public class ReferenceIndexingCommand extends ProjectCommand {

	@Override
	public String getName() {
		return "reference-index";
	}

	@Override
	public String getDescription() {
		return "Tells Bowtie to build indexes for both references, CtoT and GtoA";
	}

	@Override
	public void executeImpl(CLIApplication app, Project project, Map<Option, String> parameters) throws Exception {

		BowtieAlignment al = new BowtieAlignment(project);

		int v = Integer.parseInt(parameters.get(this.findOption("v")));
		if (v != 1 && v != 2) {
			throw new IllegalArgumentException("bowtie version must be 1 or 2");
		}

		for (Reference ref : project.getReferences()) {
			if (v == 1) {
				al.buildBowtieIndex(ref);
			} else { //bowtie 2
				al.buildBowtie2Index(ref);
			}
		}
	}

	@Override
	protected List<Option> createOptions() {
		List<Option> toret = super.createOptions();

		toret.add(new DefaultValuedOption("bowtie-version", "v",
				"bowtie version to use (valid options are 1 or 2)"
				, "2"));

		return toret;
	}

}
