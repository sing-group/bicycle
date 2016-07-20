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

import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;

public class BicycleApplication extends CLIApplication {

	@Override
	protected String getDescription() {
			Scanner sc = new Scanner(getClass().getResourceAsStream("header.txt"));
			StringBuilder builder = new StringBuilder();
			while (sc.hasNextLine()) {
				builder.append(sc.nextLine()+"\n");
			}
			return builder.toString();
	}
	
	@Override
	public List<Command> buildCommands() {
		List<Command> commands = new LinkedList<Command>();
		commands.add(new CreateProjectCommand());
		commands.add(new ReferenceBisulfitationCommand());
		commands.add(new ReadsBisulfitationCommand());
		commands.add(new ReferenceIndexingCommand());
		commands.add(new BowtieAlignmentCommand());
		commands.add(new MethylationAnalysisCommand());
		commands.add(new DifferentialMethylationAnalysisCommand());		
		return commands;
	}

	@Override
	protected String getApplicationName() {
		return "bicycle";
	}
	
	@Override
	protected String getApplicationCommand() {
		return "bicycle";
	}
	
	
}
