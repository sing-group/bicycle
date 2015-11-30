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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Date;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import es.cnio.bioinfo.bicycle.Project;

public abstract class ProjectCommand extends AbstractCommand{


	@Override
	protected List<Option> createOptions() {
		List<Option> toret = new LinkedList<Option>();
		
		toret.add(new Option("project-directory", "p", 
				"project directory. Use command create-project to create a new project", false, true));
		
		return toret;
	}
	
	@Override
	public final void execute(CLIApplication app, Map<Option, String> parameters) throws Exception {
		Project project = Project.readFromDirectory(new File(parameters.get(this.findOption("p"))));
		writeExecutionLog(app, project);
		executeImpl(app, project, parameters);
		
	}

	public static void writeExecutionLog(CLIApplication app, Project p) throws FileNotFoundException{
		File log = new File(p.getProjectDirectory()+File.separator+"executions.log");
		PrintStream out = new PrintStream(new FileOutputStream(log, true));
		out.println(new Date());
		out.print(app.getApplicationCommand());
		for(String s: app.commandline){
			out.print(" "+s);
		}
		out.println();
		out.println();
		
	}
	protected abstract void executeImpl(CLIApplication app, Project p, Map<Option, String> parameters) throws Exception;

}
