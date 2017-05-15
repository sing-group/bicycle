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

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

public abstract class AbstractCommand implements Command {

	protected List<Option> options = new LinkedList<Option>();

	public AbstractCommand() {
		this.options = createOptions();
	}

	@Override
	public List<Option> getOptions() {
		return Collections.unmodifiableList(options);
	}


	protected Option findOption(String name) {
		for (Option option : this.getOptions()) {
			if (option.getParamName().equalsIgnoreCase(name) || option.getShortName().equalsIgnoreCase(name)) {
				return option;
			}
		}
		return null;
	}

	protected abstract List<Option> createOptions(); //factory method

}
