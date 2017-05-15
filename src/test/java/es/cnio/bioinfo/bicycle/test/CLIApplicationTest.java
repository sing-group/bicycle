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

package es.cnio.bioinfo.bicycle.test;

import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.junit.Before;
import org.junit.Test;

import es.cnio.bioinfo.bicycle.cli.CLIApplication;
import es.cnio.bioinfo.bicycle.cli.Command;
import es.cnio.bioinfo.bicycle.cli.Option;

public class CLIApplicationTest extends CLITest {

	@Before
	public void resetStreams() {
		CLIApplicationTest.systemerr.reset();
		CLIApplicationTest.systemout.reset();
	}


	CLIApplication mockapp = new CLIApplication() {

		@Override
		protected List<Command> buildCommands() {
			List<Command> toret = new LinkedList<Command>();
			toret.add(new Command() {
				Option n = new Option("numerator", "n", "the numerator", false, true);
				Option d = new Option("denominator", "d", "the denominator", false, true);
				Option i = new Option("int", "i", "compute as integers", true, false);

				@Override
				public String getName() {
					return "divide";
				}

				@Override
				public String getDescription() {
					return "divides two numbers";
				}

				@Override
				public List<Option> getOptions() {
					return Arrays.asList(new Option[]{n, d, i});


				}

				@Override
				public void execute(CLIApplication app, Map<Option, String> parameters) {
					float numerator = Float.parseFloat(parameters.get(this.n));
					float denominator = Float.parseFloat(parameters.get(this.d));
					boolean asInt = false;
					if (parameters.containsKey(i)) {
						asInt = true;
					}
					float result = numerator / denominator;
					if (asInt) {
						System.out.println("Result is: " + ((int) result));

					} else {
						System.out.println("Result is: " + result);
					}

				}

			});
			return toret;
		}

		@Override
		protected String getApplicationName() {
			return "Mock application";
		}

		@Override
		protected String getApplicationCommand() {
			return "mock";
		}

	};

	@Test
	public void testMainHelp() {
		mockapp.run(new String[]{});
		assertTrue(getStdErr().contains("divide"));
	}

	@Test
	public void testSpecificHelp() {
		mockapp.run(new String[]{"help", "divide"});
		assertTrue(getStdErr().contains("--numerator"));
	}

	@Test
	public void testRun() {
		mockapp.run(new String[]{"divide", "-n", "5", "-d", "2", "-i"});
		assertTrue(getStdOut().contains("Result is: 2" + System.getProperty("line.separator")));
	}

	@Test
	public void testMandatoryParameters() {
		try {
			mockapp.run(new String[]{"divide", "-n", "5"});
		} catch (IllegalArgumentException e) {
			assertTrue(e.getMessage().contains("is mandatory"));
		}
	}

	@Test
	public void testRequiredValueParameters() {
		try {
			mockapp.run(new String[]{"divide", "-n"});
		} catch (IllegalArgumentException e) {
			assertTrue(e.getMessage().contains("requires"));
		}
	}

}
