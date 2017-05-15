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

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;


public abstract class CLIApplication {
	private static Logger logger = Logger.getLogger(CLIApplication.class.getSimpleName());

	protected List<Command> commands = buildCommands();
	private HashMap<String, Command> commandsByName = new HashMap<String, Command>();

	protected abstract List<Command> buildCommands(); //factory method

	protected abstract String getApplicationName();

	protected abstract String getApplicationCommand();

	protected String[] commandline;
	protected Command command;
	protected Map<Option, String> parameters;

	protected String getDescription() {
		return "";
	}

	public CLIApplication() {
		for (Command c : commands) {
			commandsByName.put(c.getName().toUpperCase(), c);
		}
	}

	public void run(String[] args) {
		this.commandline = args;

		if (args.length == 0) {
			printHelp();
		} else if (args[0].equalsIgnoreCase("help")) {
			printWellcome();
			if (args.length >= 2) {
				Command command = commandsByName.get(args[1].toUpperCase());
				if (command != null)
					printCommandHelp(command);
				else {
					System.err.println("Command " + args[1] + " not found");
					printHelp();
				}
			}
		} else if (args[0].equalsIgnoreCase("fullhelp")) {
			printWellcome();
			printFullHelp();
		} else { //run the command
			this.command = commandsByName.get(args[0].toUpperCase());

			if (this.command != null) {
				String[] noFirst = new String[args.length - 1];
				System.arraycopy(args, 1, noFirst, 0, noFirst.length);
				try {
					this.parameters = parseCommand(this.command, noFirst);

					this.beforeRun();
					this.command.execute(this, this.parameters);
				} catch (ParsingException e) {
					logger.log(Level.SEVERE, "Error parsing command", e);
					printCommandHelp(this.command);
				} catch (Exception e) {
					logger.log(Level.SEVERE, "Error during execution", e);
					e.printStackTrace();
				}
			} else {
				logger.log(Level.SEVERE, "Command " + args[0] + " not found");
				printHelp();
			}

		}
	}


	private Option findOption(Command command, String name) {
		for (Option option : command.getOptions()) {
			if (option.getParamName().equalsIgnoreCase(name) || option.getShortName().equalsIgnoreCase(name)) {
				return option;
			}
		}
		return null;
	}

	private HashMap<Option, String> parseCommand(Command command, String[] arguments) throws ParsingException {
		HashMap<Option, String> values = new HashMap<Option, String>();

		command.getOptions();
		Option currentOption = null;

		for (String token : arguments) {
			if (token.startsWith("-")) {

				if (currentOption != null && !values.containsKey(currentOption) && currentOption.requiresValue()) {
					throw new IllegalArgumentException("option " + currentOption.getParamName() + " requires a value");
				} else if (currentOption != null && !values.containsKey(currentOption) && !currentOption.requiresValue
						()) {
					values.put(currentOption, null);
				}

				String optionName = null;
				if (token.charAt(1) == '-') {
					//starts with --
					optionName = token.substring(2);
				} else {
					//starts with -
					optionName = token.substring(1);
				}
				Option option = findOption(command, optionName);
				if (option == null) {
					throw new ParsingException("option " + optionName + " not found");
				} else {
					currentOption = option;
				}
			} else {
				if (currentOption == null) {
					throw new ParsingException("unable to parse. You should specify an option before a value");
				} else {
					if (values.containsKey(currentOption)) {
						throw new ParsingException("option " + currentOption.getParamName() + " was already " +
								"specified");
					} else {
						values.put(currentOption, token);
					}
				}
			}

		}
		//if there is the last option with no required value
		if (currentOption != null && !values.containsKey(currentOption) && currentOption.requiresValue()) {
			throw new ParsingException("option " + currentOption.getParamName() + " requires a value");
		} else if (currentOption != null && !values.containsKey(currentOption) && !currentOption.requiresValue()) {
			values.put(currentOption, null);
		}

		//validate mandatory arguments and put defaults
		for (Option option : command.getOptions()) {
			if (!option.isOptional() && !values.containsKey(option) && (!(option instanceof DefaultValuedOption))) {
				throw new ParsingException("option " + option.getParamName() + " is mandatory");
			}

			//put default-values if not specified before
			if (!values.containsKey(option) && ((option instanceof DefaultValuedOption))) {
				values.put(option, ((DefaultValuedOption) option).getDefaultValue());
			}
		}
		return values;
	}

	private void printCommandHelp(Command command) {
		System.err.println("Command " + command.getName());
		printUsage(command);
		for (Option option : command.getOptions()) {
			System.err.println(
					"\t--" + option.getParamName() + "/-" + option.getShortName()
							+ "\n\t\t" + option.getDescription()
							+ ((option instanceof DefaultValuedOption) ? " (default: " + ((DefaultValuedOption)
							option).getDefaultValue() + ")" : ""));
		}
	}

	private void printUsage(Command command) {
		System.err.print("usage: " + this.getApplicationCommand() + " " + command.getName());
		for (Option option : command.getOptions()) {
			if (option.isOptional()) {
				System.err.print(" [");
			} else {
				System.err.print(" ");
			}
			System.err.print("-" + option.getShortName());
			if (option.requiresValue()) {
				System.err.print(" <" + option.getParamName() + ">");
			}

			if (option.isOptional()) {
				System.err.print("]");
			}
		}
		System.err.println();

	}

	private void printHelp() {
		printWellcome();
		System.err.println("usage: " + this.getApplicationCommand() + " <command> [options]");

		System.err.println("where <command> is one of:");
		for (Command option : commands) {
			System.err.println("\t" + option.getName() + "\n\t\t" + option.getDescription());
		}

		System.err.println("Write '" + this.getApplicationCommand() + " help <command>' to see command-specific help");
	}

	private void printFullHelp() {
		System.err.println("usage: " + this.getApplicationCommand() + " <command> [options]");

		System.err.println("where <command> is one of:");
		for (Command option : commands) {
			System.err.println("\t" + option.getName() + "\n\t\t" + option.getDescription());
		}
		System.err.println();
		for (Command command : commands) {
			printCommandHelp(command);
			System.err.println();
		}
	}

	private void printWellcome() {
		System.err.println("Welcome to " + this.getApplicationName());
		System.err.println(this.getDescription());
	}

	public String[] getCommandline() {
		return commandline;
	}

	//hooks
	public void beforeRun() {
	}

}

class ParsingException extends Exception {
	public ParsingException(String message) {
		super(message);
	}
}
