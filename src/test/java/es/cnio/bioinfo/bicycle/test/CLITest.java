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

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import org.junit.BeforeClass;

public abstract class CLITest {
	protected static ByteArrayOutputStream systemout = new ByteArrayOutputStream();
	protected static ByteArrayOutputStream systemerr = new ByteArrayOutputStream();
	@BeforeClass
	public static void hijackStreams(){
		System.setOut(new PrintStream(new OutputStream(){
			PrintStream out = System.out;
			@Override
			public void write(int b) throws IOException {
				systemout.write(b);
				out.write(b);				
			}
			
		}));
		System.setErr(new PrintStream(new OutputStream(){
			PrintStream err = System.err;
			@Override
			public void write(int b) throws IOException {
				systemerr.write(b);
				err.write(b);				
			}
			
		}));
	}
	
	protected String getStdOut(){
		return new String(CLITest.systemout.toByteArray());
	}
	protected String getStdErr(){
		return new String(CLITest.systemerr.toByteArray());
	}

}
