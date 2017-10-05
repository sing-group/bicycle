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

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.BufferUnderflowException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import es.cnio.bioinfo.bicycle.gatk.Strand;
import net.sf.picard.filter.SamRecordFilter;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.SamLocusIterator;
import net.sf.picard.util.SamLocusIterator.LocusInfo;
import net.sf.picard.util.SamLocusIterator.RecordAndOffset;
import net.sf.samtools.SAMFileReader;


/**
 * Performs certain useful operations.
 *
 * @author Osvaldo Gra&ntilde;a, Bioinformatics Unit (CNIO)<br>Copyright 2010 Osvaldo Gra&ntilde;a<br>Distributed
 * under the terms of the GNU General Public License. (Free Software Foundation)
 * @version 2.0, december 2006
 */
public final class Tools {
//Version 1.0, july 2003, EVAtools
// Version 2.0, december 2006


	/**
	 * @param d Receives a double number and returns a string with the number in a format like the printf %5.3f of C++
	 */
	public static String changeFormatOfNumbers(double d) {
		NumberFormat nf = NumberFormat.getInstance();
		nf.setMaximumFractionDigits(2);
		nf.setMinimumFractionDigits(2);

		// compruebo si me entra un NaN, entonces lo cambio por cero
		Double D = new Double(d);
		if (D.isNaN()) d = 0;

		// a veces ocurre que imprime la cadena con ',' en lugar de con '.' para indicar la separación de los
		// decimales,
		// esto no debe ocurrir, por tanto me aseguro de ello
		String cadena = nf.format(d);
		String parseada = new String("");

		for (int i = 0; i < cadena.length(); i++) {
			// si no es una coma lo copio
			if (cadena.charAt(i) != ',') parseada += cadena.charAt(i);
				// si es una coma la cambio por un punto
			else parseada += '.';
		}

		return (parseada);
	}

	/**
	 * Deletes a file
	 */
	public static void deleteFile(String filePath, String file) {
		File f = new File(filePath, file);
		try {
			f.delete();
		} catch (NullPointerException e) {
			System.out.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		}
	}

	/**
	 * Deletes a file. The input argument includes also the complete path to the file
	 */
	public static void deleteFile(String file) {
		File f = new File(file);
		try {
			f.delete();
		} catch (NullPointerException e) {
			System.out.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		}
	}

	/**
	 * Deletes previous results that have expired
	 */
	public static void deletePreviousResults(String directory, PrintWriter p) {
		// recojo la fecha actual y la convierto al formato con el que salvo los archivos
		String fechaYHoraActuales = Tools.getCurrentDate();
		String[] aux = fechaYHoraActuales.split("GMT");
		fechaYHoraActuales = aux[0];
		fechaYHoraActuales = fechaYHoraActuales.replace(' ', '_');
		fechaYHoraActuales = fechaYHoraActuales.replace('+', '_');
		fechaYHoraActuales = fechaYHoraActuales.replace(':', '_');
		String[] romper = fechaYHoraActuales.split("[_]");
		String horaActualString = romper[3];
		Integer IHoraActual = new Integer(horaActualString);
		int iHoraActual = IHoraActual.intValue();
		String minutosActualesString = romper[4];
		Integer IMinutosActuales = new Integer(minutosActualesString);
		int iMinutosActuales = IMinutosActuales.intValue();
		String DiaActualString = romper[2];
		Integer IDiaActual = new Integer(DiaActualString);


		int iDiaActual = IDiaActual.intValue();


		// leo el directorio donde guardo los resultados
		String[] ls = Tools.ls(directory);
		int i;
		for (i = 0; i < ls.length; i++) {
			// pongo un try-catch porque se puede generar una excepción si el archivo que encuentra no se rompe en las
			// partes
			// esperadas al array, como es el caso del archivo 'index.html', pero como ese archivo (u otros distintos
			// que pueda haber)
			// no he de borrarlo, simplemente capturo la excepción y si se produce ésta simplemente sigo a otro
			// fichero, no requiere
			// avisos de estado ni nada.
			try {
				String archivo = ls[i];
				String[] partes = archivo.split("[_]");
				String horaFicheroString = partes[3];
				Integer IHoraFichero = new Integer(horaFicheroString);
				int iHoraFichero = IHoraFichero.intValue();
				String minutosFicheroString = partes[4];
				Integer IMinutosFichero = new Integer(minutosFicheroString);
				int iMinutosFichero = IMinutosFichero.intValue();
				String DiaFicheroString = partes[2];
				Integer IDiaFichero = new Integer(DiaFicheroString);
				int iDiaFichero = IDiaFichero.intValue();

				// si el fichero es del día anterior o más
				if (iDiaFichero < iDiaActual) {
					File f = new File(directory, ls[i]);
					try {
						f.delete();
					} catch (NullPointerException e) {
						System.out.println("\nexception: " + e.getMessage() + "\n");
						e.printStackTrace();
					}
				}
				// si el fichero es de la hora anterior (me da igual el día), ha expirado y lo borro
				else if (iHoraFichero < iHoraActual) {
					File f = new File(directory, ls[i]);
					try {
						f.delete();
					} catch (NullPointerException e) {
						System.out.println("\nexception: " + e.getMessage() + "\n");
						e.printStackTrace();
					}
				}
				// si el fichero ha expirado en minutos, lo borro (me da igual la hora o el día)
				else if (iMinutosFichero + 15 < iMinutosActuales) {
					File f = new File(directory, ls[i]);
					try {
						f.delete();
					} catch (NullPointerException e) {
						System.out.println("\nexception: " + e.getMessage() + "\n");
						e.printStackTrace();
					}
				}
			} catch (ArrayIndexOutOfBoundsException e) {
				System.out.println("\nexception: " + e.getMessage() + "\n");
				e.printStackTrace();
			}
		} //for(i=0;i<ls.length;i++)

	}

	/**
	 * Returns if the file 'child' exists (file or directory)
	 */
	public static boolean doesThisFileExist(String child) {
		boolean result = false;
		File f = new File(child);
		result = f.exists();
		return (result);
	}

	/**
	 * Executes a subprocess with an output stream
	 */
	public static int executeProcessWithOutputStream(String StringForOutputStream, String command, PrintWriter p) {
		int resultado = 0;

		try {
			byte[] bytesArray = StringForOutputStream.getBytes();
			Runtime runtime = Runtime.getRuntime();
			Process process = runtime.exec(command);
			BufferedOutputStream bos = new BufferedOutputStream(process.getOutputStream());
			bos.flush();
			bos.write(bytesArray);


			bos.flush();
			bos.close();

			process.waitFor();
			resultado = process.exitValue();

		} catch (java.lang.InterruptedException e) {
			System.out.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		} catch (java.io.IOException e) {
			System.out.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		}

		return (resultado);
	}

	/**
	 * Executes a subprocess with an output stream
	 */
	public static int executeProcessWithOutputStream(String StringForOutputStream, String command) {
		int resultado = 0;

		try {
			byte[] bytesArray = StringForOutputStream.getBytes();
			Runtime runtime = Runtime.getRuntime();
			Process process = runtime.exec(command);
			BufferedOutputStream bos = new BufferedOutputStream(process.getOutputStream());
			bos.flush();
			bos.write(bytesArray);

			//resultado=process.exitValue();
			bos.flush();
			bos.close();

			process.waitFor();
			resultado = process.exitValue();

		} catch (java.lang.InterruptedException e) {
			System.out.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		} catch (java.io.IOException e) {
			System.out.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		}

		return (resultado);
	}

	/**
	 * Executes a subprocess with an output stream
	 */
	public static int executeProcessWithInputStream(String command) {
		// NO LO USO PORQUE TARDA MUCHISIMO
		int resultado = 0;

		try {
			Runtime runtime = Runtime.getRuntime();
			Process process = runtime.exec(command);
			InputStream inStream = process.getInputStream();
			InputStreamReader inReader = new InputStreamReader(inStream);
			BufferedReader in = new BufferedReader(inReader);

			String data = "";
			String output;

			while (true) {
				output = in.readLine();
				if (output == null) break;
				else data += output;
			}

			in.close();
			inReader.close();
			inStream.close();

			System.out.println("datos:\n" + data);

			process.waitFor();
			resultado = process.exitValue();

		} catch (java.lang.InterruptedException e) {
			System.out.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		} catch (java.io.IOException e) {
			System.out.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		}

		return (resultado);
	}

	/**
	 * Executes a subprocess without waiting for its complete execution
	 */
	public static int executeProcessNoWait(String string) {
		int resultado = 0;

		try {
			Runtime runtime = Runtime.getRuntime();
			Process process = runtime.exec(string);
			resultado = process.exitValue();
		} catch (java.io.IOException e) {
			System.out.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		}

		return (resultado);
	}

	/**
	 * Executes a subprocess waiting for its complete execution
	 */
	public static int executeProcessWait(String[] string) {
		int resultado = 0;

		try {
			Runtime runtime = Runtime.getRuntime();
			Process process = runtime.exec(string);
			try {
				process.waitFor();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				System.out.println("\nexception: " + e.getMessage() + "\n");
				e.printStackTrace();
			}
			resultado = process.exitValue();
		} catch (java.io.IOException e) {
			System.out.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		}

		return (resultado);
	}

	/**
	 * Executes a subprocess without waiting for its complete execution
	 */
	public static int executeProcessNoWait(String string, PrintWriter p) {
		int resultado = 0;

		try {
			Runtime runtime = Runtime.getRuntime();
			Process process = runtime.exec(string);
			resultado = process.exitValue();
		} catch (java.io.IOException e) {
			System.out.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		}

		return (resultado);
	}

	/**
	 * Executes a subprocess waiting for its complete execution
	 */
	public static int executeProcessWait(String string, PrintWriter p) {
		int resultado = 0;

		try {
			Runtime runtime = Runtime.getRuntime();
			Process process = runtime.exec(string);

			try {
				process.waitFor();
				resultado = process.exitValue();
				if (resultado != 0) {
					p.print(process.getErrorStream());
				}
			} catch (java.lang.InterruptedException e) {
				p.println("\nexception: " + e.getMessage() + "\n");
				e.printStackTrace();
			}
		} catch (java.io.IOException e) {
			p.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		}

		return (resultado);
	}

	/**
	 * Executes a subprocess waiting for its complete execution
	 */
	public static int executeProcessWait(String[] command, PrintWriter p) {
		int resultado = 0;

		try {
			Runtime runtime = Runtime.getRuntime();
			Process process = runtime.exec(command);

			try {
				process.waitFor();
				resultado = process.exitValue();
				if (resultado != 0) {
					p.println("executed command:");
					for (int i = 0; i < command.length; i++) p.print(command[i] + " ");
					p.print("\n\nexecution result: ");
					p.println(resultado);
					p.println();

					BufferedReader reader = new BufferedReader(new InputStreamReader(process.getErrorStream()));

					String errorLine;
					while ((errorLine = reader.readLine()) != null) p.println(errorLine);
				}
			} catch (java.lang.InterruptedException e) {
				p.println("\nexception: " + e.getMessage() + "\n");
				e.printStackTrace(p);
			}
		} catch (java.io.IOException e) {
			p.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace(p);
		}

		return (resultado);
	}

	private static void redirect(final InputStream in, final OutputStream out) {
		new Thread() {
			@Override
			public void run() {

				byte[] buffer = new byte[1024];

				int readed = -1;

				try {
					while ((readed = in.read(buffer)) != -1) {
						synchronized (in) {
							out.write(buffer, 0, readed);
						}
					}
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}.start();

	}


	/**
	 * Executes a command line program waiting for its complete execution. The output stream recovers the standard
	 * output of the executed program. It also recovers
	 * the possible errors in the same file<br>
	 * It serves as the linux console command ">", to redirect the program output to a file. it returns 0 if the
	 * command was successfully executed, otherwise returns 1.<br>
	 * In case of problems a message is displayed in the error output.
	 *
	 * @param string
	 * @param stdout
	 * @param stderr
	 * @return int
	 */
	public static int executeProcessWait(String[] string, OutputStream stdout, OutputStream stderr) {
		int resultado = 0;

		try {
			Runtime runtime = Runtime.getRuntime();
			Process process = runtime.exec(string);

			// creo dos hilos de ejecución, uno para el estándar output
			// y otro para el error output
			if (stdout != null) {
				redirect(process.getInputStream(), stdout);
			}
			if (stderr != null) {
				redirect(process.getErrorStream(), stderr);
			}


			try {
				process.waitFor();
				resultado = process.exitValue();
			} catch (java.lang.InterruptedException e) {
				System.out.println("\nexception: " + e.getMessage() + "\n");
				e.printStackTrace();
			}
		} catch (java.io.IOException e) {
			System.out.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		}

		return (resultado);
	}

	public static Process executeProcess(String[] string, OutputStream stdout, OutputStream stderr) {

		try {
			Runtime runtime = Runtime.getRuntime();
			Process process = runtime.exec(string);

			// creo dos hilos de ejecución, uno para el estándar output
			// y otro para el error output
			if (stdout != null) {
				redirect(process.getInputStream(), stdout);
			}
			if (stderr != null) {
				redirect(process.getErrorStream(), stderr);
			}


			return process;
		} catch (java.io.IOException e) {
			throw new RuntimeException(e);

		}


	}

	/**
	 * Executes a command line program waiting for its complete execution. It returns 0 if the command was
	 * successfully executed, otherwise returns 1.<br>
	 * In case of problems a message is displayed in the error output.
	 *
	 * @param command
	 * @return
	 */
	public static int executeProcessWait(String command) {
		try {
			Runtime run = Runtime.getRuntime();
			Process pid = run.exec(command);
			try {
				pid.waitFor();
				if (pid.exitValue() != 0) {
					return (pid.exitValue());
				}

			} catch (java.lang.InterruptedException e) {
				System.out.println("\nexception: " + e.getMessage() + "\n");
				e.printStackTrace();
			}
		} catch (java.io.IOException e) {
			System.out.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		}

		return (0);
	}

	/**
	 * Gets current date
	 */
	public static String getCurrentDate() {
		long l = System.currentTimeMillis();
		Date date = new Date(l);
		return (date.toString());
	}

	/**
	 * Returns the current week of the current year
	 */
	public static int getCurrentWeek() {
		Calendar calendar = Calendar.getInstance();
		int weekNumber = calendar.get(Calendar.WEEK_OF_YEAR);

		// Si la primera semana del año no es aquella que contiene al primer dia del año del primer mes del año,
		// getMinimalDaysInFirstWeek retorna !=1, por tanto el número de semana se devuelve 1
		if (calendar.getMinimalDaysInFirstWeek() != 1 && weekNumber > 52) return (1);
		else if (calendar.getMinimalDaysInFirstWeek() != 1 && weekNumber < 52) return (weekNumber + 1);
		// Si la primera semana del año es aquella que contiene al primer dia del año del primer mes del año,
		// getMinimalDaysInFirstWeek retorna 1, por tanto el número de semana se devuelve tal cual esta
		return (weekNumber);
	}

	/**
	 * Returns the current year
	 */
	public static int getCurrentYear() {
		Calendar calendar = Calendar.getInstance();
		int year = calendar.get(Calendar.YEAR);
		return (year);
	}

	/**
	 * Returns if the directory 'child' is empty
	 */
	public static boolean isThisDirectoryEmpty(String child) {
		boolean resultado = false;
		String datos[] = ls(child);
		if (datos == null || datos.length == 0) resultado = true;
		return (resultado);
	}

	/**
	 * ls of the child directory
	 */
	public static String[] ls(String child) {
		File f = new File(child);
		return (f.list());
	}

	/**
	 * ls of the child directory. All the files in the directory that ends with the filter are returned.
	 *
	 * @param child  directory
	 * @param filter
	 */
	public static String[] ls(String child, String filter) {
		File f = new File(child);
		String[] results = f.list();
		ArrayList<String> files = new ArrayList<String>(0);

		// reduzco la lista inicial a un ArrayList
		for (int i = 0; i < results.length; i++) if (results[i].endsWith(filter)) files.add(results[i]);
		// la recupero de nuevo en un array de cadenas
		results = new String[files.size()];
		for (int i = 0; i < files.size(); i++) results[i] = files.get(i);

		return (results);
	}

	/**
	 * Creates a new directory 'child'
	 */
	public static boolean mkdirs(String child) {
		File f = new File(child);
		return (f.mkdirs());
	}

	/**
	 * Saves data
	 */
	public static void saveData(String path, String fichero, String contenido, PrintWriter p) {
		try {
			byte[] is = contenido.getBytes();
			FileOutputStream fileoutputstream = new FileOutputStream(path + "/" + fichero);
			try {
				fileoutputstream.write(is);
			} catch (IOException e) {
				System.out.println("\nexception: " + e.getMessage() + "\n");
				e.printStackTrace();
			}
			fileoutputstream.close();
		} catch (IOException e) {
			System.out.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		}
	}

	/**
	 * Saves data
	 */
	public static void saveData(String path, String fichero, String contenido) {
		try {
			byte[] is = contenido.getBytes();
			FileOutputStream fileoutputstream = new FileOutputStream(path + "/" + fichero);
			try {
				fileoutputstream.write(is);
			} catch (IOException e) {
				System.out.println("\nexception: " + e.getMessage() + "\n");
				e.printStackTrace();
			}
			fileoutputstream.close();
		} catch (IOException e) {
			System.out.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		}
	}

	/**
	 * Saves data appending
	 */
	public static void saveAppending(String path, String fichero, String contenido) {
		try {
			byte[] is = contenido.getBytes();
			FileOutputStream fileoutputstream = new FileOutputStream(path + "/" + fichero, true);
			try {
				fileoutputstream.write(is);
			} catch (IOException e) {
				System.out.println("\nexception: " + e.getMessage() + "\n");
				e.printStackTrace();
			}
			fileoutputstream.close();
		} catch (IOException e) {
			System.out.println("\nexception: " + e.getMessage() + "\n");
			e.printStackTrace();
		}
	}

	/**
	 * Sets a name for the current process
	 */
	public static String setNameForCurrentProcess(String remoteMachine) {
		String fecha = getCurrentDate();
		// elimno desde GMT en adelante
		String[] aux = fecha.split("GMT");
		fecha = aux[0] + "_GMT";
		fecha = fecha.replace(' ', '_');
		fecha = fecha.replace('+', '_');
		fecha = fecha.replace(':', '_');
		String proceso = fecha + "_" + remoteMachine;
		// retorno el nombre asignado al fichero con los datos del servidor	
		return (proceso);
	}

	/**
	 * Sets a name for the current process
	 */
	public static String setNameForCurrentProcess() {
		String fecha = getCurrentDate();
		// elimno desde GMT en adelante
		String[] aux = fecha.split("GMT");
		fecha = aux[0] + "_GMT";
		fecha = fecha.replace(' ', '_');
		fecha = fecha.replace('+', '_');
		fecha = fecha.replace(':', '_');

		// retorno el nombre asignado al fichero	
		return (fecha);

	}

	/**
	 * Returns date and time
	 */
	public static String getDateAndTime() {
		// recojo la fecha actual y la convierto al formato con el que salvo los archivos
		String fechaYHoraActuales = Tools.getCurrentDate();
		String[] aux = fechaYHoraActuales.split("GMT");
		fechaYHoraActuales = aux[0];
		fechaYHoraActuales = fechaYHoraActuales.replace(' ', '_');
		fechaYHoraActuales = fechaYHoraActuales.replace('+', '_');

		return (fechaYHoraActuales);
	}

	/**
	 * Returns the content of a file as a String. The input argument is the file name (including the path)
	 */
	public static String readFile(String fileName) {
		File readFile = new File(fileName);
		String contenido = new String();


		try {
			FileInputStream fs = new FileInputStream(readFile);
			long l = readFile.length();
			byte[] content = new byte[(int) l];
			fs.read(content);
			fs.close();
			contenido = new String(content);

		} catch (FileNotFoundException e) {
			System.out.println("\n[Error]: unable to find " + fileName);
			e.printStackTrace();
		} catch (IOException e) {
			System.out.println("\n[Error]: " + fileName + "\n" + e.getMessage() + "\n");
			e.printStackTrace();
		}

		return (contenido);
	}

	/**
	 * Copies the content of one file to another and after removes the original one
	 *
	 * @param from
	 * @param to
	 */
	public static void moveOneFileToAnother(File from, File to) {

		if (to.exists()) {
			to.delete();

		}
		if (!from.renameTo(to)) {
			throw new IllegalArgumentException("Delete: " + from.toString() + " deletion failed");
		}
		//empiezo copiando el contenido de from al fichero to
		/*InputStream FROM = null;
		try {
			FROM = new FileInputStream(from);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		OutputStream TO = null;
		try {
			TO = new FileOutputStream(to);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// transfiero los bytes desde FROM hasta TO
		byte[] buf = new byte[1024];
		int len;
		try {
			while((len = FROM.read(buf))>0) TO.write(buf, 0, len);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			FROM.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			TO.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		// ahora borro el fichero original
		boolean success=from.delete();
		if (!success) throw new IllegalArgumentException("Delete: "+from.toString()+" deletion failed");
		*/

	}


	private static ThreadLocal<HashMap<Strand, SamLocusIterator>> previousIterators = new ThreadLocal<>();
	private static ThreadLocal<HashMap<Strand, Interval>> previousIntervals = new ThreadLocal<>();
	private static Map<SamLocusIterator, LocusInfo> currentLocus = Collections.synchronizedMap(new HashMap<SamLocusIterator, LocusInfo>());

	private static final int INTERVAL_WINDOW = 10000000;

	public static List<RecordAndOffset> getReadsForLocus(Strand strand,
														 String contig, int pos, boolean filterDuplicates,
														 List<SAMFileReader> readers,
														 List<SamRecordFilter> filters, char trimbase) {
		List<RecordAndOffset> toret = new LinkedList<RecordAndOffset>();
		for (SAMFileReader reader : readers) {

			//try to find previous iterator
			HashMap<Strand, SamLocusIterator> threadIterators = previousIterators.get();
			if (threadIterators == null) {
				threadIterators = new HashMap<>();
				previousIterators.set(threadIterators);

				HashMap<Strand, Interval> threadPreviousIntervals = new HashMap<>();
				previousIntervals.set(threadPreviousIntervals);
			}

			SamLocusIterator iterator = threadIterators.get(strand);
			Interval interval = null;
			if (iterator != null) {
				interval = previousIntervals.get().get(strand);
				if (!interval.getSequence().equals(contig) || interval.getEnd() < pos || interval.getStart() > pos) {

					//System.out.println("consuming iterator with interval "+interval+" want to go to: "+pos);
					while (iterator.hasNext()) {

						iterator.next(); //consume the iterator (it throws exception if we create another before
						// running out this one!)
					}
					currentLocus.remove(iterator);
					iterator = null;

				}
			}
			if (iterator == null) {

				int ADAPTED_WINDOW = INTERVAL_WINDOW;
				while (iterator == null && ADAPTED_WINDOW > 0) {

					IntervalList site = new IntervalList(reader.getFileHeader());
					interval = new Interval(contig, pos, Math.min(reader.getFileHeader().getSequence(contig)
							.getSequenceLength(), pos + ADAPTED_WINDOW)); //don't be cowboy here, using 1M gives
					// BufferUnderflow
					site.add(interval);

					iterator = createSamLocusIterator(reader, site, filters);

					try {
						iterator.iterator();//call to avoid NPE
					} catch (BufferUnderflowException e) {
						iterator = null;

						ADAPTED_WINDOW /= 2;
						System.err.println("Reducing window due to BufferUnderflow to " + ADAPTED_WINDOW);
					}
				}
				if (iterator == null) {
					throw new RuntimeException("Cant create iterator due to buffer underflow exceptions");
				}
			}
			//put the iterator in cache
			try {
				previousIterators.get().put(strand, iterator);
				previousIntervals.get().put(strand, interval);
			} catch (NullPointerException e) {
				System.out.println("null pointer");
			}
			LocusInfo li = null;

			//check if the previous iteration has reached exactly the position we are aiming at
			if (currentLocus.get(iterator) != null && currentLocus.get(iterator).getPosition() == pos) {
				li = currentLocus.get(iterator);
			} else {
				while (iterator.hasNext()) {
					li = iterator.next();

					if (li.getPosition() >= pos) {
						break;
					}
				}
				currentLocus.put(iterator, li);
			}

			if (li != null && li.getPosition() == pos) {
				//System.out.println("yes");
				Collection<RecordAndOffset> ros = li.getRecordAndPositions();

				if (filterDuplicates) {
					HashMap<Integer, RecordAndOffset> uniqueReads = new HashMap<Integer, RecordAndOffset>();
					for (RecordAndOffset read : ros) {
						//System.err.println("is clonal: "+read.getRecord().getReadName());
						if (strand.isNegative()) {
							if (!uniqueReads.containsKey(read.getRecord().getAlignmentEnd()) ||
									uniqueReads.get(read.getRecord().getAlignmentEnd()).getRecord().getMappingQuality
											() < read.getRecord().getMappingQuality()) {
								uniqueReads.put(read.getRecord().getAlignmentEnd(), read);
							}/*else{
								System.err.println("filtering duplicate: "+read.getRecord().getReadName());
							}*/
						} else {
							if (!uniqueReads.containsKey(read.getRecord().getAlignmentStart()) ||
									uniqueReads.get(read.getRecord().getAlignmentStart()).getRecord()
											.getMappingQuality() < read.getRecord().getMappingQuality()) {
								uniqueReads.put(read.getRecord().getAlignmentStart(), read);
							}/*else{
								System.err.println("filtering duplicate: "+read.getRecord().getReadName());
							}*/
						}
					}
					ros = new LinkedList<RecordAndOffset>();
					ros = uniqueReads.values();
				}

				for (RecordAndOffset ro : ros) {
					if (ro.getReadBase() != trimbase) {
						toret.add(ro);
					}
				}

			}
		}

		return toret;
	}

	private static SamLocusIterator createSamLocusIterator(
			SAMFileReader reader,
			IntervalList site, List<SamRecordFilter> filters) {
		SamLocusIterator iterator;
		iterator = new SamLocusIterator(reader, site, true);


		iterator.setSamFilters(filters);

		iterator.setEmitUncoveredLoci(false);
		return iterator;
	}

} // fin de clase
