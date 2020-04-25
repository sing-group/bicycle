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

import static es.cnio.bioinfo.bicycle.FastqSplitter.splitfastq;
import static java.io.File.createTempFile;
import static java.util.Arrays.asList;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

import es.cnio.bioinfo.bicycle.FileSequenceInputStream;

public class FasqSplitterTest {

  private static String FASTQ_1 = "@sequence1name\nsequence1\n+aaa\nquality1\n@sequence2name\nsequence2\n+\nquality2";
  private static String FASTQ_2 = "@sequence3name\nsequence3\n+\n@quality3";
  private static String FASTQ_3 = "@sequence4name\nsequence4\n+\n@quality4\n@sequence5name\nsequence5\n+\nquality5";
  
  private static File FASTQ_FILE_1;
  private static File FASTQ_FILE_2;
  private static File FASTQ_FILE_3;

  @BeforeClass
  public static void prepareFiles() throws IOException {
    FASTQ_FILE_1 = createTempFile("fastqexample", ".fastq");
    FASTQ_FILE_1.deleteOnExit();
    PrintStream ps1 = new PrintStream(new FileOutputStream(FASTQ_FILE_1));

    FASTQ_FILE_2 = createTempFile("fastqexample", ".fastq");
    FASTQ_FILE_2.deleteOnExit();
    PrintStream ps2 = new PrintStream(new FileOutputStream(FASTQ_FILE_2));

    FASTQ_FILE_3 = createTempFile("fastqexample", ".fastq");
    FASTQ_FILE_3.deleteOnExit();
    PrintStream ps3 = new PrintStream(new FileOutputStream(FASTQ_FILE_3));

    ps1.println(FASTQ_1);
    ps2.println(FASTQ_2);
    ps3.println(FASTQ_3);

    ps1.close();
    ps2.close();
    ps3.close();
  }
  
  @Test
  public void testBasicSplitting() throws IOException {
    List<BufferedReader> iss = splitfastq(asList(new File[]{FASTQ_FILE_1, FASTQ_FILE_2, FASTQ_FILE_3}), 2);
    assertEquals(2, iss.size());
    assertSplitIsAdjusted(iss);
    
    iss = splitfastq(asList(new File[]{FASTQ_FILE_1, FASTQ_FILE_2, FASTQ_FILE_3}), 20);
    assertEquals(20, iss.size());
    assertSplitIsAdjusted(iss);

    iss = splitfastq(Arrays.asList(new File[]{FASTQ_FILE_1, FASTQ_FILE_2, FASTQ_FILE_3}), 1);
    assertEquals(1, iss.size());
    assertSplitIsAdjusted(iss);
  }

  @Test
  public void testFileSequenceInputStream() throws Exception {
    int maxRead = 125;
    InputStream is = new FileSequenceInputStream(asList(new File[]{FASTQ_FILE_1, FASTQ_FILE_2, FASTQ_FILE_3}), maxRead, 0);
    try (BufferedReader reader = new BufferedReader(new InputStreamReader(is))) {
    
    String line = null;
    int readedChars = 0;
    while ((line = reader.readLine()) != null) {
      readedChars += line.length() + 1;
    }
    assertEquals(maxRead, readedChars - 1);
    }
  }
  private void assertSplitIsAdjusted(List<BufferedReader> iss) throws IOException {
    for (BufferedReader is : iss) {
      String line = null;
      int lineCount = 0;
      while ((line = is.readLine()) != null) {
        if (lineCount % 4 == 0) {
          assertTrue(line.startsWith("@"));
        }
        lineCount ++;
        
      }
      assertTrue(lineCount % 4 == 0);
    }
  }
}
