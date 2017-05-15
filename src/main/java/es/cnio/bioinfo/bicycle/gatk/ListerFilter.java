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

package es.cnio.bioinfo.bicycle.gatk;

import java.util.LinkedList;
import java.util.List;
import java.util.regex.Pattern;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.PileupElementFilter;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import net.sf.samtools.SAMRecord;

public class ListerFilter extends ReadFilter {

	@Argument(doc = "control genome for error computation", required = false)
	public int trimUntil = 4;

	public static boolean trim = false;

	@Argument(doc = "remove bad bisulfited", required = false)
	public boolean removeBad = false;

	@Argument(doc = "remove ambigous read using tagged reads with ZA flag", required = false)
	public boolean removeAmbiguous = false;

	public static char TRIMMED_BASE = 'X';

	public Pattern watson = Pattern.compile(".*[Cc][^Gg].*[Cc][^Gg].*[Cc][^Gg].*[Cc][^Gg].*");
	public Pattern crick = Pattern.compile(".*[^Cc][Gg].*[^Cc][Gg].*[^Cc][Gg].*[^Cc][Gg].*");

	private GenomeAnalysisEngine engine;

	@Argument(doc = "ignore bases with less depth of coverage")
	public static int mindepth = 1;


	@Override
	public void initialize(GenomeAnalysisEngine engine) {
		this.engine = engine;
	}

	@Override
	public boolean filterOut(SAMRecord record) {
		//remove reads with two or more alignments
		if (record.getAttribute("XM") != null && !record.getAttribute("XM").toString().equals("1")) {
			//System.err.println("filtering record with two possible alignments");
			return true;
		}

		//remove ambiguous
		if (removeAmbiguous) {
			//if (record.getAttribute("ZA")!=null)System.out.println("attribute: "+record.getAttribute("ZA").toString
			// ());
			if (record.getAttribute("ZA") != null && record.getAttribute("ZA").toString().equals("Y")) {
				return true;
			}
		}

		//trim to x mismatch (trim must be before bad bisulfited)
		if (trim) {
			trim(record, trimUntil);
		}
		//bad bisulfited filter
		if (removeBad) {
			if (!record.getReadNegativeStrandFlag() && watson.matcher(record.getReadString()).find()
					||
					record.getReadNegativeStrandFlag() && crick.matcher(record.getReadString()).find()) {
				return true;
			}
		}
		return false;
	}


	public static ReadBackedPileup applyFilters(ReadBackedPileup pileup) {
		if (pileup == null) {
			return null;
		}
		ReadBackedPileup toret = pileup.getFilteredPileup(new PileupElementFilter() {

			@Override
			public boolean allow(PileupElement arg0) {
				for (PileupElementFilter filter : pileupfilters) {
					if (!filter.allow(arg0)) {
						return false;
					}
				}
				return true;
			}

		});

		if (toret.getBases().length < mindepth) {
			return null;
		}
		return toret;

	}

	static List<PileupElementFilter> pileupfilters;

	static {

		pileupfilters = new LinkedList<PileupElementFilter>();
		pileupfilters.add(new PileupElementFilter() {

			@Override
			public boolean allow(PileupElement arg0) {

				if (trim) {

					if (arg0.getBase() == ListerFilter.TRIMMED_BASE) {
						return false;
					}
				}
				return true;
			}

		});

	}

	private void trim(SAMRecord record, int trimMismatches) {
		record.setAttribute("XT", "true");
		String sequenceCT = record.getReadString();

		int mismatches = (Integer) record.getAttribute("NM");
		if (mismatches >= trimMismatches) {
			String MDString = (String) record.getAttribute("MD");
			String[] MDStringTokens = MDString.split("[ACTGN]");

			if (record.getReadNegativeStrandFlag()) {
				//is in the reverse
				int newStart = sequenceCT.length();
				for (int m = MDStringTokens.length - 1; m >= MDStringTokens.length - trimMismatches; m--) {
					newStart -= Integer.parseInt(MDStringTokens[m]) + 1;
				}
				newStart++;
				sequenceCT = sequenceCT.substring(newStart);


				record.setReadString(record.getReadString().substring(0, newStart).replaceAll(".", "" + TRIMMED_BASE)
						+ record.getReadString().substring(newStart));
			} else {
				int newLength = 0;

				for (int m = 0; m < trimMismatches; m++) {
					newLength += Integer.parseInt(MDStringTokens[m]) + 1;
				}
				newLength--; //ignore the lastmismatch				
				sequenceCT = sequenceCT.substring(0, newLength);

				record.setReadString(record.getReadString().substring(0, newLength) + record.getReadString().substring
						(newLength).replaceAll(".", "" + TRIMMED_BASE));
			}
		}

	}

	@Override
	public String toString() {
		return "remove ambiguous reads: " + this.removeAmbiguous + ", remove non-correctly bisulfite-converted reads: " +
				"" + this.removeBad + ", trim to 'x' mismatch: " + this.trim + ((this.trim) ? " x=" + this.trimUntil : "");
	}

}
