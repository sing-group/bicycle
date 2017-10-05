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

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.PileupElementFilter;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.TextCigarCodec;

public class ListerFilter extends ReadFilter {

	@Argument(doc = "control genome for error computation", required = false)
	public int trimUntil = 4;

	public static boolean trim = false;

	@Argument(doc = "remove bad bisulfited", required = false)
	public boolean removeBad = false;

	@Argument(doc = "remove ambigous read using tagged reads with ZA flag", required = false)
	public boolean removeAmbiguous = false;

	@Argument(doc = "keep only reads with more than one alignment", required = false)
	public boolean onlyWithOneAlignment = false;

	public static char TRIMMED_BASE = 'X';

	public Pattern watson = Pattern.compile(".*[Cc][^Gg].*[Cc][^Gg].*[Cc][^Gg].*[Cc][^Gg].*");
	public Pattern crick = Pattern.compile(".*[^Cc][Gg].*[^Cc][Gg].*[^Cc][Gg].*[^Cc][Gg].*");

	private GenomeAnalysisEngine engine;

	@Argument(doc = "ignore bases with less depth of coverage")
	public static int mindepth = 1;


	// counters
	private long unmmapedReadsCounter = 0;
	private long withMoreThanOneAlignmentCounter = 0;
	private long ambiguousReadCounter = 0;
	private long trimmedCounter = 0;
	private long badBisulfitedCounter = 0;
	private long processedReadsCounter = 0;

	private ThreadLocal<Boolean> freezeCountersInThread = new ThreadLocal<>();

	public void freezeCountersInThread() {
		this.freezeCountersInThread.set(true);
	}
	public void unfreezeCountersInThread() {
		this.freezeCountersInThread.set(false);
	}

	public void resetCounters() {
		this.processedReadsCounter = 0;
		this.unmmapedReadsCounter = 0;
		this.withMoreThanOneAlignmentCounter = 0;
		this.ambiguousReadCounter = 0;
		this.trimmedCounter = 0;
		this.badBisulfitedCounter = 0;
	}

	public long getUnmmapedReadsCounter() {
		return unmmapedReadsCounter;
	}

	public long getWithMoreThanOneAlignmentCounter() {
		return withMoreThanOneAlignmentCounter;
	}

	public long getAmbiguousReadCounter() {
		return ambiguousReadCounter;
	}

	public long getBadBisulfitedCounter() {
		return badBisulfitedCounter;
	}

	public long getProcessedReadsCounter() {
		return processedReadsCounter;
	}

	public long getTrimmedCounter() {
		return trimmedCounter;
	}

	@Override
	public void initialize(GenomeAnalysisEngine engine) {
		this.engine = engine;
	}

	private Set<List<StackTraceElement>> stackTraces = new HashSet<>();
	@Override
	public boolean filterOut(SAMRecord record) {

		if (this.freezeCountersInThread.get() == null) {
			this.freezeCountersInThread.set(false);
		}

		//remove reads with two or more alignments

		if (this.freezeCountersInThread.get() == false) this.processedReadsCounter ++;
		if (record.getReadUnmappedFlag()) {

			if (this.freezeCountersInThread.get() == false) this.unmmapedReadsCounter++;
			System.err.println("UNMAPPED READ!");
			System.exit(1);
			return false;
		}

		if (onlyWithOneAlignment) {
			// bowtie 1
			if (record.getHeader().getProgramRecord("Bowtie") != null) {
				if (record.getAttribute("XM") != null
						&& ((Integer) record.getAttribute("XM")) > 1) {
					if (this.freezeCountersInThread.get() == false) this.withMoreThanOneAlignmentCounter++;
					return true;
				}
			}
			// bowtie 2
			if (record.getHeader().getProgramRecord("bowtie2") != null) {
				if (record.getAttribute("XS") != null) {
					if (this.freezeCountersInThread.get() == false) this.withMoreThanOneAlignmentCounter++;
					return true;
				}
			}
		}

		//remove ambiguous
		if (removeAmbiguous) {
			if (record.getAttribute("ZA") != null && record.getAttribute("ZA").toString().equals("Y")) {
				if (this.freezeCountersInThread.get() == false) this.ambiguousReadCounter++;
				return true;
			}
		}

		//trim to x mismatch (trim must be before bad bisulfited)
		if (trim) {
			boolean trimmed = trim(record, trimUntil);
			if (trimmed && this.freezeCountersInThread.get() == false) this.trimmedCounter++;
		}
		//bad bisulfited filter
		if (removeBad) {
			if (record.getReadGroup().getId().equals("WATSON") && watson.matcher(record.getReadString()).find()
					||
					record.getReadGroup().getId().equals("CRICK") && crick.matcher(record.getReadString()).find()) {
				if (this.freezeCountersInThread.get() == false) this.badBisulfitedCounter++;
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

	private boolean trim(SAMRecord record, int trimMismatches) {
		record.setAttribute("XT", "true");
		String sequenceCT = record.getReadString();

		int mismatches = (Integer) record.getAttribute("NM");

		if (mismatches >= trimMismatches) {

			boolean isReverse = record.getReadNegativeStrandFlag();


			int currentMismatches = 0;
			int trimPoint = isReverse ? sequenceCT.length() : 0;

			CigarIterator cigarIterator =
					new CigarIterator(record.getCigarString(), isReverse);

			MDTagIterator mdIterator = new MDTagIterator((String) record.getAttribute("MD"), isReverse);

			CigarOperator lastCigarOperator = null;
			while (cigarIterator.hasNext() && currentMismatches < trimMismatches) {
				CigarOperator operator = cigarIterator.next();

				switch (operator) {
					case M:
						MDTagOperator mdOperator = mdIterator.next();
						if (mdOperator == MDTagOperator.VARIANT) {
							currentMismatches++;
							if (currentMismatches < trimMismatches) {
								trimPoint = trimPoint + (isReverse ? -1 : 1);
							}
						} else if (mdOperator == MDTagOperator.SEQUENCE_MATCH) {
							// do nothing
							trimPoint = trimPoint + (isReverse ? -1 : 1);
						} else {
							//should not happen
							return false; // do not trim
						}

						break;
					case I:
						currentMismatches++;
						if (currentMismatches < trimMismatches) {
							trimPoint = trimPoint + (isReverse ? -1 : 1);
						}
						break;
					case D:
						currentMismatches++;
						mdIterator.next();
						break;
					default:
						return false; //do not trim this read, since it has non-supported CIGAR operations
				}
				lastCigarOperator = operator;
			}

			if (currentMismatches >= trimMismatches) {
				if (isReverse) {
					//trimPoint ++; //trim the last mismatch
					record.setReadString(record.getReadString().substring(0, trimPoint).replaceAll(".", "" +
							TRIMMED_BASE)
							+ record.getReadString().substring(trimPoint));
				} else {
					//trimPoint --; //trim the last mismatch
					record.setReadString(record.getReadString().substring(0, trimPoint) + record.getReadString()
							.substring
									(trimPoint).replaceAll(".", "" + TRIMMED_BASE));
				}
				return true;
			}

		}
		return false;

	}

	@Override
	public String toString() {

		double ambiguousRatio = (double)this.getAmbiguousReadCounter()/(double)this.getProcessedReadsCounter();
		double onlyWithOneAlignmentRatio = (double)this.getWithMoreThanOneAlignmentCounter()/(double)this.getProcessedReadsCounter();
		double badBisulfitedRatio = (double)this.getBadBisulfitedCounter()/(double)this
				.getProcessedReadsCounter();
		double trimmedRatio = (double)this.getTrimmedCounter()/(double)this.getProcessedReadsCounter();

		return "Mapped reads processed: "+this.getProcessedReadsCounter()+"" +
				", remove ambiguous reads: " + this.removeAmbiguous +
					(this.removeAmbiguous ? " (" + this.getAmbiguousReadCounter() + " removed " +
						"("+asPercent(ambiguousRatio)+"))"
						: "") +
				", remove with more than one alignment: " + this.onlyWithOneAlignment +
					(this.onlyWithOneAlignment ? " (" + this.getWithMoreThanOneAlignmentCounter() + " removed " +
						"("+asPercent(onlyWithOneAlignmentRatio)+"))"
						: "") +
				", remove non-correctly bisulfite-converted reads: " + this.removeBad +
					(this.removeBad ? " (" + this.getBadBisulfitedCounter() + " removed " +
							"("+asPercent(badBisulfitedRatio)+"))"
						: "") +
				", trim to 'x' mismatch: " + this.trim +
					(this.trim ? " x=" + this.trimUntil + " " + this.getTrimmedCounter() + " trimmed " +
						"("+asPercent(trimmedRatio)+")"
						:"");
	}

	private String asPercent(double decimal) {
		DecimalFormat fmt = new DecimalFormat("0.00");

		return fmt.format(decimal * 100d) + "%";
	}

	protected static class CigarIterator implements Iterator<CigarOperator> {

		private final boolean reverse;
		private List<CigarElement> cigarElements;
		private int currentCigarElementsPos;

		private int repetitionsOfCurrentOperation = 0;
		private CigarOperator currentOperator = null;

		public CigarIterator(String cigar, boolean reverse) {
			this.cigarElements = new ArrayList(TextCigarCodec.getSingleton().decode(cigar).getCigarElements());
			this.reverse = reverse;
			this.currentCigarElementsPos = this.reverse ? cigarElements.size() : -1;
			this.advanceCigarElement();
		}

		@Override
		public boolean hasNext() {
			return repetitionsOfCurrentOperation > 0
					|| (this.currentCigarElementsPos >= 0 && this.currentCigarElementsPos < this.cigarElements.size());
		}

		@Override
		public CigarOperator next() {

			CigarOperator result = this.currentOperator;
			this.repetitionsOfCurrentOperation--;
			if (repetitionsOfCurrentOperation == 0 && (this.currentCigarElementsPos >= 0 && this
					.currentCigarElementsPos < this.cigarElements.size())) {
				advanceCigarElement();
			}

			return result;
		}

		private void advanceCigarElement() {
			this.currentCigarElementsPos = this.currentCigarElementsPos + (this.reverse ? -1 : 1);
			if (this.currentCigarElementsPos >= 0 && this.currentCigarElementsPos < this.cigarElements.size()) {
				CigarElement element = this.cigarElements.get(this.currentCigarElementsPos);
				this.repetitionsOfCurrentOperation = element.getLength();
				this.currentOperator = element.getOperator();
			}
		}
	}

	protected enum MDTagOperator {
		DELETION, SEQUENCE_MATCH, VARIANT
	}

	protected static abstract class MDTagOperation {
		public abstract int getLength();

		public abstract MDTagOperator getOperator();
	}

	protected static class SequenceMatchOperation extends MDTagOperation {
		int length;

		public SequenceMatchOperation(int length) {
			this.length = length;
		}

		@Override
		public int getLength() {
			return length;
		}

		@Override
		public MDTagOperator getOperator() {
			return MDTagOperator.SEQUENCE_MATCH;
		}
	}

	protected static class SequenceVariantOperation extends MDTagOperation {
		private char allele;

		public SequenceVariantOperation(char alelle) {
			this.allele = allele;
		}

		@Override
		public int getLength() {
			return 1;
		}

		@Override
		public MDTagOperator getOperator() {
			return MDTagOperator.VARIANT;
		}
	}

	protected static class DeletionOperation extends MDTagOperation {
		private String deletedBases;

		public DeletionOperation(String deletedBases) {
			this.deletedBases = deletedBases;
		}

		public String getDeletedBases() {
			return deletedBases;
		}

		@Override
		public int getLength() {
			return deletedBases.length();
		}

		@Override
		public MDTagOperator getOperator() {
			return MDTagOperator.DELETION;
		}
	}

	protected static class MDTagIterator implements Iterator<MDTagOperator> {
		static final Pattern mdPat = Pattern.compile("\\G(?:([0-9]+)|([ACTGNactgn])|(\\^[ACTGNactgn]+))");

		private final boolean reverse;
		private List<MDTagOperation> mdOperations;
		private int currentMDOperationsPos;

		private int repetitionsOfCurrentOperation = 0;
		private MDTagOperator currentOperator = null;

		private List<MDTagOperation> decodeMDTag(String mdString) {
			Matcher matcher = mdPat.matcher(mdString);

			ArrayList result = new ArrayList<>();
			while (matcher.find()) {
				if (matcher.group(1) != null) {
					// match
					int length = Integer.parseInt(matcher.group(1));
					if (length > 0) {
						result.add(new SequenceMatchOperation(Integer.parseInt(matcher.group(1))));
					}
				} else if (matcher.group(2) != null) {
					// variant
					result.add(new SequenceVariantOperation(matcher.group(2).charAt(0)));
				} else if (matcher.group(3) != null) {
					// deletion
					result.add(new DeletionOperation(matcher.group(3).substring(1)));
				}
			}

			return result;
		}

		public MDTagIterator(String mdString, boolean reverse) {
			this.mdOperations = decodeMDTag(mdString);
			this.reverse = reverse;
			this.currentMDOperationsPos = this.reverse ? mdOperations.size() : -1;
			this.advanceMDOperation();
		}

		@Override
		public boolean hasNext() {
			return repetitionsOfCurrentOperation > 0
					|| (this.currentMDOperationsPos >= 0 && this.currentMDOperationsPos < this.mdOperations.size());
		}

		@Override
		public MDTagOperator next() {

			MDTagOperator result = this.currentOperator;
			this.repetitionsOfCurrentOperation--;
			if (repetitionsOfCurrentOperation == 0 && (this.currentMDOperationsPos >= 0 && this
					.currentMDOperationsPos < this.mdOperations.size())) {
				advanceMDOperation();
			}

			return result;
		}

		private void advanceMDOperation() {
			this.currentMDOperationsPos = this.currentMDOperationsPos + (this.reverse ? -1 : 1);
			if (this.currentMDOperationsPos >= 0 && this.currentMDOperationsPos < this.mdOperations.size()) {
				MDTagOperation element = this.mdOperations.get(this.currentMDOperationsPos);
				this.repetitionsOfCurrentOperation = element.getLength();
				this.currentOperator = element.getOperator();
			}
		}
	}

}
