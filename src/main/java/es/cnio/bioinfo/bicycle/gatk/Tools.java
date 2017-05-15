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

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.PileupElementFilter;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import net.sf.picard.filter.SamRecordFilter;
import net.sf.picard.util.SamLocusIterator.RecordAndOffset;
import net.sf.samtools.SAMFileReader;


public class Tools {


	///cache files (readers are per thread, in order to avoid synchronization)
	private HashMap<Thread, HashMap<Strand, ArrayList<SAMFileReader>>> perThreadReaders = new HashMap<Thread,
			HashMap<Strand, ArrayList<SAMFileReader>>>();


	public Collection<RecordAndOffset> getReadsForLocus(GenomeAnalysisEngine toolkit, Strand strand, String contig,
														int pos, boolean filterDuplicates) {
		HashMap<Strand, ArrayList<SAMFileReader>> perStrandReaders = perThreadReaders.get(Thread.currentThread());

		if (perStrandReaders == null) {

			perStrandReaders = new HashMap<Strand, ArrayList<SAMFileReader>>();
			perThreadReaders.put(Thread.currentThread(), perStrandReaders);
		}

		ArrayList<SAMFileReader> readers = perStrandReaders.get(strand);
		if (readers == null) {

			readers = new ArrayList<SAMFileReader>();
			for (File f : getFilesForStrand(toolkit, strand)) {
				SAMFileReader reader = new SAMFileReader(f, true);
				readers.add(reader);

			}
			perStrandReaders.put(strand, readers);
		}

		LinkedList<SamRecordFilter> filters = new LinkedList<SamRecordFilter>(toolkit.getFilters());

		return es.cnio.bioinfo.bicycle.Tools.getReadsForLocus(strand, contig, pos, filterDuplicates, readers,
				filters, ListerFilter.TRIMMED_BASE);
	}


	private HashMap<Strand, Collection<File>> strandFilesCache = new HashMap<Strand, Collection<File>>();

	private Collection<File> getFilesForStrand(GenomeAnalysisEngine toolkit, Strand strand) {
		Collection<File> toret = strandFilesCache.get(strand);
		if (toret == null) {
			toret = new LinkedList<File>();
			for (SAMReaderID readId : toolkit.getReadsDataSource().getReaderIDs()) {
				File file = toolkit.getReadsDataSource().getSAMFile(readId);
				if (strand.getFileRegExp().matcher(file.getName()).find()) {
					toret.add(file);
				}
			}
			strandFilesCache.put(strand, toret);
		}

		return toret;
	}


	private String getReadGroupsPerStrand(Strand strand, AlignmentContext arg2) {
		for (String group : arg2.getBasePileup().getReadGroups()) {
			if (strand.name().equals(group)) {

				return group;
			}
		}
		return null;
	}

	private int sumQuality(GATKSAMRecord read) {
		int sum = 0;
		for (byte b : read.getBaseQualities()) {
			sum += (int) b;
		}
		return sum;
	}

	public ReadBackedPileup getReadsForStrand(final Strand strand, final AlignmentContext alnContext, boolean
			filterDuplicates) {

		final String strandGroup = getReadGroupsPerStrand(strand, alnContext);

		if (strandGroup != null) {

			if (filterDuplicates) {

				ReadBackedPileup reads = alnContext.getBasePileup().getPileupForReadGroup(strandGroup);
				final HashMap<Integer, PileupElement> uniqueReads = new HashMap<Integer, PileupElement>();

				Iterator<PileupElement> it = reads.iterator();

				while (it.hasNext()) {
					PileupElement element = it.next();

					if (strand.isNegative()) {
						if (!uniqueReads.containsKey(element.getRead().getAlignmentEnd()) ||
								sumQuality(uniqueReads.get(element.getRead().getAlignmentEnd()).getRead()) <
										sumQuality(element.getRead())) {
							uniqueReads.put(element.getRead().getAlignmentEnd(), element);
						}
						/*else{
							System.err.println("REMOVING clonal: "+element.getRead());
						}*/

					} else {
						if (!uniqueReads.containsKey(element.getRead().getAlignmentStart()) ||
								sumQuality(uniqueReads.get(element.getRead().getAlignmentStart()).getRead()) <
										sumQuality(element.getRead())) {
							uniqueReads.put(element.getRead().getAlignmentStart(), element);
						}/*else{
							System.err.println("REMOVING clonal: "+element.getRead());
						}*/
					}
				}
				return alnContext.getBasePileup().getFilteredPileup(new PileupElementFilter() {

					@Override
					public boolean allow(PileupElement arg0) {
						//filter per read group
						if (!arg0.getRead().getReadGroup().getId().equals(strandGroup)) {
							return false;
						}
						if (strand.isNegative()) {
							PileupElement value = uniqueReads.get(arg0.getRead().getAlignmentEnd());
							if (value == arg0) {
								return true;
							} else {
								return false;
							}
						} else {
							PileupElement value = uniqueReads.get(arg0.getRead().getAlignmentStart());
							if (value == arg0) {
								return true;
							} else {
								return false;
							}
						}
					}

				});
			} else {


				return alnContext.getBasePileup().getFilteredPileup(new PileupElementFilter() {

					@Override
					public boolean allow(PileupElement arg0) {
						if (arg0.getRead().getReadGroup().getId().equals(strandGroup)) {
							return true;
						}
						return false;
					}

				});
			}
		} else {
			return null;
		}

	}

}
