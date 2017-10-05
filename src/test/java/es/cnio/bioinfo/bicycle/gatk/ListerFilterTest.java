package es.cnio.bioinfo.bicycle.gatk;

import static junit.framework.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import junit.framework.Assert;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMRecord;

public class ListerFilterTest {


	private ListerFilter filter = null;
	private SAMRecord aSAMRecord = null;

	@Test
	public void testTrim() {
		aSAMRecord = createSAMRecord("T", "1M", "A", 1, false);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("X", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("T", "1M", "1", 0, false);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("T", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTATTAT", "7M", "2A2A1", 2, false);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("TTXXXXX", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTATTAT", "7M", "2A2A1", 2, false);
		filter = createTrimmingFilterUntil(2);
		filter.filterOut(aSAMRecord);
		assertEquals("TTATTXX", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTATTAT", "7M", "2A2A1", 2, false);
		filter = createTrimmingFilterUntil(3);
		filter.filterOut(aSAMRecord);
		assertEquals("TTATTAT", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("ATATTAT", "7M", "A1A2A1", 3, false);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("XXXXXXX", aSAMRecord.getReadString());
	}

	@Test

	public void testTrimReverse() {
		aSAMRecord = createSAMRecord("T", "1M", "A", 1, true);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("X", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("T", "1M", "1", 0, true);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("T", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTATTAT", "7M", "2A2A1", 2, true);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("XXXXXXT", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTATTAT", "7M", "2A2A1", 2, true);
		filter = createTrimmingFilterUntil(2);
		filter.filterOut(aSAMRecord);
		assertEquals("XXXTTAT", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTATTAT", "7M", "2A2A1", 2, true);
		filter = createTrimmingFilterUntil(3);
		filter.filterOut(aSAMRecord);
		assertEquals("TTATTAT", aSAMRecord.getReadString());
	}

	@Test
	public void testTrimWithInsertions() {
		aSAMRecord = createSAMRecord("TTTTTTTT", "4M1I3M", "7", 2, false);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("TTTTXXXX", aSAMRecord.getReadString());


		aSAMRecord = createSAMRecord("TTTTTTTTT", "4M2I3M", "7", 2, false);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("TTTTXXXXX", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTTTTTTTT", "4M2I3M", "7", 2, false);
		filter = createTrimmingFilterUntil(2);
		filter.filterOut(aSAMRecord);
		assertEquals("TTTTTXXXX", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTTTTTATT", "4M2I3M", "4A2", 3, false);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("TTTTXXXXX", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTTTTTATT", "2M1I1M2I3M", "4A2", 3, false);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("TTXXXXXXX", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTTTTTATT", "2M1I1M2I3M", "4A2", 3, false);
		filter = createTrimmingFilterUntil(2);
		filter.filterOut(aSAMRecord);
		assertEquals("TTTTXXXXX", aSAMRecord.getReadString());
	}

	@Test
	public void testTrimWithInsertionsReverse() {
		aSAMRecord = createSAMRecord("TTTTTTTT", "4M1I3M", "7", 2, true);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("XXXXXTTT", aSAMRecord.getReadString());


		aSAMRecord = createSAMRecord("TTTTTTTTT", "4M2I3M", "7", 2, true);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("XXXXXXTTT", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTTTTTTTT", "4M2I3M", "7", 2, true);
		filter = createTrimmingFilterUntil(2);
		filter.filterOut(aSAMRecord);
		assertEquals("XXXXXTTTT", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTTTTTATT", "4M2I3M", "4A2", 3, true);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("XXXXXXXTT", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTTTTTATT", "2M1I1M2I3M", "4A2", 3, true);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("XXXXXXXTT", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTTTTTATT", "2M1I1M2I3M", "4A2", 3, true);
		filter = createTrimmingFilterUntil(2);
		filter.filterOut(aSAMRecord);
		assertEquals("XXXXXXATT", aSAMRecord.getReadString());
	}

	@Test
	public void testTrimWithDeletions() {
		aSAMRecord = createSAMRecord("TTTTTTT", "4M1D3M", "4^A3", 1, false);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("TTTTXXX", aSAMRecord.getReadString());


		aSAMRecord = createSAMRecord("TTTTTTT", "4M2D3M", "4^AA3", 2, false);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("TTTTXXX", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTTTTTT", "4M2D3M", "4^AA3", 2, false);
		filter = createTrimmingFilterUntil(2);
		filter.filterOut(aSAMRecord);
		assertEquals("TTTTXXX", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTTTATT", "4M2D3M", "4^AA1A2", 3, false);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("TTTTXXX", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTTATT", "2M1D1M2D3M", "2^A1^AA2", 3, false);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("TTXXXX", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTTTTT", "2M1D1M2D3M", "2^A1^AA3", 3, false);
		filter = createTrimmingFilterUntil(2);
		filter.filterOut(aSAMRecord);
		assertEquals("TTTXXX", aSAMRecord.getReadString());
	}

	@Test
	public void testTrimWithDeletionsReverse() {
		aSAMRecord = createSAMRecord("TTTTTTT", "4M1D3M", "4^A3", 1, true);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("XXXXTTT", aSAMRecord.getReadString());


		aSAMRecord = createSAMRecord("TTTTTTT", "4M2D3M", "4^AA3", 2, true);
		filter = createTrimmingFilterUntil(1);
		filter.filterOut(aSAMRecord);
		assertEquals("XXXXTTT", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTTTTTT", "4M2D3M", "4^AA3", 2, true);
		filter = createTrimmingFilterUntil(2);
		filter.filterOut(aSAMRecord);
		assertEquals("XXXXTTT", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TATTTTT", "4M2D3M", "1A2^AA3", 3, true);
		filter = createTrimmingFilterUntil(3);
		filter.filterOut(aSAMRecord);
		assertEquals("XXTTTTT", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTTTATT", "4M2D3M", "4^AA1A2", 3, true);
		filter = createTrimmingFilterUntil(2);
		filter.filterOut(aSAMRecord);
		assertEquals("XXXXATT", aSAMRecord.getReadString());

		aSAMRecord = createSAMRecord("TTTTTT", "3M2D1M1D2M", "3^AA1^A2", 3, true);
		filter = createTrimmingFilterUntil(2);
		filter.filterOut(aSAMRecord);
		assertEquals("XXXTTT", aSAMRecord.getReadString());
	}

	private SAMRecord createSAMRecord(String read, String cigar, String md, int NM, boolean reverse) {
		SAMFileHeader header = new SAMFileHeader();
		List<SAMProgramRecord> PRs = new ArrayList<>();
		PRs.add(new SAMProgramRecord("Bowtie"));
		header.setProgramRecords(PRs);

		SAMRecord record = new SAMRecord(header);
		record.setReadString(read);
		record.setAttribute("NM", NM); //number of mismatches
		record.setCigarString(cigar);
		record.setAttribute("MD", md);

		if (reverse) {
			record.setReadNegativeStrandFlag(true);
		}

		return record;
	}

	private ListerFilter createTrimmingFilterUntil(int mismatches) {
		ListerFilter filter = new ListerFilter();
		filter.trim = true;
		filter.trimUntil = mismatches;
		return filter;
	}
}
