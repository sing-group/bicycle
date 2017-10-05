package es.cnio.bioinfo.bicycle.gatk;


import static junit.framework.Assert.assertEquals;

import java.io.IOException;

import org.junit.Test;

import es.cnio.bioinfo.bicycle.gatk.ListerFilter.CigarIterator;
import net.sf.samtools.CigarOperator;

public class CigarIteratorTest {

	@Test
	public void testPerfectMatch() {
		CigarIterator iterator = new CigarIterator("3M", false);

		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
	}

	@Test
	public void testPerfectMatchReverse() {
		CigarIterator iterator = new CigarIterator("3M", true);

		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
	}

	@Test
	public void testInsertions() {
		CigarIterator iterator = new CigarIterator("3M2I1M", false);

		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
		assertEquals(CigarOperator.INSERTION, iterator.next());
		assertEquals(CigarOperator.INSERTION, iterator.next());
		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
	}

	@Test
	public void testInsertionsReverse() {
		CigarIterator iterator = new CigarIterator("3M2I1M", true);

		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
		assertEquals(CigarOperator.INSERTION, iterator.next());
		assertEquals(CigarOperator.INSERTION, iterator.next());
		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
	}

	@Test
	public void testDeletions() {
		CigarIterator iterator = new CigarIterator("3M2D1M", false);

		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
		assertEquals(CigarOperator.DELETION, iterator.next());
		assertEquals(CigarOperator.DELETION, iterator.next());
		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
	}

	@Test
	public void testDeletionsReverse() {
		CigarIterator iterator = new CigarIterator("3M2D1M", true);

		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
		assertEquals(CigarOperator.DELETION, iterator.next());
		assertEquals(CigarOperator.DELETION, iterator.next());
		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
		assertEquals(CigarOperator.MATCH_OR_MISMATCH, iterator.next());
	}

}
