package es.cnio.bioinfo.bicycle.gatk;


import static junit.framework.Assert.assertEquals;

import org.junit.Test;

import es.cnio.bioinfo.bicycle.gatk.ListerFilter.MDTagIterator;
import es.cnio.bioinfo.bicycle.gatk.ListerFilter.MDTagOperator;
import net.sf.samtools.CigarOperator;

public class MDIteratorTest {

	@Test
	public void testPerfectMatch() {
		MDTagIterator iterator = new MDTagIterator("3", false);

		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());

	}

	@Test
	public void testPerfectMatchReverse() {
		MDTagIterator iterator = new MDTagIterator("3", true);

		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
	}

	@Test
	public void testVariants() {
		MDTagIterator iterator = new MDTagIterator("3A1", false);

		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.VARIANT, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
	}


	@Test
	public void testVariantsReverse() {
		MDTagIterator iterator = new MDTagIterator("3A1", true);

		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.VARIANT, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
	}


	@Test
	public void testConsecutiveVariants() {
		MDTagIterator iterator = new MDTagIterator("3A0C1", false);

		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.VARIANT, iterator.next());
		assertEquals(MDTagOperator.VARIANT, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
	}

	@Test
	public void testDeletions() {
		MDTagIterator iterator = new MDTagIterator("1^ACG0T2", false);

		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.DELETION, iterator.next());
		assertEquals(MDTagOperator.DELETION, iterator.next());
		assertEquals(MDTagOperator.DELETION, iterator.next());
		assertEquals(MDTagOperator.VARIANT, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
	}

	@Test
	public void testDeletionsReverse() {
		MDTagIterator iterator = new MDTagIterator("1^ACG0T2", true);

		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
		assertEquals(MDTagOperator.VARIANT, iterator.next());
		assertEquals(MDTagOperator.DELETION, iterator.next());
		assertEquals(MDTagOperator.DELETION, iterator.next());
		assertEquals(MDTagOperator.DELETION, iterator.next());
		assertEquals(MDTagOperator.SEQUENCE_MATCH, iterator.next());
	}
}
