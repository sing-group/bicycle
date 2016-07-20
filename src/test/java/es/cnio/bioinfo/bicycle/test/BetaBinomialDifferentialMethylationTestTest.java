package es.cnio.bioinfo.bicycle.test;

import org.junit.Assert;
import org.junit.Test;

import es.cnio.bioinfo.bicycle.operations.BetaBinomialDifferentialMethylationTest;

public class BetaBinomialDifferentialMethylationTestTest {

	@Test
	public void testPValue() {
		 int[] treatmentCytosines = new int[]{3, 2, 1};
		 int[] treatmentDepth = new int[]{6, 8, 8};
		 int[] controlCytosines = new int[]{7, 9, 3};
		 int[] controlDepth = new int[]{8, 9, 3};
		

		BetaBinomialDifferentialMethylationTest test = new BetaBinomialDifferentialMethylationTest();

		double pVal = test.getPvalue(treatmentCytosines, treatmentDepth, controlCytosines, controlDepth);
		//System.err.println(pVal);
		Assert.assertEquals(0.00302370, pVal, 1e-5);
		
	}
	
	@Test
	public void testNoDepth() {
		 int[] treatmentCytosines = new int[]{3, 2, 1};
		 int[] treatmentDepth = new int[]{6, 8, 8};
		 int[] controlCytosines = new int[]{0, 9, 3};
		 int[] controlDepth = new int[]{0, 9, 3};
		

		BetaBinomialDifferentialMethylationTest test = new BetaBinomialDifferentialMethylationTest();
		
		double pVal = test.getPvalue(treatmentCytosines, treatmentDepth, controlCytosines, controlDepth);
		//System.err.println(pVal);
		Assert.assertEquals(0.00579775, pVal, 1e-8);
		
	}
}
