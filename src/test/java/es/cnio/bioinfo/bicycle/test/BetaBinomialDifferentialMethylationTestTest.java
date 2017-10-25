package es.cnio.bioinfo.bicycle.test;

import org.junit.Assert;
import org.junit.Test;

import es.cnio.bioinfo.bicycle.operations.BetaBinomialDifferentialMethylationTest;

/**
 * This test case ensures that bicycle produces the same results as MethylSig for this example:
 * Input files:
 * ==> case1.txt <==
 * chrBase chr     base    strand  coverage        freqC   freqT
 * chr21.9764539   chr21   9764539 R       6      50   50
 * chr21.9764540   chr21   9764540 R       6      50   50
 *
 * ==> case2.txt <==
 * chrBase chr     base    strand  coverage        freqC   freqT
 * chr21.9764539   chr21   9764539 R       8      25.00   75.00
 * chr21.9764540   chr21   9764540 R       8      25.00   75.00
 *
 * ==> case3.txt <==
 * chrBase chr     base    strand  coverage        freqC   freqT
 * chr21.9764539   chr21   9764539 R       8      12.5   87.5
 * chr21.9764540   chr21   9764540 R       8      12.5   87.5
 *
 * ==> control1.txt <==
 * chrBase chr     base    strand  coverage        freqC   freqT
 * chr21.9764540   chr21   9764540 R       8      87.5   12.5

 * ==> control2.txt <==
 * chrBase chr     base    strand  coverage        freqC   freqT
 * chr21.9764539   chr21   9764539 R       9      100   0
 * chr21.9764540   chr21   9764540 R       9      100   0
 *
 * ==> control3.txt <==
 * chrBase chr     base    strand  coverage        freqC   freqT
 * chr21.9764539   chr21   9764539 R       3      100   0
 * chr21.9764540   chr21   9764540 R       3      100   0
 *
 * ==> sites.bed <==
 * chr21	9764538	9764540	promoter1	1	.
 * chr10	10	20	promoter2	1	.
 *
 *
 * In order to get p-values with methylSig for these two sites, you have to run
 *  fileList <- c("case1.txt", "case2.txt", "case3.txt", "control1.txt", "control2.txt", "control3.txt");
 *  sample.id <- c("CASE1", "CASE2", "CASE3", "CONTROL1", "CONTROL2", "CONTROL3");
 *  treatment <- c(1,1,1,0,0,0)
 *  meth <- methylSigReadData(fileList, sample.ids = sample.id, treatment = treatment, destranded=FALSE, minCount = 0);
 *  methylSigCalc(meth, min.per.group=c(2,2));
 *
 *  # for regions
 *  tfbsInfo <- getTFBSInfo("sites.bed")
 *  methPromoters <- methylSigTileTFBS(meth, tfbsInfo);
 *  methylSigCalc(methPromoters, min.per.group=c(2,2));
 */
public class BetaBinomialDifferentialMethylationTestTest {

	/**
	 * Tests that the p-value for the cytosine chr21.9764540 is the one given by MethylSig
	 */
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


	/**
	 * Tests that the p-value for the cytosine chr21.9764539 is the one given by MethylSig
	 */
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


	/**
	 * Tests that the p-value for the region "promoter1" is the one given by MethylSig
	 */
	@Test
	public void testRegion() {
		int[] treatmentCytosines = new int[]{6, 4, 2};
		int[] treatmentDepth = new int[]{12, 16, 16};
		int[] controlCytosines = new int[]{7, 18, 6};
		int[] controlDepth = new int[]{8, 18, 6};

		BetaBinomialDifferentialMethylationTest test = new BetaBinomialDifferentialMethylationTest();

		double pVal = test.getPvalue(treatmentCytosines, treatmentDepth, controlCytosines, controlDepth);
		//System.err.println(pVal);
		Assert.assertEquals(0.001166, pVal, 1e-5);

	}
}
