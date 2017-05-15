package es.cnio.bioinfo.bicycle.operations;

import static es.uvigo.ei.sing.math.ArrayUtils.sumArray;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BracketingNthOrderBrentSolver;
import org.apache.commons.math3.analysis.solvers.UnivariateSolver;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.special.Gamma;

import es.uvigo.ei.sing.math.ArrayUtils;

/**
 * Computes differential methylation between 2 groups (treatment and control)
 * by using a beta-binomial model with a log-likelihood ratio test.
 * <p>
 * The method is the same as in MethylSig. Citation:
 * Park Y, Figueroa ME, Rozek LS, Sartor MA. MethylSig: a whole genome DNA
 * methylation analysis pipeline. Bioinformatics. 2014;30(17):2414-2422.
 * doi:10.1093/bioinformatics/btu339.
 *
 * @author lipido
 */
public class BetaBinomialDifferentialMethylationTest {


	public double getPvalue(int[] treatmentCytosines, int[] treatmentDepth,
							int[] controlCytosines, int[] controlDepth) {

		List<Integer> validControl = new ArrayList<>();
		List<Integer> validTreatment = new ArrayList<>();

		for (int i = 0; i < treatmentDepth.length; i++) {
			if (treatmentDepth[i] > 0) {
				validTreatment.add(i);
			}
		}

		for (int i = 0; i < controlDepth.length; i++) {
			if (controlDepth[i] > 0) {
				validControl.add(i);
			}
		}

		int[] validTreatmentCytosines = new int[validTreatment.size()];
		int[] validTreatmentDepth = new int[validTreatment.size()];
		int[] validControlCytosines = new int[validControl.size()];
		int[] validControlDepth = new int[validControl.size()];

		int i = 0;
		for (Integer validIndex : validTreatment) {
			validTreatmentCytosines[i] = treatmentCytosines[validIndex];
			validTreatmentDepth[i] = treatmentDepth[validIndex];
			i++;
		}
		i = 0;
		for (Integer validIndex : validControl) {
			validControlCytosines[i] = controlCytosines[validIndex];
			validControlDepth[i] = controlDepth[validIndex];
			i++;
		}

		return getPvalueValidSamples(validTreatmentCytosines, validTreatmentDepth, validControlCytosines,
				validControlDepth);
	}

	private double getPvalueValidSamples(int[] treatmentCytosines, int[] treatmentDepth,
										 int[] controlCytosines, int[] controlDepth) {


//		System.out.println(getMethylationValuesString(treatmentCytosines, treatmentDepth,
//				  controlCytosines, controlDepth));

		double theta = mleTheta(treatmentCytosines, treatmentDepth,
				controlCytosines, controlDepth);

//		System.err.println("MLE theta: "+theta);

		int[] cytosines = concat(treatmentCytosines, controlCytosines);
		int[] depth = concat(treatmentDepth, controlDepth);


		double mu = mleMu(cytosines, depth, theta);
//		System.err.println("MLE mu: "+mu);
		double muTreatment = mleMu(treatmentCytosines, treatmentDepth, theta);
//		System.err.println("MLE muTreatment: "+muTreatment);
		double muControl = mleMu(controlCytosines, controlDepth, theta);
//		System.err.println("MLE muControl: "+muControl);

		double likelihoodRatio = computeLikelihoodRatio(treatmentCytosines, treatmentDepth, controlCytosines,
				controlDepth, theta, mu, muTreatment, muControl);

//		System.err.println("Log likelihood ratio "+likelihoodRatio);
		return computePValue(likelihoodRatio, treatmentCytosines.length + controlCytosines.length);
	}

	private String getMethylationValuesString(int[] treatmentCytosines, int[] treatmentDepth, int[] controlCytosines,
											  int[] controlDepth) {
		StringBuilder methylationValuesSB = new StringBuilder();

		for (int i = 0; i < treatmentCytosines.length; i++) {
			methylationValuesSB.append(treatmentCytosines[i] + "/" + treatmentDepth[i] + " ");
		}

		methylationValuesSB.append(" vs. ");

		for (int i = 0; i < controlCytosines.length; i++) {
			methylationValuesSB.append(controlCytosines[i] + "/" + controlDepth[i] + " ");
		}

		return methylationValuesSB.toString().trim();
	}


	private double computePValue(double likelihoodRatio, int sampleSize) {

		return new TDistribution(sampleSize).cumulativeProbability(-1.0 * Math.sqrt(Math.max(0.0d, likelihoodRatio)))
				* 2;
	}

	private int[] concat(int[] head, int[] tail) {
		int[] toret = new int[head.length + tail.length];

		System.arraycopy(head, 0, toret, 0, head.length);
		System.arraycopy(tail, 0, toret, head.length, tail.length);

		return toret;
	}

	private double computeLikelihoodRatio(int[] treatmentCytosines, int[] treatmentDepth, int[] controlCytosines,
										  int[] controlDepth, double theta, double mu, double muTreatment, double
												  muControl) {

		int[] cytosines = concat(treatmentCytosines, controlCytosines);
		int[] depth = concat(treatmentDepth, controlDepth);


		return 2.0d * (
				logLikelihood(treatmentCytosines, treatmentDepth, muTreatment, theta) +
						logLikelihood(controlCytosines, controlDepth, muControl, theta)
						-
						logLikelihood(cytosines, depth, mu, theta));
	}

	private double logLikelihood(int[] cytosines, int[] depth, double mu, double theta) {

		double logLikelihood = 0.0d;

		for (int j = 0; j < cytosines.length; j++) {
			logLikelihood += Gamma.logGamma(mu * theta + cytosines[j]) + Gamma.logGamma((1 - mu) * theta + depth[j] -
					cytosines[j]) +
					Gamma.logGamma(theta) - Gamma.logGamma(mu * theta) - Gamma.logGamma((1 - mu) * theta) - Gamma
					.logGamma(depth[j] + theta);
		}
		return logLikelihood;
	}

	private double mleMu(int[] cytosines, int[] depth, double theta) {
		final UnivariateFunction function = new UnivariateFunction() {
			private double muEst;

			{
				muEst = (double) sumArray(cytosines) / (double) sumArray(depth);
//				System.err.println("muEst "+ muEst);
			}

			@Override
			public double value(double mu) {

				return muDerivative(mu, cytosines, depth, theta, muEst);
			}
		};

		final double relativeAccuracy = 1.0e-12;
		final double absoluteAccuracy = 1.0e-8;
		final int maxOrder = 5;
		UnivariateSolver solver = new BracketingNthOrderBrentSolver(relativeAccuracy, absoluteAccuracy, maxOrder);
		if (ArrayUtils.sumArray(cytosines) == ArrayUtils.sumArray(depth)) {
			return 1.0d - 1e-10;
		} else if (ArrayUtils.sumArray(cytosines) == 0) {
			return 0.0d + 1e-10;
		} else return solver.solve(1000, function, 0.0d + 1e-10, 1.0d - 1e-10);

	}

	private double mleTheta(int[] treatmentCytosines, int[] treatmentDepth, int[] controlCytosines, int[]
			controlDepth) {
		final UnivariateFunction function = new UnivariateFunction() {

			@Override
			public double value(double theta) {
				// TODO Auto-generated method stub
				return thetaDerivative(theta, treatmentCytosines, treatmentDepth, controlCytosines, controlDepth);
			}
		};

		final double relativeAccuracy = 1.0e-6;
		final double absoluteAccuracy = 1.0e-3;
		final int maxOrder = 5;
		UnivariateSolver solver = new BracketingNthOrderBrentSolver(relativeAccuracy, absoluteAccuracy, maxOrder);

		double[] interval = new double[]{0.001, 1e6};

		double estPhi = 0.0d;
//		System.err.println("thetaDerivative("+interval[0]+") = "+function.value(interval[0]));
//		System.err.println("thetaDerivative("+interval[1]+") = "+function.value(interval[1]));

		if (Double.isNaN(function.value(interval[1])) || function.value(interval[1]) >= 0) {
			estPhi = interval[1];

		} else if (function.value(interval[0]) <= 0) {
			estPhi = interval[0];
		} else {
			estPhi = solver.solve(1000, function, interval[0], interval[1]);
		}

		return estPhi;
	}

	private double thetaDerivative(double theta, int[] treatmentCytosines, int[] treatmentDepth, int[]
			controlCytosines, int[] controlDepth) {
//		System.err.println("computing derivative for theta "+theta);
		double muTreatmentEst = (double) sumArray(treatmentCytosines) / (double) sumArray(treatmentDepth);
		double muControlEst = (double) sumArray(controlCytosines) / (double) sumArray(controlDepth);

		double derivative = 0.0d;

		//treatment
		for (int j = 0; j < treatmentCytosines.length; j++) {
			double currentValue = muTreatmentEst * Gamma.digamma(muTreatmentEst * theta + treatmentCytosines[j]) +
					(1 - muTreatmentEst) * Gamma.digamma((1 - muTreatmentEst) * theta + treatmentDepth[j] -
							treatmentCytosines[j]) +
					Gamma.digamma(theta) -
					muTreatmentEst * Gamma.digamma(muTreatmentEst * theta)
					- (1 - muTreatmentEst) * Gamma.digamma((1 - muTreatmentEst) * theta)
					- Gamma.digamma(treatmentDepth[j] + theta);
//			System.err.println("current value "+j+" "+currentValue);
			derivative += currentValue;
		}

		//control
		for (int j = 0; j < controlCytosines.length; j++) {
			double currentValue = muControlEst * Gamma.digamma(muControlEst * theta + controlCytosines[j]) +
					(1 - muControlEst) * Gamma.digamma((1 - muControlEst) * theta + controlDepth[j] -
							controlCytosines[j]) +
					Gamma.digamma(theta) -
					muControlEst * Gamma.digamma(muControlEst * theta)
					- (1 - muControlEst) * Gamma.digamma((1 - muControlEst) * theta)
					- Gamma.digamma(controlDepth[j] + theta);

//			System.err.println("control current value "+j+" "+currentValue);
			derivative += currentValue;
		}

		return derivative;
	}

	private double muDerivative(double mu, int[] cytosines, int[] depth, double thetaEst, double muEst) {

		double derivative = 0.0d;
//		System.err.print("computing derivative for mu "+mu+", muest "+muEst+". ");
		for (int j = 0; j < cytosines.length; j++) {
			derivative += Gamma.digamma(muEst * thetaEst + cytosines[j]) - Gamma.digamma((1 - mu) * thetaEst +
					depth[j] - cytosines[j])
					- Gamma.digamma(mu * thetaEst) + Gamma.digamma((1 - mu) * thetaEst);
		}
//		System.err.println("computing mu derivative "+mu+ " "+derivative);
		return derivative;
	}

}
