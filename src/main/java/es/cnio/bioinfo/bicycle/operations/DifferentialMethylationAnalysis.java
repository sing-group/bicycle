package es.cnio.bioinfo.bicycle.operations;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.logging.Logger;

import org.apache.commons.lang.ArrayUtils;

import es.cnio.bioinfo.bicycle.MethylationCall;
import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Reference;
import es.cnio.bioinfo.bicycle.Sample;
import es.cnio.bioinfo.bicycle.gatk.GPFilesReader;
import es.cnio.bioinfo.pileline.core.Interval;
import es.cnio.bioinfo.pileline.core.IntervalsIndex;
import es.cnio.bioinfo.pileline.core.IntervalsIndexFactory;
import es.uvigo.ei.sing.math.statistical.corrections.FDRCorrection;

public class DifferentialMethylationAnalysis {
	private static final Logger logger = Logger.getLogger(DifferentialMethylationAnalysis.class.getName());
	
	private Project project;
	
	private MethylationAnalysis ma;
	
	private File outputFile;

	
	public DifferentialMethylationAnalysis(MethylationAnalysis ma) {
		System.out.println(ma);
		
		this.ma = ma;
		
	}
	
	public void analyzeDifferentialMethylationByRegions(
			Reference reference,
			List<Sample> treatmentSamples,
			List<Sample> controlSamples,
			File bedFile) throws NumberFormatException, IOException
	{
		logger.info("Performing differential methylation at region level (DMR)");
		logger.info("reference: "+reference.getReferenceFile());
		logger.info("treament samples: "+treatmentSamples);
		logger.info("control samples: "+controlSamples);
		logger.info("regions: "+bedFile);
		
		
		File outputFile = getDifferentiallyMethylatedRegionsFile(reference, treatmentSamples, controlSamples, bedFile);
		
		this.outputFile = outputFile;

		File tempFile = new File(outputFile.getAbsolutePath()+".temp");
		tempFile.deleteOnExit();
		
		PrintStream outTemp = new PrintStream(tempFile);
		writeOutputHeadersByRegion(reference, treatmentSamples, controlSamples, outTemp);
		
		Map<Interval, MethylationCounts> regionCounts = computeMethylationCountsByRegion(reference, treatmentSamples,
				controlSamples, bedFile);
		
		List<Double> pValues = new LinkedList<Double>();		
		

		//****************************************************//
		//modified by Osvaldo, 9Apr2017
		//adds methylation average values for treatment and control samples.
		//Additionally it computes a log2FC(treatment average/control average)
		
		
		//OLD way
		/*logger.info("Computing DMRs...");
		for (Interval interval : regionCounts.keySet()) {
			
			MethylationCounts regionCount = regionCounts.get(interval);
			
			outTemp.print(interval.getData()+"\t"+interval.getSequence()+"\t"+interval.getStart()+"\t"+interval.getStop()+"\t");
			outTemp.print(regionCount.toString());
			
			double pValue = this.computePValueForBase(regionCount);
			pValues.add(pValue);
			
			outTemp.println(pValue);			
		}*/
		
		
		//NEW way		
		int nTreatmentSamples=treatmentSamples.size();
		int nControlSamples=controlSamples.size();
		
		
		logger.info("Computing DMRs...");
		for (Interval interval : regionCounts.keySet()) {
			
			MethylationCounts regionCount = regionCounts.get(interval);
			
			outTemp.print(interval.getData()+"\t"+interval.getSequence()+"\t"+interval.getStart()+"\t"+interval.getStop()+"\t");
			outTemp.print(regionCount.toString());
			
			//added by Osvaldo
			double treatmentAverage=0.0;
			double controlAverage=0.0;
			int counter=0;
			//splits the methylation line by tab to get individual methylation ratios
			String[] tokens=regionCount.toString().split("\t");
			
			/*Two possible ways of calculating it:
			 * (1) calculating the ratio of methylation for each sample of each condition (numerator/denominator),
			 * and then calculating the average of methylation of each condition as the sum of its sample values
			 * divided by the number of samples.
			 * treatment average=((sample1 numerator/sample1 denominator) + (sample 2 numerator/sample 2 denominator)) /2
			 * control average=((sample1 numerator/sample1 denominator) + (sample 2 numerator/sample 2 denominator)) /2 
			 * Finally, it computes the log2FC (treatment average / control average)
			 * This way is implemented below
			 */
			
			// FIRST WAY DISCARDED !!!
			/*for(String token:tokens){// for each methylation ratio
				
				//treatment samples
				if(counter<nTreatmentSamples){
					//separates values that come represented as a methylation ratio
					String[] subToken=token.split("/");
					double numerator=new Double(subToken[0]).doubleValue();
					double denominator=new Double(subToken[1]).doubleValue();
					treatmentAverage+=new Double(numerator/denominator).doubleValue();
					counter++;
					
				}else{//control samples
					//separates values that come represented as a methylation ratio
					String[] subToken=token.split("/");
					double numerator=new Double(subToken[0]).doubleValue();
					double denominator=new Double(subToken[1]).doubleValue();
					controlAverage+=new Double(numerator/denominator).doubleValue();					
				}			
			}
			
			treatmentAverage=treatmentAverage/nTreatmentSamples;
			controlAverage=controlAverage/nControlSamples;
			
			double log2FC=Math.log(treatmentAverage/controlAverage)/Math.log(2);
			outTemp.print(treatmentAverage+"\t");
			outTemp.print(controlAverage+"\t");
			outTemp.print(log2FC+"\t");
			*/
			/*##################################
			 * END of type (1) calculation
			 */
			
			
			/*(2) another different way:
			 * for each condition, it sums the numerator values of the samples belonging to it and computes
			 * the average of the numerator values.
			 * It then sums the denominator values of the samples of the same condition and computes the average
			 * of the denominator values.
			 * treatment numerator=(sample1 numerator + sample 2 numerator) / 2
			 * treatment denominator=(sample1 denominator + sample 2 denominator) / 2
			 * treatment average=treatment numerator/treatment denominator
			 * 
			 * A similar calculaton is done for the control samples. Finally, it computes the
			 * log2FC (treatment average / control average)
			 * This second way is implemented below
			 */
			
			// SECOND WAY USED !!!
			//calculates average values of methylation for treatment and control samples
			double treatmentNumerator=0.0;
			double treatmentDenominator=0.0;
			double controlNumerator=0.0;
			double controlDenominator=0.0;
			treatmentAverage=0.0;
			controlAverage=0.0;
			counter=0;
			
			//calculates average values of methylation for treatment and control samples
			for(String token:tokens){// for each methylation ratio
				
				//treatment samples
				if(counter<nTreatmentSamples){
					//separates values that come represented as a methylation ratio
					String[] subToken=token.split("/");
					treatmentNumerator+=new Double(subToken[0]).doubleValue();
					treatmentDenominator+=new Double(subToken[1]).doubleValue();
					counter++;
					
				}else{//control samples
					//separates values that come represented as a methylation ratio
					String[] subToken=token.split("/");
					controlNumerator+=new Double(subToken[0]).doubleValue();
					controlDenominator+=new Double(subToken[1]).doubleValue();					
				}			
			}
			//computes averages
			treatmentNumerator=new Double(treatmentNumerator/nTreatmentSamples).doubleValue();
			treatmentDenominator=new Double(treatmentDenominator/nTreatmentSamples).doubleValue();
			controlNumerator=new Double(controlNumerator/nControlSamples).doubleValue();
			controlDenominator=new Double(controlDenominator/nControlSamples).doubleValue();
			
			treatmentAverage=treatmentNumerator/treatmentDenominator;
			controlAverage=controlNumerator/controlDenominator;
			
			double log2FC=Math.log(treatmentAverage/controlAverage)/Math.log(2);
			outTemp.print(treatmentAverage+"\t");
			outTemp.print(controlAverage+"\t");
			outTemp.print(log2FC+"\t");			
			/*##################################
			 * END of type (2) calculation
			 */
			
				
			double pValue = this.computePValueForBase(regionCount);
			pValues.add(pValue);
			
			outTemp.println(pValue);
		}
		// ENDS modification by Osvaldo
		//****************************************************//
		
		logger.info("[OK]");
		
		logger.info("Writing output file with adjusted p-values");
		createOutputFileWithAdjustedPValues(tempFile, pValues);
		logger.info("[OK]");
	}

	public File getDifferentiallyMethylatedRegionsFile(Reference reference, List<Sample> treatmentSamples,
			List<Sample> controlSamples, File bedFile) {
		String baseOutputFileName = computeOutputFileBaseName(reference, treatmentSamples, controlSamples);
		File outputFile = new File(baseOutputFileName+"."+bedFile.getName()+"-DMR.tsv");
		return outputFile;
	}

	public Map<Interval, MethylationCounts> computeMethylationCountsByRegion(Reference reference,
			List<Sample> treatmentSamples, List<Sample> controlSamples, File bedFile)
			throws IOException, FileNotFoundException {
		Map<Interval, MethylationCounts> regionCounts = new HashMap<>();
		IntervalsIndex bedIndex = IntervalsIndexFactory.createIntervalsIndex(bedFile.getName(), bedFile, 1, 2, 3, 1, false);
		
		int treatmentIndex = 0;
		for (Sample s: treatmentSamples) {
			Scanner scanner = new Scanner(ma.getMethylcytosinesFile(reference, s));
			while (scanner.hasNextLine()) {
				String line = scanner.nextLine();
				if (line.startsWith("#")) continue;
				MethylationCall call = MethylationCall.unmarshall(line);
				for (Interval interval : asIterable(bedIndex.getOverlappingIntervals(call.getContig(), (int) call.getPosition(), (int) call.getPosition()))){
					MethylationCounts counts = regionCounts.get(interval);
					if (counts == null) {
						counts = new MethylationCounts(treatmentSamples.size(), controlSamples.size());
						regionCounts.put(interval, counts);
					}
					counts.incrementTreatmentCytosines(treatmentIndex, call.getCytosines());
					counts.incrementTreatmentDepth(treatmentIndex, call.getCTdepth());
				}
			}
			scanner.close();
			treatmentIndex ++;
		}
		
		int controlIndex = 0;
		for (Sample s: controlSamples) {
			Scanner scanner = new Scanner(ma.getMethylcytosinesFile(reference, s));
			while (scanner.hasNextLine()) {
				String line = scanner.nextLine();
				if (line.startsWith("#")) continue;
				MethylationCall call = MethylationCall.unmarshall(line);
				for (Interval interval : asIterable(bedIndex.getOverlappingIntervals(call.getContig(), (int) call.getPosition(), (int) call.getPosition()))){
					MethylationCounts counts = regionCounts.get(interval);
					if (counts == null) {
						counts = new MethylationCounts(treatmentSamples.size(), controlSamples.size());
						regionCounts.put(interval, counts);
					}
					counts.incrementControlCytosines(controlIndex, call.getCytosines());
					counts.incrementControlDepth(controlIndex, call.getCTdepth());
				}
			}
			scanner.close();
			controlIndex ++;
		}
		return regionCounts;
	}
	
	
	private static <T> Iterable<T> asIterable(Iterator<T> t) {
		return new Iterable<T>(){
			@Override
			public Iterator<T> iterator() {
				return t;
			} 
		};
	}

	public void analyzeDifferentialMethylationByBase(
			Reference reference,
			List<Sample> treatmentSamples,
			List<Sample> controlSamples) throws IOException
	{
		
		logger.info("Performing differential methylation at base level (DMC)");
		logger.info("reference: "+reference.getReferenceFile());
		logger.info("treament samples: "+treatmentSamples);
		logger.info("control samples: "+controlSamples);
		
		File outputFile = getDifferentiallyMethylatedCytosinesFile(reference, treatmentSamples, controlSamples);
		this.outputFile = outputFile;
		
		File tempFile = new File(outputFile.getAbsolutePath()+".temp");
		tempFile.deleteOnExit();
		
		PrintStream outTemp = new PrintStream(tempFile);
		
		writeOutputHeadersByBase(reference, treatmentSamples, controlSamples, outTemp);

		BufferedReader[] sampleFiles = new BufferedReader[treatmentSamples.size() + controlSamples.size()];
		int i = 0;
		for (Sample s : treatmentSamples) {
			sampleFiles[i++] = new BufferedReader(new FileReader(ma.getMethylcytosinesFile(reference, s)));
		}
		for (Sample s : controlSamples) {
			sampleFiles[i++] = new BufferedReader(new FileReader(ma.getMethylcytosinesFile(reference, s)));
		}
		
		
		GPFilesReader reader = new GPFilesReader(reference.getSequenceNames(), sampleFiles);
		
		String currentSeq = null;
		long currentPos = -1;
		Map<Sample, MethylationCall> currentBaseCalls = new HashMap<>();
		
		String line = null;
		
		List<Double> pValues = new LinkedList<>();
		
		logger.info("Computing DMCs...");
		while ((line = reader.readLine()) != null) {
			if (line.startsWith("#"))
				continue;
			
			String[] tokens = line.split("\t");
			String lineSeq = tokens[0];
			long linePos = Long.parseLong(tokens[1]);

			if (currentSeq != null && (!lineSeq.equals(currentSeq) || linePos != currentPos)) {
				processBase(currentSeq, currentPos, treatmentSamples, controlSamples, currentBaseCalls, pValues, outTemp);
				
				currentBaseCalls.clear();
			}
			
			currentSeq = lineSeq;
			currentPos = linePos;
			Sample lineSample = 
					getSample(treatmentSamples, controlSamples, 
							reader.getLastLineReaderIndex());
			currentBaseCalls.put(lineSample, MethylationCall.unmarshall(line));
			
		}
		
		//process the last base
		if (currentSeq != null) {
			processBase(currentSeq, currentPos, treatmentSamples, controlSamples, currentBaseCalls, pValues, outTemp);
		}
		outTemp.close();
		logger.info("[OK]");
		
		logger.info("Writing output file with adjusted p-values");
		createOutputFileWithAdjustedPValues(tempFile, pValues);
		logger.info("[OK]");
	}

	public File getDifferentiallyMethylatedCytosinesFile(Reference reference, List<Sample> treatmentSamples,
			List<Sample> controlSamples) {
		String baseOutputFileName = computeOutputFileBaseName(reference, treatmentSamples, controlSamples);
		File outputFile = new File(baseOutputFileName+".DMC.tsv");
		return outputFile;
	}

	private void createOutputFileWithAdjustedPValues(File tempFile, List<Double> pValues) throws FileNotFoundException {
		String line;
		//adjust p-values
		double[] qValues = adjustPValues(pValues);
		
		//postprocess file adding q-values
		PrintStream out = new PrintStream(this.outputFile);
		
		Scanner scanner = new Scanner(tempFile);
		int lineNo = 0;
		while (scanner.hasNextLine()) {
			line = scanner.nextLine();
			if (!line.startsWith("#")) {
				out.println(line+"\t"+qValues[lineNo++]);
			} else {
				out.println(line);
			}
		}
		
		out.close();
	}

	private void writeOutputHeadersByBase(Reference reference, List<Sample> treatmentSamples,
			List<Sample> controlSamples, PrintStream outTemp)
			throws FileNotFoundException {
		outTemp.print("#SEQ\tPOS\tCONTEXT");
		
		writeSampleNames(treatmentSamples, controlSamples, outTemp);
		
		//modified by Osvaldo 11Apr2017
		//outTemp.println("\tp-value\tq-value");
		outTemp.println("\ttreatment average\tcontrol average\tlog2FC(treament/control)\tp-value\tq-value");
		
	}
	
	private void writeOutputHeadersByRegion(Reference reference, List<Sample> treatmentSamples,
			List<Sample> controlSamples, PrintStream outTemp)
			throws FileNotFoundException {
		outTemp.print("#region\tsequence\tstart\tstop");
		
		writeSampleNames(treatmentSamples, controlSamples, outTemp);
		
		//modified by Osvaldo 9Apr2017
		//outTemp.println("\tp-value\tq-value");
		outTemp.println("\ttreatment average\tcontrol average\tlog2FC(treament/control)\tp-value\tq-value");
	}

	public void writeSampleNames(List<Sample> treatmentSamples, List<Sample> controlSamples, PrintStream outTemp) {
		for (Sample s : treatmentSamples) {
			outTemp.print("\t"+s.getName()+" (treatment)");
		}
		for (Sample s : controlSamples) {
			outTemp.print("\t"+s.getName()+" (control)");
		}
	}

	private String computeOutputFileBaseName(Reference reference, List<Sample> treatmentSamples,
			List<Sample> controlSamples) {
		StringBuilder fileNameSB = new StringBuilder();
		fileNameSB.append(this.ma.getProject().getOutputDirectory()+File.separator+reference.getReferenceFile().getName());
		for (Sample sample : treatmentSamples) {
			fileNameSB.append("_"+sample.getName());
		}
		fileNameSB.append("__VS__");
		for (Sample sample : controlSamples) {
			fileNameSB.append("_"+sample.getName());
		}
		
		String baseOutputFileName = fileNameSB.toString();
		return baseOutputFileName;
	}

	private double[] adjustPValues(List<Double> pValues) {
		double[] pValuesArray = ArrayUtils.toPrimitive(pValues.toArray(new Double[0]));
		
		double[] qValues = null;
		FDRCorrection correction = new FDRCorrection();
		
		try {
			qValues = correction.correct(pValuesArray);
		} catch (InterruptedException e) {
			//should not see this
			throw new RuntimeException(e);
		}
		return qValues;
	}

	private void processBase(String currentSeq, long currentPos, List<Sample> treatmentSamples,
			List<Sample> controlSamples, Map<Sample, MethylationCall> currentBaseCalls, List<Double> pValues, PrintStream outTemp) {
		outTemp.print(currentSeq+"\t"+currentPos+"\t");
		outTemp.print(currentBaseCalls.get(controlSamples.get(0)).getContext().toString()+"\t");
		MethylationCounts counts = computeMethylationCounts(treatmentSamples, controlSamples, currentBaseCalls);
		outTemp.print(counts.toString());
	
		//****************************************************//
		//Added by Osvaldo, 11Apr2017
		//To include average methylation values for treatment and control, and log2FC(average treatment/average control)
		//SECOND WAY USED, as in analyzeDifferentialMethylationByRegions()
		//calculates average values of methylation for treatment and control samples
		int nTreatmentSamples=treatmentSamples.size();
		int nControlSamples=controlSamples.size();
		double treatmentNumerator=0.0;
		double treatmentDenominator=0.0;
		double controlNumerator=0.0;
		double controlDenominator=0.0;
		double treatmentAverage=0.0;
		double controlAverage=0.0;
		int counter=0;
		String[] tokens=counts.toString().split("\t");
		
		//calculates average values of methylation for treatment and control samples
		for(String token:tokens){// for each methylation ratio
			
			//treatment samples
			if(counter<nTreatmentSamples){
				//separates values that come represented as a methylation ratio
				String[] subToken=token.split("/");
				treatmentNumerator+=new Double(subToken[0]).doubleValue();
				treatmentDenominator+=new Double(subToken[1]).doubleValue();
				counter++;
				
			}else{//control samples
				//separates values that come represented as a methylation ratio
				String[] subToken=token.split("/");
				controlNumerator+=new Double(subToken[0]).doubleValue();
				controlDenominator+=new Double(subToken[1]).doubleValue();					
			}			
		}
		//computes averages
		treatmentNumerator=new Double(treatmentNumerator/nTreatmentSamples).doubleValue();
		treatmentDenominator=new Double(treatmentDenominator/nTreatmentSamples).doubleValue();
		controlNumerator=new Double(controlNumerator/nControlSamples).doubleValue();
		controlDenominator=new Double(controlDenominator/nControlSamples).doubleValue();
		
		treatmentAverage=treatmentNumerator/treatmentDenominator;
		controlAverage=controlNumerator/controlDenominator;
		
		double log2FC=Math.log(treatmentAverage/controlAverage)/Math.log(2);
		outTemp.print(treatmentAverage+"\t");
		outTemp.print(controlAverage+"\t");
		outTemp.print(log2FC+"\t");			
		/*##################################
		 * END of type (2) calculation
		 */
				
		double pValue = computePValueForBase(counts);
		outTemp.print(pValue);
		outTemp.println();
		
		pValues.add(pValue);
		
	}

	
	private class MethylationCounts {
		List<Integer> treatmentCytosines = new ArrayList<>();
		List<Integer> treatmentDepth = new ArrayList<>();
		List<Integer> controlCytosines = new ArrayList<>();
		List<Integer> controlDepth = new ArrayList<>();
		
		public MethylationCounts(int treatmentSize, int controlSize) {
			for (int i = 0; i < treatmentSize; i++) {
				treatmentCytosines.add(0);
				treatmentDepth.add(0);
			}
			for (int i = 0; i < controlSize; i++) {
				controlCytosines.add(0);
				controlDepth.add(0);
			}
		}
		
		public MethylationCounts(List<Integer> treatmentCytosines, List<Integer> treatmentDepth,
				List<Integer> controlCytosines, List<Integer> controlDepth) {
			super();
			this.treatmentCytosines = treatmentCytosines;
			this.treatmentDepth = treatmentDepth;
			this.controlCytosines = controlCytosines;
			this.controlDepth = controlDepth;
		}
		
		public void incrementTreatmentCytosines(int sampleIndex, int increment) {
			this.treatmentCytosines.set(sampleIndex, this.treatmentCytosines.get(sampleIndex) + increment);
		}
		
		public void incrementTreatmentDepth(int sampleIndex, int increment) {
			this.treatmentDepth.set(sampleIndex, this.treatmentDepth.get(sampleIndex) + increment);
		}

		public void incrementControlCytosines(int sampleIndex, int increment) {
			this.controlCytosines.set(sampleIndex, this.controlCytosines.get(sampleIndex) + increment);
		}
		
		public void incrementControlDepth(int sampleIndex, int increment) {
			this.controlDepth.set(sampleIndex, this.controlDepth.get(sampleIndex) + increment);
		}
		
		@Override
		public String toString() {
			StringBuilder methylationCountsSB = new StringBuilder();
			for (int i = 0; i < treatmentCytosines.size(); i++) {
				methylationCountsSB.append(treatmentCytosines.get(i)+"/"+treatmentDepth.get(i)+"\t");
			}
			
			for (int i = 0; i < controlCytosines.size(); i++) {
				methylationCountsSB.append(controlCytosines.get(i)+"/"+controlDepth.get(i)+"\t");
			}
			
			return methylationCountsSB.toString();
		}
		
		
	}
	
	private MethylationCounts computeMethylationCounts(List<Sample> treatmentSamples,
			List<Sample> controlSamples, Map<Sample, MethylationCall> baseCalls) {
		List<Integer> treatmentCytosines = new ArrayList<>();
		List<Integer> treatmentDepth = new ArrayList<>();
		List<Integer> controlCytosines = new ArrayList<>();
		List<Integer> controlDepth = new ArrayList<>();
		
		
		for (Sample sample: treatmentSamples) {
			if (baseCalls.containsKey(sample)) {
				treatmentCytosines.add(baseCalls.get(sample).getCytosines());
				treatmentDepth.add(baseCalls.get(sample).getCTdepth());
			}
		}
		
		for (Sample sample: controlSamples) {
			if (baseCalls.containsKey(sample)) {
				controlCytosines.add(baseCalls.get(sample).getCytosines());
				controlDepth.add(baseCalls.get(sample).getCTdepth());
			}
		}
		
		return new MethylationCounts(treatmentCytosines, treatmentDepth, controlCytosines, controlDepth);
	}
	private double computePValueForBase(MethylationCounts counts) {
		
		int[] treatmentCytosinesArray = ArrayUtils.toPrimitive(counts.treatmentCytosines.toArray(new Integer[]{}));
		int[] treatmentDepthArray = ArrayUtils.toPrimitive(counts.treatmentDepth.toArray(new Integer[]{}));
		int[] controlCytosinesArray = ArrayUtils.toPrimitive(counts.controlCytosines.toArray(new Integer[]{}));
		int[] controlDepthArray = ArrayUtils.toPrimitive(counts.controlDepth.toArray(new Integer[]{}));
		
		BetaBinomialDifferentialMethylationTest test = new BetaBinomialDifferentialMethylationTest();
		
		try {
			return test.getPvalue(treatmentCytosinesArray, treatmentDepthArray, controlCytosinesArray, controlDepthArray);
			
		} catch(Exception e) {
			return Double.NaN;
		}
	}


	private Sample getSample(List<Sample> treatmentSamples, List<Sample> controlSamples, int index) {
		if (index < treatmentSamples.size()) {
			return treatmentSamples.get(index);
		} else {
			return controlSamples.get(index - treatmentSamples.size());
		}
	}
	
	public Project getProject() {
		return project;
	}
	
	
	
}
