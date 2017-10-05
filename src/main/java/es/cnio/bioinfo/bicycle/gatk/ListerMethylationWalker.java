/*

Copyright 2012 Daniel Gonzalez Pe��a, Osvaldo Gra��a


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

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.RandomAccessFile;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.UUID;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;
import org.broad.tribble.Feature;
import org.broad.tribble.bed.BEDFeature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.IntervalBinding;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Downsample;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import es.cnio.bioinfo.bicycle.MethylationCall;
import net.sf.picard.filter.SamRecordFilter;
import net.sf.picard.util.SamLocusIterator.RecordAndOffset;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;


class MethylationFilePair {
	private static int instancecount = 0;
	private File watsonFile;
	private File crickFile;

	//private PrintStream watsonStream, crickStream;

	public File getWatsonFile() {
		return watsonFile;
	}

	public File getCrickFile() {
		return crickFile;
	}

	public MethylationFilePair(File watson, File crick) {
		instancecount++;
		this.watsonFile = watson;
		this.crickFile = crick;
		watsonPvals.put(Context.CG, new HashMap<Double, Integer>());
		watsonPvals.put(Context.CHG, new HashMap<Double, Integer>());
		watsonPvals.put(Context.CHH, new HashMap<Double, Integer>());
		crickPvals.put(Context.CG, new HashMap<Double, Integer>());
		crickPvals.put(Context.CHG, new HashMap<Double, Integer>());
		crickPvals.put(Context.CHH, new HashMap<Double, Integer>());


	}

	public MethylationFilePair(File watson, File crick, Map<Context, Map<Double, Integer>> watsonPvals, Map<Context,
			Map<Double, Integer>> crickPvals) {
		this(watson, crick);
		this.watsonPvals = watsonPvals;
		this.crickPvals = crickPvals;
	}

	private static int PRINTS_BEFORE_FLUSH = 50;
	private int printWatsonCount = 0;
	private int printCrickCount = 0;
	private StringBuffer watsonBuffer = new StringBuffer();
	private StringBuffer crickBuffer = new StringBuffer();

	private Map<Context, Map<Double, Integer>> watsonPvals = Collections.synchronizedMap(new HashMap<>());
	private Map<Context, Map<Double, Integer>> crickPvals = Collections.synchronizedMap(new HashMap<>());

	public void pushCall(MethylationCall call) {
		Map<Double, Integer> pvals = null;
		if (call.getStrand() == Strand.WATSON) {
			this.printWatson(call.marshall() + "\n");
			pvals = watsonPvals.get(call.getContext());
		} else {
			this.printCrick(call.marshall() + "\n");
			pvals = crickPvals.get(call.getContext());

		}

		Integer currentCount = pvals.get(call.getPval());
		if (currentCount == null) {
			currentCount = 0;
		}
		currentCount++;
		pvals.put(call.getPval(), currentCount);
	}

	public Map<Context, Map<Double, Integer>> getWatsonPvals() {
		return watsonPvals;
	}

	public Map<Context, Map<Double, Integer>> getCrickPvals() {
		return crickPvals;
	}

	private void printWatson(String str) {
		watsonBuffer.append(str);
		printWatsonCount++;
		if (printWatsonCount == PRINTS_BEFORE_FLUSH) {
			flushWatson();
		}
	}

	private void flushWatson() {
		PrintStream stream;
		try {
			stream = new PrintStream(new BufferedOutputStream(new FileOutputStream(watsonFile, true)));
			stream.print(watsonBuffer.toString());
			watsonBuffer.setLength(0);
			stream.flush();
			stream.close();
			printWatsonCount = 0;
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}
	}

	private void printCrick(String str) {
		crickBuffer.append(str);
		printCrickCount++;
		if (printCrickCount == PRINTS_BEFORE_FLUSH) {
			flushCrick();
		}
	}

	private void flushCrick() {
		PrintStream stream;
		try {
			stream = new PrintStream(new BufferedOutputStream(new FileOutputStream(crickFile, true)));
			stream.print(crickBuffer.toString());
			crickBuffer.setLength(0);
			stream.flush();
			stream.close();
			printCrickCount = 0;
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}
	}
	/*private PrintStream getWatsonStream() throws FileNotFoundException{
		if (watsonStream == null){
			watsonStream = new PrintStream(new BufferedOutputStream(new FileOutputStream(watsonFile)));
		}
		return this.watsonStream;
	}
	
	private PrintStream getCrickStream() throws FileNotFoundException{
		if (crickStream == null){
			crickStream = new PrintStream(new BufferedOutputStream(new FileOutputStream(crickFile)));
		}
		return this.crickStream;
	}*/

	public void flush() {
		flushWatson();
		flushCrick();
	}

	public void close() {
		instancecount--;
		flush();
		this.watsonBuffer = null;
		this.crickBuffer = null;
		//System.out.println("open reduces (MethylationFilePair): "+instancecount);
	}
}

@By(DataSource.READS)
@Reference(window = @Window(start = -2, stop = 2))
@Downsample(by = DownsampleType.NONE)
public class ListerMethylationWalker extends LocusWalker<List<MethylationCall>, MethylationFilePair> implements
		TreeReducible<MethylationFilePair> {


	@Argument(doc = "control genome for error computation. If parameter errorrate is also provided, this contig will " +
			"be skipped in methylcytosine call", required = false)
	public String controlGenome = "";

	@Argument(doc = "use a fixed error rate form <rate_watson>,<rate_crick>", required = false)
	public String errorRate = "";

	@Argument(doc = "FDR", required = false)
	public double FDR = 0.01;

	@Argument(doc = "correctNonCG", required = false)
	public boolean correctNonCG = false;

	@Argument(doc = "outdir", required = false)
	public File outdir = new File("./");

	@Argument(doc = "methylationwatsonfile", required = false)
	public File methylationwatsonfile = null;
	@Argument(doc = "methylationcrickfile", required = false)
	public File methylationcrickfile = null;

	@Argument(doc = "summaryfile", required = false)
	public File summaryfile = null;

	@Argument(doc = "methylcytosinesfile", required = false)
	public File methylcytosinesfile = null;

	@Argument(doc = "methylcytosinesvcffile", required = false)
	public File methylcytosinesvcffile = null;

	@Argument(doc = "remove clonal reads", required = false)
	public boolean removeClonal = false;

	@Argument(doc = "trim", required = false)
	public boolean trim = false;

	@Input(fullName = "annotation", shortName = "annotation", doc = "BED files to annotate methylcytosines", required
			= false)
	public List<RodBinding<BEDFeature>> beds = new ArrayList<RodBinding<BEDFeature>>();

	@Output
	public PrintStream out;

	private ContigBisulfiteError error;

	private HashMap<Strand, File> methylationFiles = new HashMap<Strand, File>();
	//private HashMap<Strand, PrintStream> methylationFilesOuts = new HashMap<Strand, PrintStream>();

	private HashMap<Strand, HashMap<Context, Double>> cutOffs = new HashMap<Strand, HashMap<Context, Double>>();

	private Tools tools = new Tools();
	private ListerFilter listerFilter;

	@Override
	public List<MethylationCall> map(RefMetaDataTracker metadata, ReferenceContext refContext, AlignmentContext
			alignmentContext) {
		if (refContext.getLocus().getContig().equals(controlGenome)) {
			return null;
		}

		if (refContext.getBase() == Strand.WATSON.getCytosineBase() || refContext.getBase() == Strand.CRICK
				.getCytosineBase()) {
			List<String> annotations = new LinkedList<String>();
			for (RodBinding<BEDFeature> binding : this.beds) {
				String annotation = "";
				List<BEDFeature> features = metadata.getValues(binding);

				if (features.size() == 0) {
					annotation = "N/A";
				} else {
					boolean first = true;
					for (BEDFeature feature : features) {

						if (!first) {
							annotation += "|";
						}
						annotation += feature.getName();
						first = false;
					}
				}
				annotations.add(annotation);

			}


			if (refContext.getBase() == Strand.WATSON.getCytosineBase()) { //WATSON
				List<MethylationCall> call = computeMethylationCall(refContext, alignmentContext, Strand.WATSON,
						Strand.CRICK, this.error, annotations);
				return call;
			} else { //CRICK
				List<MethylationCall> call = computeMethylationCall(refContext, alignmentContext, Strand.CRICK, Strand
						.WATSON, this.error, annotations);
				return call;

			}
		}
		return null;

	}

	@Override
	public void onTraversalDone(MethylationFilePair result) {
		super.onTraversalDone(result);
		result.flush();
		if (this.getToolkit().getArguments().numberOfThreads > 1) {
			//System.out.println("renaming "+result.CRICK+ " to "+this.methylationcrickfile);
			result.getCrickFile().renameTo(this.methylationcrickfile);
			//System.out.println("renaming "+result.WATSON+ " to "+this.methylationwatsonfile);
			result.getWatsonFile().renameTo(this.methylationwatsonfile);
		}

		String details = computePValCutOffs(result);

		try {
			GlobalMethylationStatistics stats = writeMethylCytosines();

			printSummary(stats, details);

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}

	}

	private String getOutputFilesPrefix() {
		StringBuilder toret = new StringBuilder();
		for (SAMReaderID readId : getToolkit().getReadsDataSource().getReaderIDs()) {
			toret.append(getToolkit().getReadsDataSource().getSAMFile(readId).getName());
		}
		return toret.toString();
	}

	private GlobalMethylationStatistics writeMethylCytosines() throws FileNotFoundException {

		PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(getMethylcytosinesfile())));
		PrintStream outvcf = new PrintStream(new BufferedOutputStream(new FileOutputStream(getMethylcytosinesVCFfile()
		)));

		try {
			BufferedReader wReader = new BufferedReader(new FileReader(this.methylationFiles.get(Strand.WATSON)));
			BufferedReader cReader = new BufferedReader(new FileReader(this.methylationFiles.get(Strand.CRICK)));

			List<String> sortedSequenceNames = toSequenceNames(super.getMasterSequenceDictionary());
			GPFilesReader reader = new GPFilesReader(sortedSequenceNames, wReader, cReader);

			String line = null;

			GlobalMethylationStatistics stats = new GlobalMethylationStatistics();

			writeMethylcytosinesHeader(out);
			writeVCFHeader(outvcf);
			while ((line = reader.readLine()) != null) {
				MethylationCall call = MethylationCall.unmarshall(line);

				writeMehylcytosinesRecord(out, call, stats);
				writeVCFRecord(outvcf, call);
			}
			out.close();
			outvcf.close();
			return stats;

		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

	}

	private List<String> toSequenceNames(SAMSequenceDictionary masterSequenceDictionary) {
		String[] sequenceNames = new String[masterSequenceDictionary.getSequences().size()];

		for (SAMSequenceRecord record : masterSequenceDictionary.getSequences()) {
			sequenceNames[record.getSequenceIndex()] = record.getSequenceName();
		}
		return Arrays.asList(sequenceNames);
	}

	private void writeMethylcytosinesHeader(PrintStream out) {
		out.print(MethylationCall.getMarshallHeader());
		for (RodBinding<BEDFeature> binding : this.beds) {
			out.print("\t" + binding.getName());
		}
		out.println("\tSTATUS");

	}

	private void writeMehylcytosinesRecord(PrintStream out, MethylationCall call,
										   GlobalMethylationStatistics stats) {
		double cutOff = this.cutOffs.get(call.getStrand()).get(call.getContext());
		call.setCutOff(cutOff);
		stats.add(call);

		call.marshall(out);
		if (call.getPval() < cutOff) {
			out.println("\tMETHYLATED");
		} else {
			out.println("\tUNMETHYLATED");
		}
	}


	private void writeVCFHeader(PrintStream out) {
		//header
		out.println("#fileformat=VCFv4.1");
		out.println("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
		out.println("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
		out.println("##INFO=<ID=CTDP,Number=1,Type=Integer,Description=\"CorT Depth\">");
		out.println("##INFO=<ID=CD,Number=1,Type=Integer,Description=\"Cytosine Depth\">");

		//modified (osvaldo, 3jan2016)
		//out.println("##INFO=<ID=PER,Number=1,Type=Float,Description=\"Methylation percentage\">");		
		out.println("##INFO=<ID=BS,Number=1,Type=Float,Description=\"Beta Score\">");

		out.println("##INFO=<ID=PU,Number=1,Type=Float,Description=\"Readed bases at this position\">");
		out.println("##INFO=<ID=CO,Number=1,Type=Flag,Description=\"Corrected, i.e., this CG is derived from a non-GC " +
				"to GC correction\">");
		out.println("##INFO=<ID=AC,Number=1,Type=Flag,Description=\"Added by correction, i.e., this CG is added due to" +
				" a correction from non-CG to CG in the opposite strand\">");
		out.println("##INFO=<ID=STR,Number=1,Type=String,Description=\"Strand Aligment\">");

		for (RodBinding<BEDFeature> binding : this.beds) {
			out.println("##INFO=<ID=" + binding.getName() + ",Number=1,Type=String,Description=\"" + binding.getName()
					+ " annotation\">");
		}

		out.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
	}

	private void writeVCFRecord(PrintStream out, MethylationCall call) {
		double cutOff = this.cutOffs.get(call.getStrand()).get(call.getContext());

		out.print(call.getContig() + "\t");
		out.print(call.getPosition() + "\t");
		out.print(call.getContext() + "\t");
		out.print("C\t");
		if (call.getPval() < cutOff) {
			out.print("C\t");
		} else {
			out.print(".\t");
		}
		out.print(".\t.\t");

		//info
		out.print("NS=1;");
		out.print("DP=" + call.getDepth() + ";");
		out.print("CTDP=" + call.getCTdepth() + ";");
		out.print("CD=" + call.getCytosines() + ";");

		//modified (osvaldo, 3jan2016)
		//out.print("PER="+new DecimalFormat("###.##").format(100*(double)call.getCytosines()/(double)call.getDepth())
		// +";");
		out.print("BS=" + new DecimalFormat("#.#######").format(call.getBetaScore()) + ";");

		out.print("PU=" + call.getPileup() + ";");
		if (call.isCorrectedFromNonCG()) {
			out.print("CO;");
		}
		if (call.isAddedByCorrection()) {
			out.print("AC;");
		}

		out.print("STR=" + (call.getStrand() == Strand.WATSON ? "+" : "-") + ";");

		for (int i = 0; i < this.beds.size(); i++) {
			out.print(this.beds.get(i).getName() + "=" + call.getAnnotations().get(i) + ";");
		}
		out.println();
	}

	private File getMethylcytosinesfile() {
		if (this.methylcytosinesfile == null) {
			return new File(this.outdir + "/" + getOutputFilesPrefix() + ".methylcytosines");
		} else {
			return this.methylcytosinesfile;
		}
	}

	private File getMethylcytosinesVCFfile() {
		if (this.methylcytosinesvcffile == null) {
			return new File(this.outdir + "/" + getMethylcytosinesfile().getName() + ".vcf");
		} else {
			return this.methylcytosinesvcffile;
		}
	}

	private void printSummary(GlobalMethylationStatistics stats, String cutoffDetails)
			throws FileNotFoundException {
		PrintStream summary = new PrintStream(new FileOutputStream(getSummaryFile()));
		summary.println("====METHYLATION RESULTS=======================================================");
		summary.println("File: " + getSummaryFile().getName());
		summary.println("Date: " + new Date());
		summary.println();

		summary.println("====ANALYSIS PARAMETERS=======================================================");
		summary.println(" Correct non-CG: " + this.correctNonCG);


		summary.println(" Filters:" + (listerFilter == null ? "\n" : "\n  " + listerFilter.toString().replace(",", "\n" +
				" ")));
		summary.println("  remove clonal reads: " + this.removeClonal);
		summary.println(" FDR threshold: " + this.FDR);
		summary.println();
		summary.println("====ERROR ESTIMATION AND SIGNIFICANCE ADJUSTMENTS=============================");
		//error rates
		summary.print(" Error rates (");
		if (!this.errorRate.equals("")) {
			summary.println("fixed):");
		} else if (!this.controlGenome.equals("")) {
			summary.println("from control genome: " + this.controlGenome + "):");
		}
		summary.println("  " + this.error.toString().replaceAll("\n", "\n  "));
		summary.println("  p-value cutoffs: " + this.cutOffs);
		summary.println();
		summary.println("====METHYLATION ANALYSIS RESULTS==============================================");
		//statistics		
		summary.println(stats);

		//cut-off details
		summary.println("Cut-off computation details:\n" + cutoffDetails);

	}

	private ListerFilter getListerFilter() {
		//has lister filters?
		ListerFilter listerFilter = null;
		for (SamRecordFilter filter : this.getToolkit().getFilters()) {
			if (filter instanceof ListerFilter) {
				listerFilter = (ListerFilter) filter;
				break;
			}
		}
		return listerFilter;
	}

	private File getSummaryFile() {
		if (this.summaryfile == null) {
			return new File(this.outdir + "/" + getOutputFilesPrefix() + ".summary");
		} else {
			return this.summaryfile;
		}
	}


	@Override
	public MethylationFilePair reduce(List<MethylationCall> arg0, MethylationFilePair arg1) {
		if (arg0 != null) {
			if (arg1 == null) {
				File outwatson = null, outcrick = null;
				if (this.getToolkit().getArguments().numberOfThreads > 1) {
					try {
						outwatson = new File(this.outdir + File.separator + this.methylationwatsonfile.getName() +
								UUID.randomUUID() + ".methylation");
						outwatson.deleteOnExit();
						new FileOutputStream(outwatson).close(); //touch

						outcrick = new File(this.outdir + File.separator + this.methylationcrickfile.getName() + UUID
								.randomUUID() + ".methylation");
						outcrick.deleteOnExit();
						new FileOutputStream(outcrick).close(); //touch
					} catch (FileNotFoundException e) {
						e.printStackTrace();
					} catch (IOException e) {
						throw new RuntimeException(e);
					}
				} else {
					outwatson = this.methylationwatsonfile;
					outcrick = this.methylationcrickfile;
				}
				arg1 = new MethylationFilePair(outwatson, outcrick);
			}
			for (MethylationCall call : arg0) {
				arg1.pushCall(call);

			}


		}

		return arg1;
	}

	@Override
	public MethylationFilePair reduceInit() {
		return null;
	}

	@Override
	public MethylationFilePair treeReduce(MethylationFilePair arg0, MethylationFilePair arg1) {
		//System.err.println("treereduce: "+arg0+" "+arg1);
		if (arg0 == null && arg1 != null) {
			arg1.flush();
			return arg1;
		}
		if (arg1 == null && arg0 != null) {
			arg0.flush();
			return arg0;
		}
		if (arg1 == null && arg0 == null) {
			return null;
		}


		try {

			File outwatson = new File(this.outdir + File.separator + this.methylationwatsonfile.getName() + UUID
					.randomUUID() + ".methylation");
			outwatson.deleteOnExit();
			File outcrick = new File(this.outdir + File.separator + this.methylationcrickfile.getName() + UUID
					.randomUUID() + ".methylation");
			outcrick.deleteOnExit();

			arg0.flush();
			arg1.flush();
			arg0.close();
			arg1.close();
			appendFiles(arg0.getWatsonFile(), arg1.getWatsonFile(), outwatson);
			appendFiles(arg0.getCrickFile(), arg1.getCrickFile(), outcrick);

			Map<Context, Map<Double, Integer>> watsonPvals = mergePvals(arg0.getWatsonPvals(), arg1.getWatsonPvals());
			Map<Context, Map<Double, Integer>> crickPvals = mergePvals(arg0.getCrickPvals(), arg1.getCrickPvals());
			return new MethylationFilePair(outwatson, outcrick, watsonPvals, crickPvals);

		} catch (IOException e) {
			throw new RuntimeException(e);
		}

	}

	private Map<Context, Map<Double, Integer>> mergePvals(
			Map<Context, Map<Double, Integer>> pvals1,
			Map<Context, Map<Double, Integer>> pvals2) {

		Map<Context, Map<Double, Integer>> toret = new HashMap<Context, Map<Double, Integer>>();
		toret.put(Context.CG, new HashMap<Double, Integer>());
		toret.put(Context.CHG, new HashMap<Double, Integer>());
		toret.put(Context.CHH, new HashMap<Double, Integer>());

		for (Context context : Context.values()) {
			HashMap<Double, Integer> toretPvals = new HashMap<Double, Integer>();
			toret.put(context, toretPvals);
			Set<Double> allPvals = new HashSet<Double>();
			allPvals.addAll(pvals1.get(context).keySet());
			allPvals.addAll(pvals2.get(context).keySet());
			for (Double pval : allPvals) {
				int count = 0;

				Integer count1 = pvals1.get(context).get(pval);
				if (count1 == null) {
					count1 = 0;
				}
				Integer count2 = pvals2.get(context).get(pval);
				if (count2 == null) {
					count2 = 0;
				}
				count = count1 + count2;
				toretPvals.put(pval, count);
			}
		}

		return toret;
	}

	@Override
	public void initialize() {
		super.initialize();
		this.listerFilter = this.getListerFilter();
		ListerFilter.trim = this.trim;

		if (this.controlGenome.equals("") && this.errorRate.equals("")) {
			throw new RuntimeException("Please provide at least --controlgenome or --erorrate");
		}
		if (!this.controlGenome.equals("") && this.errorRate.equals("")) {


			GATKArgumentCollection arguments = this.getToolkit().getArguments();
			List<IntervalBinding<Feature>> oldIntervals = arguments.intervals;


			arguments.intervals = new LinkedList<IntervalBinding<Feature>>();

			arguments.intervals.add(new IntervalBinding(controlGenome));


			out.println("computing error");

			ComputeErrorFromContig walker = (ComputeErrorFromContig) this.getToolkit().getWalkerByName
					("ComputeErrorFromContig");
			walker.tools = this.tools;
			walker.removeClones = this.removeClonal;


			this.getToolkit().setWalker(walker);
			this.getToolkit().execute();
			//this.getToolkit().execute(arguments, walker, this.getToolkit().getFilters());

			this.error = walker.getLastTraversalResult();

			arguments.intervals = oldIntervals;

			if (this.listerFilter != null) {
				this.listerFilter.resetCounters();
			}


		} else {
			final double WATSON_ERROR = Double.parseDouble(this.errorRate.split(",")[0]);
			final double CRICK_ERROR = Double.parseDouble(this.errorRate.split(",")[1]);

			this.error = new ContigBisulfiteError() {
				@Override
				public BisulfiteError getError(final Strand strand, Context context) {
					return new BisulfiteError() {
						@Override
						public double getError() {
							if (strand == Strand.WATSON) {
								return WATSON_ERROR;
							} else if (strand == Strand.CRICK) {
								return CRICK_ERROR;
							}
							throw new RuntimeException("Incompatible strand " + strand);
						}

					};

				}

				@Override
				public String toString() {
					StringBuilder toret = new StringBuilder();
					for (Strand s : Strand.values()) {
						toret.append(s + " = {");
						boolean first = true;
						for (Context c : Context.values()) {
							if (!first) toret.append(", ");
							else first = false;
							toret.append(c + " = " + getError(s, c).getError());
						}
						toret.append("} ");
					}
					return toret.toString();
				}

			};
		}

		out.println("Error computed " + this.error);
		for (Strand strand : Strand.values()) {
			File file = getMethylationfile(strand);
			methylationFiles.put(strand, file);
		}
	}

	private File getMethylationfile(Strand strand) {
		if (strand == Strand.WATSON && this.methylationwatsonfile != null) {
			return this.methylationwatsonfile;
		}
		if (strand == Strand.CRICK && this.methylationcrickfile != null) {
			return this.methylationcrickfile;
		}
		return new File(this.outdir.getAbsolutePath() + "/" + getOutputFilesPrefix() + "_" + strand + ".methylation");
	}

	private List<MethylationCall> computeMethylationCall(
			ReferenceContext refContext,
			AlignmentContext alignmentContext,
			Strand strand,
			Strand oppositeStrand,
			ContigBisulfiteError error, List<String> annotations) {
		List<MethylationCall> toret = new LinkedList<MethylationCall>();

		Context context = strand.getContext(refContext, alignmentContext.getPosition());
		if (context == null) {
			return null;
		}
		ReadBackedPileup reads = ListerFilter.applyFilters(tools.getReadsForStrand(strand, alignmentContext,
				removeClonal));


		if (reads != null) {
			
			/* DEBUG */
			//System.out.println(strand+" "+alignmentContext.getContig()+":"+alignmentContext.getPosition()+":"+new
			// String(reads.getBases()));

			int mCCount = 0;

			int depth = reads.depthOfCoverage();
			int CTdepth = reads.depthOfCoverage();
			StringBuilder pileupb = new StringBuilder();

			for (byte base : reads.getBases()) {
				pileupb.append((char) base);
				if (base == strand.getCytosineBase()) {
					mCCount++;
				}
				if (base != strand.getCytosineBase() && base != strand.getThymineBase()) {
					CTdepth--;
				}
			}
			double pval = computePval(strand, context, error, mCCount, depth);

			//added (osvaldo, 3jan2016)
			double CRatio = (double) mCCount / (double) depth;

			//modified (osvaldo, 31dec2015): adds CRatio argument
			MethylationCall call = new MethylationCall(alignmentContext.getContig(), alignmentContext.getPosition(),
					strand, context, pval, depth, CTdepth, mCCount, pileupb.toString(), false, false, annotations, new
					Double(CRatio));
			toret.add(call);

			//perform nonCG to CG correction
			if (toret != null && context != Context.CG && this.correctNonCG) {

				int downstreamPosition = (int) strand.downstream(alignmentContext.getPosition());
				// use raw Picard api to see the downstream bases

				if (this.listerFilter != null) {
					this.listerFilter.freezeCountersInThread();
				}
				Collection<RecordAndOffset> strandReads = tools.getReadsForLocus(this.getToolkit(), strand,
						alignmentContext.getContig(), downstreamPosition, removeClonal);

				if (this.listerFilter != null){
					this.listerFilter.unfreezeCountersInThread();
				}
				int strandGCount = 0;
				int strandDepth = strandReads.size();
				double strandGRatio = 0;
				char guanineBase = strand.getGuanineBase();
				for (RecordAndOffset ro : strandReads) {
					if (ro.getReadBase() == guanineBase) {
						strandGCount++;
					}
				}
				strandGRatio = (double) strandGCount / (double) strandDepth;

				if (strandGRatio >= 0.2d) {

					if (this.listerFilter != null) {
						this.listerFilter.freezeCountersInThread();
					}
					Collection<RecordAndOffset> oppositeReads = tools.getReadsForLocus(this.getToolkit(),
							oppositeStrand, alignmentContext.getContig(), downstreamPosition, removeClonal);
					if (this.listerFilter != null){
						this.listerFilter.unfreezeCountersInThread();
					}
					/* DEBUG */
					/*StringBuilder s = new StringBuilder();
					for (RecordAndOffset ro : oppositeReads){
						s.append((char)ro.getReadBase());
					}
					System.out.println(oppositeStrand+" "+alignmentContext.getContig()+":"+downstreamPosition+":"+s
					.toString());
					*/
					if (!oppositeReads.isEmpty()) {
						//both G ratios above 20%?


						int oppositeCCount = 0;
						int oppositeDepth = oppositeReads.size();
						int oppositeCTdepth = oppositeDepth;
						StringBuilder oppositePileup = new StringBuilder();
						double oppositeCRatio = 0;
						for (RecordAndOffset ro : oppositeReads) {
							oppositePileup.append((char) ro.getReadBase());
							if (ro.getReadBase() == oppositeStrand.getCytosineBase()) {
								oppositeCCount++;
							}
							if (ro.getReadBase() != oppositeStrand.getCytosineBase() && ro.getReadBase() !=
									oppositeStrand.getThymineBase()) {
								oppositeCTdepth--;
							}
						}
						oppositeCRatio = (double) oppositeCCount / (double) oppositeDepth;

						if (oppositeCRatio >= 0.2) {
							call.correctToCG();

							//modified (osvaldo, 31dec2015): adds oppositeCRatio argument
							//add another one to the opposite strand
							MethylationCall oppositeCall = new MethylationCall(alignmentContext.getContig(),
									downstreamPosition, oppositeStrand, Context.CG, computePval(oppositeStrand,
									Context.CG, this.error, oppositeCCount, oppositeDepth), oppositeDepth,
									oppositeCTdepth, oppositeCCount, oppositePileup.toString(), false, true,
									annotations, oppositeCRatio);
							if (downstreamPosition < alignmentContext.getPosition()) {
								toret.add(0, oppositeCall); //prepend
							} else {
								toret.add(oppositeCall); //append
							}
						}

					} else if (Math.abs(strandGCount - call.getCytosines()) <= 2) {
						call.correctToCG();

					}
				}

			}
		}
		return toret;

	}

	private double computePval(Strand strand, Context context,
							   ContigBisulfiteError error, int mCCount, int depth) {
		BinomialDistribution binomial = new BinomialDistributionImpl(depth, error.getError(strand, context).getError
				());

		double pval = 1.0;

		try {
			pval = (mCCount == 0) ? 1.0d : (1.0d - binomial.cumulativeProbability(mCCount - 1));
		} catch (MathException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
		return pval;
	}


	private String computePValCutOffs(final MethylationFilePair results) {


		class CutOffThread extends Thread {
			double[] cutoffs;
			Strand strand;
			Throwable error;
			StringBuffer details = new StringBuffer();

			public CutOffThread(Strand strand) {
				this.strand = strand;
			}

			public void run() {
				try {
					cutoffs = computePValCutoffsFile(strand, results, details);
					//System.out.println("cutoffs: "+cutoffs);
				} catch (Exception e) {
					e.printStackTrace();
					this.error = e;
				}
			}
		}
		;

		CutOffThread watsonT = new CutOffThread(Strand.WATSON);
		CutOffThread crickT = new CutOffThread(Strand.CRICK);

		watsonT.start();
		crickT.start();

		try {

			for (CutOffThread thread : new CutOffThread[]{watsonT, crickT}) {
				thread.join();
				if (thread.error != null) {
					throw new RuntimeException(thread.error);
				}

				out.println("p-val cutoffs computed for strand " + thread.strand + " : " + Arrays.toString(thread
						.cutoffs));

				HashMap<Context, Double> cutOffs = new HashMap<Context, Double>();

				cutOffs.put(Context.CG, thread.cutoffs[0]);
				cutOffs.put(Context.CHG, thread.cutoffs[1]);
				cutOffs.put(Context.CHH, thread.cutoffs[2]);

				this.cutOffs.put(thread.strand, cutOffs);


			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		return "WATSON\n" + watsonT.details + "\nCRICK\n" + crickT.details;
	}

	private double[] computePValCutoffsFile(Strand strand, MethylationFilePair results, StringBuffer computingDetails)
			throws FileNotFoundException, IOException {
		File outFile = null;
		if (strand == Strand.WATSON) {
			outFile = results.getWatsonFile();
		} else {
			outFile = results.getCrickFile();
		}

		double[] positiveRate = {0d, 0d, 0d};
		double[] cutoffs = {1, 1, 1};
		boolean[] needAdjust = {true, true, true};

		int[] cCount = countCs(strand, results);


//			int[] mcCount = {0,0,0};

		int iteration = 0;
		while (needAdjust[0] || needAdjust[1] || needAdjust[2]) {

			iteration++;
			int[] mcCount = countMCs(strand, results, cutoffs);
				
				/*
				BufferedReader outIn = new BufferedReader(new FileReader(outFile));
				
				String line = null;
				while ((line = outIn.readLine())!=null){
					MethylationCall call = MethylationCall.unmarshall(line); 
					//String[] tokens = line.split("\t");
					if (call.getContext()==Context.CG){
						cCount[0]++;
						if (call.getPval()<cutoffs[0]){
							mcCount[0]++;
						}
					}else if (call.getContext()==Context.CHG){
						cCount[1]++;
						if (call.getPval()<cutoffs[1]){
							mcCount[1]++;
						}
					}else if (call.getContext()==Context.CHH){
						cCount[2]++;
						if (call.getPval()<cutoffs[2]){
							mcCount[2]++;
						}
					}
				}
				outIn.close();
				*/

			double[] _positiveRate = {((double) mcCount[0] / (double) cCount[0]) * 100, ((double) mcCount[1] /
					(double) cCount[1]) * 100, ((double) mcCount[2] / (double) cCount[2]) * 100};
			computingDetails.append("Iteration " + iteration + ", M: " + Arrays.toString(cutoffs) + " %mC: " + Arrays
					.toString(_positiveRate) + "\n");

			if (Double.isNaN(_positiveRate[0])) {
				cutoffs[0] = 0; // ??
				needAdjust[0] = false;
			} else if (_positiveRate[0] != positiveRate[0]) {
				cutoffs[0] = FDR * _positiveRate[0] / (100d - _positiveRate[0]);
				positiveRate[0] = _positiveRate[0];
			} else {
				needAdjust[0] = false;
			}


			if (Double.isNaN(_positiveRate[1])) {
				cutoffs[1] = 0; // ??
				needAdjust[1] = false;
			} else if (_positiveRate[1] != positiveRate[1]) {
				cutoffs[1] = FDR * _positiveRate[1] / (100d - _positiveRate[1]);
				positiveRate[1] = _positiveRate[1];
			} else {
				needAdjust[1] = false;
			}

			if (Double.isNaN(_positiveRate[2])) {
				cutoffs[2] = 0; // ??
				needAdjust[2] = false;
			} else if (_positiveRate[2] != positiveRate[2]) {
				cutoffs[2] = FDR * _positiveRate[2] / (100d - _positiveRate[2]);
				positiveRate[2] = _positiveRate[2];
			} else {
				needAdjust[2] = false;
			}


			computingDetails.append("\tneed Adjust: " + Arrays.toString(needAdjust) + "\n");
		}

		computingDetails.append("Finished p-value adjust. Result " + Arrays.toString(cutoffs) + "\n");
		return cutoffs;
	}

	private int[] countMCs(Strand strand, MethylationFilePair results, double[] cutoffs) {
		int[] toret = {0, 0, 0};
		Map<Context, Map<Double, Integer>> pvals = null;
		if (strand == Strand.WATSON) {
			pvals = results.getWatsonPvals();
		} else {
			pvals = results.getCrickPvals();
		}

		int i = 0;
		for (Context context : new Context[]{Context.CG, Context.CHG, Context.CHH}) {
			int total = 0;
			for (Double pval : pvals.get(context).keySet()) {
				if (pval < cutoffs[i]) {

					total += pvals.get(context).get(pval);
				}
			}
			toret[i] = total;
			i++;
		}
		return toret;
	}

	private int[] countCs(Strand strand, MethylationFilePair results) {
		int[] toret = {0, 0, 0};
		Map<Context, Map<Double, Integer>> pvals = null;
		if (strand == Strand.WATSON) {
			pvals = results.getWatsonPvals();
		} else {
			pvals = results.getCrickPvals();
		}

		int i = 0;
		for (Context context : new Context[]{Context.CG, Context.CHG, Context.CHH}) {
			int total = 0;
			for (Integer count : pvals.get(context).values()) {
				total += count;
			}
			toret[i] = total;
			i++;
		}
		return toret;
	}

	private void appendFiles(File f1, File f2, File outfile) throws IOException {
		//System.err.println("appending: "+f1+", "+f2+": "+outfile);

		byte[] buffer = new byte[2048];

		BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outfile));

		//file 1
		FileInputStream fis = new FileInputStream(f1);

		BufferedInputStream in = new BufferedInputStream(fis);
		int readed = -1;

		while ((readed = in.read(buffer)) != -1) {
			out.write(buffer, 0, readed);
		}
		in.close();
		fis.close();


		//file 2
		fis = new FileInputStream(f2);
		in = new BufferedInputStream(fis);
		readed = -1;

		while ((readed = in.read(buffer)) != -1) {
			out.write(buffer, 0, readed);
		}
		in.close();
		fis.close();

		if (!f1.delete()) {
			throw new RuntimeException("Could no delete tempary file: " + f1);
		}

		if (!f2.delete()) {
			throw new RuntimeException("Could no delete tempary file: " + f2);
		}

		out.flush();
		out.close();

	}

	public static void main(String[] args) throws IOException {
		File f = new File("/mnt/lacie15t/BACKUP/lipido/NGS/LISTER/REFERENCES/phageLambda_plus_hg18.fa.fai");
		RandomAccessFile ra = new RandomAccessFile(f, "rw");
		ra.getChannel().lock(0, Long.MAX_VALUE, false);

		//a
	/*	ra.seek(0);
		ra.write("hello2\n".getBytes());	
		ra.close();
		*/
		//b
		FileOutputStream fos = new FileOutputStream(f);
		fos.write("hello\n".getBytes());
		fos.close();


	}
}
