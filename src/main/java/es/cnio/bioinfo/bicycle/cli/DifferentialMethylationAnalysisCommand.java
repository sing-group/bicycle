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

package es.cnio.bioinfo.bicycle.cli;

import java.io.File;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Reference;
import es.cnio.bioinfo.bicycle.Sample;
import es.cnio.bioinfo.bicycle.gatk.Context;
import es.cnio.bioinfo.bicycle.operations.DifferentialMethylationAnalysis;
import es.cnio.bioinfo.bicycle.operations.MethylationAnalysis;

public class DifferentialMethylationAnalysisCommand extends ProjectCommand {

	@Override
	public String getName() {
		return "analyze-differential-methylation";
	}

	@Override
	public String getDescription() {
		return "Analyzes differential methylation for treatment-control samples at base or region level";
	}

	@Override
	public void executeImpl(CLIApplication app, Project project, Map<Option, String> parameters) throws Exception {
		
		List<Sample> treatmentSamples = new LinkedList<>();
		if (parameters.containsKey(this.findOption("t"))){
			final String treatmentSamplesString = parameters.get(this.findOption("t"));
			treatmentSamples= parseSamples(project, treatmentSamplesString);
		} else {
			String validSampleNames = getValidSampleNames(project);
			throw new IllegalArgumentException("treatment samples are mandatory. Valid sample names are: "+validSampleNames);
		}
		
		List<Sample> controlSamples = new LinkedList<>();
		if (parameters.containsKey(this.findOption("c"))){
			final String controlSamplesString = parameters.get(this.findOption("c"));
			controlSamples= parseSamples(project, controlSamplesString);
		} else {
			String validSampleNames = getValidSampleNames(project);
			throw new IllegalArgumentException("control samples are mandatory. Valid sample names are: "+validSampleNames);
		}
		
		
		List<File> bedFiles = new LinkedList<File>();
		if (parameters.containsKey(this.findOption("b"))){
			String bedFilesString = parameters.get(this.findOption("b"));
			if (bedFilesString!=null){
				String[] tokens = bedFilesString.split(",");
				
				for (String token: tokens){
					File bedFile = new File(token);
					if (!bedFile.exists()){
						throw new IllegalArgumentException("BED file not found: "+bedFile);
					}
					bedFiles.add(bedFile);
				}
			}
		} 

		Set<Context> contexts = parseContextParameter(parameters.get(this.findOption("x")));

		DifferentialMethylationAnalysis dma = 
				new DifferentialMethylationAnalysis(new MethylationAnalysis(project), contexts);
		for (Reference reference: project.getReferences()) {
			dma.analyzeDifferentialMethylationByBase(reference, treatmentSamples, controlSamples);
			//by region
			for (File bedFile: bedFiles) {
				dma.analyzeDifferentialMethylationByRegions(reference, treatmentSamples, controlSamples, bedFile);
			}
		}
	}

	private Set<Context> parseContextParameter(String contextsString) {
		String[] tokens = contextsString.split(",");
		Set<Context> toret = new HashSet<>();

		for(String token: tokens) {
			try {
				toret.add(Context.valueOf(token));
			} catch(IllegalArgumentException e) {
				throw new IllegalArgumentException("Invalid context parameter, must be a comma-separated list of " +
						"contexts (CG, CHG, CHH)");
			}
		}
		return toret;
	}

	private List<Sample> parseSamples(Project project, String treatmentSamplesString) {
		List<Sample> samplesList = new LinkedList<>();
		if (treatmentSamplesString != null) {
			String[] tokens = treatmentSamplesString.split(",");
			
			for (String token: tokens) {
				boolean found = false;
				for (Sample s: project.getSamples()) {
					if (s.getName().equals(token)) {
						samplesList.add(s);
						found = true;
					}
				}
				if (!found) {
					String validSampleNames = getValidSampleNames(project);
					throw new IllegalArgumentException("sample "+token+" not found. Valid sample names are: "+validSampleNames);
				}
				
			}
		} else {
			throw new IllegalArgumentException("You must give a value for samples parameter");
		}
		return samplesList;
	}

	private String getValidSampleNames(Project project) {
		StringBuilder sb = new StringBuilder();
		for (Sample s: project.getSamples()) {
			sb.append(s.getName()+" ");
		}
		String validSampleNames = sb.toString().trim();
		return validSampleNames;
	}

	@Override
	protected List<Option> createOptions() {
		List<Option> toret = super.createOptions();
		
		toret.add(new Option("treatment-samples", "t", 
				"Comma-separated (with no spaces) list of sample names belonging to 'treatment' group", false, true));
		
		toret.add(new Option("control-samples", "c", 
				"Comma-separated (with no spaces) list of sample names belonging to 'control' group", false, true));

		toret.add(new DefaultValuedOption("context", "x",
				"Comma-separated (with no spaces) list of CpG contexts to analyze: CG, CHG or CHH. For example: CG," +
						"CHG",	"CG"));

		toret.add(new Option("region-beds", "b",
				"Comma-separated (with no spaces) list of BED files to analyze at region-level", true, true));
		
		
		return toret;
	}

}
