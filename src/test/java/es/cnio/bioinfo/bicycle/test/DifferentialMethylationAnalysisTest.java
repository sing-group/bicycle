package es.cnio.bioinfo.bicycle.test;

import static java.util.Arrays.asList;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import org.easymock.EasyMock;
import org.easymock.EasyMockRunner;
import org.easymock.Mock;
import org.easymock.MockType;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;

import es.cnio.bioinfo.bicycle.Project;
import es.cnio.bioinfo.bicycle.Reference;
import es.cnio.bioinfo.bicycle.Sample;
import es.cnio.bioinfo.bicycle.gatk.Context;
import es.cnio.bioinfo.bicycle.operations.DifferentialMethylationAnalysis;
import es.cnio.bioinfo.bicycle.operations.MethylationAnalysis;

@RunWith(EasyMockRunner.class)
public class DifferentialMethylationAnalysisTest {

	private DifferentialMethylationAnalysis dma;

	
	@Mock
	private Project project;
	
	@Mock
	private MethylationAnalysis ma;
	
	@Mock
	private Reference reference;
	
	@Mock
	private Sample controlSample1;

	@Mock
	private Sample controlSample2;
	
	@Mock
	private Sample treatmentSample1;
	
	@Mock
	private Sample treatmentSample2;
	
	private File controlSample1File = new File(Utils.getMethylcytosinesDirectory()+File.separator+"controlSample1.methylcytosines");
	private File controlSample2File = new File(Utils.getMethylcytosinesDirectory()+File.separator+"controlSample2.methylcytosines");
	private File treatmentSample1File = new File(Utils.getMethylcytosinesDirectory()+File.separator+"treatmentSample1.methylcytosines");
	private File treatmentSample2File = new File(Utils.getMethylcytosinesDirectory()+File.separator+"treatmentSample2.methylcytosines");
	
	private File regionsFile = new File(Utils.getBedsDirectory()+File.separator+"regions.bed");
	
	@Before
	public void instantiateDMA() {
		this.dma = new DifferentialMethylationAnalysis(this.ma, new HashSet<>(Arrays.asList(Context.CG)));
	}
	
	@Test
	public void simpleTest() throws IOException {
		
		File tmpDir = new File(System.getProperty("java.io.tmpdir"));
		
		EasyMock.expect(ma.getMethylcytosinesFile(reference, controlSample1)).andReturn(controlSample1File);
		EasyMock.expect(ma.getMethylcytosinesFile(reference, controlSample2)).andReturn(controlSample2File);
		EasyMock.expect(ma.getMethylcytosinesFile(reference, treatmentSample1)).andReturn(treatmentSample1File);
		EasyMock.expect(ma.getMethylcytosinesFile(reference, treatmentSample2)).andReturn(treatmentSample2File);
		
		EasyMock.expect(reference.getSequenceNames())
				.andReturn(asList("control", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
						"chr18", "chr19", "chr1", "chr20", "chr21", "chr22", "chr2", "chr3", "chr4", "chr5", "chr6",
						"chr7", "chr8", "chr9", "chrM", "chrX", "chrY"));
		EasyMock.expect(reference.getReferenceFile()).andReturn(new File("hg18.fa")).anyTimes();
		
		EasyMock.expect(controlSample1.getName()).andReturn("C1").anyTimes();
		EasyMock.expect(controlSample2.getName()).andReturn("C2").anyTimes();
		EasyMock.expect(treatmentSample1.getName()).andReturn("T1").anyTimes();
		EasyMock.expect(treatmentSample2.getName()).andReturn("T2").anyTimes();
	
		EasyMock.expect(ma.getProject()).andReturn(project).anyTimes();
		
		EasyMock.expect(project.getOutputDirectory()).andReturn(tmpDir).anyTimes();
		
		EasyMock.replay(ma, project);
		EasyMock.replay(reference);
		EasyMock.replay(controlSample1, controlSample2, treatmentSample1, treatmentSample2);
		
		
		
		List<Sample> controlSamples = asList(controlSample1, controlSample2);
		List<Sample> treatmentSamples = asList(treatmentSample1, treatmentSample2);
		
		dma.analyzeDifferentialMethylationByBase(reference, treatmentSamples, controlSamples);
		
		Assert.assertTrue(dma.getDifferentiallyMethylatedCytosinesFile(reference, treatmentSamples, controlSamples).exists());
		
		EasyMock.verify(ma);
		
	}
	@Test
	public void regionsTest() throws IOException {
		
		File tmpDir = new File(System.getProperty("java.io.tmpdir"));
		
		EasyMock.expect(ma.getMethylcytosinesFile(reference, controlSample1)).andReturn(controlSample1File);
		EasyMock.expect(ma.getMethylcytosinesFile(reference, controlSample2)).andReturn(controlSample2File);
		EasyMock.expect(ma.getMethylcytosinesFile(reference, treatmentSample1)).andReturn(treatmentSample1File);
		EasyMock.expect(ma.getMethylcytosinesFile(reference, treatmentSample2)).andReturn(treatmentSample2File);
		
		EasyMock.expect(reference.getSequenceNames())
		.andReturn(asList("control", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
				"chr18", "chr19", "chr1", "chr20", "chr21", "chr22", "chr2", "chr3", "chr4", "chr5", "chr6",
				"chr7", "chr8", "chr9", "chrM", "chrX", "chrY"));
		EasyMock.expect(reference.getReferenceFile()).andReturn(new File("hg18.fa")).anyTimes();
		
		EasyMock.expect(controlSample1.getName()).andReturn("C1").anyTimes();
		EasyMock.expect(controlSample2.getName()).andReturn("C2").anyTimes();
		EasyMock.expect(treatmentSample1.getName()).andReturn("T1").anyTimes();
		EasyMock.expect(treatmentSample2.getName()).andReturn("T2").anyTimes();
		
		EasyMock.expect(ma.getProject()).andReturn(project).anyTimes();
		
		EasyMock.expect(project.getOutputDirectory()).andReturn(tmpDir).anyTimes();
		
		EasyMock.replay(ma, project);
		EasyMock.replay(reference);
		EasyMock.replay(controlSample1, controlSample2, treatmentSample1, treatmentSample2);
		
		
		
		List<Sample> controlSamples = asList(controlSample1, controlSample2);
		List<Sample> treatmentSamples = asList(treatmentSample1, treatmentSample2);
		
		dma.analyzeDifferentialMethylationByRegions(reference, treatmentSamples, controlSamples, regionsFile);
		
		Assert.assertTrue(dma.getDifferentiallyMethylatedRegionsFile
				(reference, treatmentSamples, controlSamples, regionsFile).exists());

		EasyMock.verify(ma);
	}
}
