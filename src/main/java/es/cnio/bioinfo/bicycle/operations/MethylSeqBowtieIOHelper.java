package es.cnio.bioinfo.bicycle.operations;

import java.io.File;

import es.cnio.bioinfo.bicycle.operations.BowtieAlignment.Quals;
import es.cnio.bioinfo.bicycle.operations.BowtieAlignment.Strand;

public class MethylSeqBowtieIOHelper {

	
	
	
	public String[] getBowtieCommand(String bowtieCommand, boolean isPaired, boolean nohead, 

			Strand strand,
			File ref,
			/* bowtie params */
			final int e,
			final int l,
			final int n,
			
			
			final int chunkmbs,
			final Quals solexaQ,
			
			/*bowtie paired-end parameters*/
			final int I, 
			final int X
			)
	{
		final int M = 1; //tag if there are more than one possible alignment with XM:i:>2
		final int k = 1; //report only 1 alignment
		if (!isPaired){
			if (nohead){
				return new String[]{bowtieCommand, "-t","--chunkmbs",""+chunkmbs,"--mm",solexaQ.getParameterValue(),"-e",""+e,"-l",""+l,"-n",""+n,"-k",""+k,"-M",""+M,"-S","--sam-nohead","--sam-RG","ID:"+strand.name(),"--sam-RG","SM:"+strand.name(),"--best","--nomaqround",ref.getAbsolutePath(),"-"};
			}else{
				return new String[]{bowtieCommand, "-t","--chunkmbs",""+chunkmbs,solexaQ.getParameterValue(),"-e",""+e,"-l",""+l,"-n",""+n,"-k",""+k,"-M",""+M,"-S","--best","--sam-RG","ID:"+strand.name(),"--sam-RG","SM:"+strand.name(),"--nomaqround",ref.getAbsolutePath(),"-"};				
			}
		}else{ //paired
			if (nohead){
				return new String[]{bowtieCommand, "-t","--chunkmbs",""+chunkmbs,"--mm",solexaQ.getParameterValue(),"-e",""+e,"-l",""+l,"-n",""+n,"-k",""+k,"-M",""+M,"-S","--sam-nohead","--sam-RG","ID:"+strand.name(),"--sam-RG","SM:"+strand.name(),"--best","--nomaqround",ref.getAbsolutePath(),"-I", ""+I, "-X",""+X, "--ff", "--12", "-"};
			}else{
				return new String[]{bowtieCommand, "-t","--chunkmbs",""+chunkmbs,solexaQ.getParameterValue(),"-e",""+e,"-l",""+l,"-n",""+n,"-k",""+k,"-M",""+M,"-S","--best","--sam-RG","ID:"+strand.name(),"--sam-RG","SM:"+strand.name(),"--nomaqround",ref.getAbsolutePath(),"-I", ""+I, "-X",""+X, "--ff", "--12", "-"};	
			}
		}
	}
}
