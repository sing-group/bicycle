/*

Copyright 2012 Daniel Gonzalez Peña, Osvaldo Graña


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

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;


@Reference(window=@Window(start=-2,stop=2))
public class ComputeErrorFromContig extends LocusWalker<DefaultBisulfiteError, DefaultContigBisulfiteError>
		implements TreeReducible<DefaultContigBisulfiteError> {
	@Argument(doc="removeClones", required=false)
    public boolean removeClones=false;
	
	public Tools tools;
	private DefaultContigBisulfiteError lastTraversalResult;
	
	@Override
	public DefaultBisulfiteError map(RefMetaDataTracker arg0, ReferenceContext refContext, AlignmentContext alignmentContext) {
		//System.out.println(this.removeClones);
		if (refContext.getBase()==Strand.WATSON.getCytosineBase()){
			
			return computeError(refContext, alignmentContext, Strand.WATSON);
		}else if (refContext.getBase()==Strand.CRICK.getCytosineBase()){	
			return computeError(refContext, alignmentContext, Strand.CRICK);
		}

		return null;
	}
	private DefaultBisulfiteError computeError(ReferenceContext arg1,
			AlignmentContext arg2, Strand strand) {
		
		ReadBackedPileup strandReads = ListerFilter.applyFilters(tools.getReadsForStrand(strand, arg2, removeClones));
		
		Context context = strand.getContext(arg1, arg2.getPosition());
		if (strandReads!=null && context != null){
			
			DefaultBisulfiteError toret = new DefaultBisulfiteError(context, strand);
			
			int error=0;
			for (byte base: strandReads.getBases()){
				if (base==strand.getCytosineBase()){
					error++;
				}
			}
			
			toret.add(strandReads.getBases().length,error);	
			return toret;
		}
		return null;
	}

	
	
	
	@Override
	public DefaultContigBisulfiteError reduce(DefaultBisulfiteError arg0, DefaultContigBisulfiteError arg1) {
		if (arg0==null){
			return arg1;
		}
		
		arg1.addError(arg0.getStrand(), arg0.getContext(), arg0.getTotalReads(), arg0.getErrorReads());
		return arg1;
	}

	@Override
	public DefaultContigBisulfiteError reduceInit() {
		return new DefaultContigBisulfiteError();
	}

	@Override
	public DefaultContigBisulfiteError treeReduce(DefaultContigBisulfiteError left, DefaultContigBisulfiteError right) {
	
		for (Strand strand:Strand.values()){
			
			for (Context context: Context.values()){
				left.addError(strand, context, right.getError(strand, context).getTotalReads(), right.getError(strand, context).getErrorReads());
			}
		}
		
		return left;

	}
	
	
	
	@Override
	public void onTraversalDone(DefaultContigBisulfiteError result) {
		this.lastTraversalResult = result;
	}
	public DefaultContigBisulfiteError getLastTraversalResult() {
		return lastTraversalResult;
	}
	
	

}
