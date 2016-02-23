package es.cnio.bioinfo.bicycle;

import java.text.DecimalFormat;

public class RegionMethylation {

	private int mCG_WATSON;
	private int CGDepth_WATSON;
	private int mCHG_WATSON;
	private int CHGDepth_WATSON;
	private int mCHH_WATSON;
	private int CHHDepth_WATSON;
	private int mCG_CRICK;
	private int CGDepth_CRICK;
	private int mCHG_CRICK;
	private int CHGDepth_CRICK;
	private int mCHH_CRICK;
	private int CHHDepth_CRICK;
	private String regionName;
	
	public RegionMethylation(
			String regionName,
			int mCG_WATSON, 
			int cGDepth_WATSON, 
			int mCHG_WATSON, 
			int cHGDepth_WATSON, 
			int mCHH_WATSON,
			int cHHDepth_WATSON, 
			
			int mCG_CRICK, 
			int cGDepth_CRICK, 
			int mCHG_CRICK, 
			int cHGDepth_CRICK, 
			int mCHH_CRICK,
			int cHHDepth_CRICK) 
	{
		super();
		this.regionName = regionName;
		this.mCG_WATSON = mCG_WATSON;
		CGDepth_WATSON = cGDepth_WATSON;
		this.mCHG_WATSON = mCHG_WATSON;
		CHGDepth_WATSON = cHGDepth_WATSON;
		this.mCHH_WATSON = mCHH_WATSON;
		CHHDepth_WATSON = cHHDepth_WATSON;
		this.mCG_CRICK = mCG_CRICK;
		CGDepth_CRICK = cGDepth_CRICK;
		this.mCHG_CRICK = mCHG_CRICK;
		CHGDepth_CRICK = cHGDepth_CRICK;
		this.mCHH_CRICK = mCHH_CRICK;
		CHHDepth_CRICK = cHHDepth_CRICK;
	}
	
	public String getRegionName() {
		return regionName;
	}

	public int getmCG_WATSON() {
		return mCG_WATSON;
	}

	public int getCGDepth_WATSON() {
		return CGDepth_WATSON;
	}

	public int getmCHG_WATSON() {
		return mCHG_WATSON;
	}

	public int getCHGDepth_WATSON() {
		return CHGDepth_WATSON;
	}

	public int getmCHH_WATSON() {
		return mCHH_WATSON;
	}

	public int getCHHDepth_WATSON() {
		return CHHDepth_WATSON;
	}

	public int getmCG_CRICK() {
		return mCG_CRICK;
	}

	public int getCGDepth_CRICK() {
		return CGDepth_CRICK;
	}

	public int getmCHG_CRICK() {
		return mCHG_CRICK;
	}

	public int getCHGDepth_CRICK() {
		return CHGDepth_CRICK;
	}

	public int getmCHH_CRICK() {
		return mCHH_CRICK;
	}

	public int getCHHDepth_CRICK() {
		return CHHDepth_CRICK;
	}
	
	public int getmCG() {
		return this.mCG_WATSON + this.mCG_CRICK;
	}
	
	public int getCGDepth() {
		return this.CGDepth_WATSON + this.CGDepth_CRICK;
	}
	
	public int getmCHG() {
		return this.mCHG_WATSON + this.mCHG_CRICK;
	}
	
	public int getCHGDepth() {
		return this.CHGDepth_WATSON + this.CHGDepth_CRICK;
	}
	
	public int getmCHH() {
		return this.mCHH_WATSON + this.mCHH_CRICK;
	}
	
	public int getCHHDepth() {
		return this.CHHDepth_WATSON + this.CHHDepth_CRICK;
	}
	
	public static String getMarshallHeader() {
		final StringBuilder sb = new StringBuilder();
		
		sb.append("#Region\tmCG WATSON\tdepthCG WATSON\tCG WATSON methylation\tmCG CRICK\tdepthCG CRICK\tCG CRICK methylation");
		sb.append("\tmCG total methylation\tmCG total depth\tCG TOTAL methylation");
		sb.append("\tmCHG WATSON\tdepthCHG WATSON\tCHG WATSON methylation\tmCHG CRICK\tdepthCHG CRICK\tCHG CRICK methylation");
		sb.append("\tmCHG total methylation\tmCHG total depth\tCHG TOTAL methylation");
		sb.append("\tmCHH WATSON\tdepthCHH WATSON\tCHH WATSON methylation\tmCHH CRICK\tdepthCHH CRICK\tCHH CRICK methylation");
		sb.append("\tmCHH total methylation\tmCHH total depth\tCHH TOTAL methylation");
		
		return sb.toString();
		
	}
	
	public String marshall() {
		final StringBuilder sb = new StringBuilder(this.regionName);
		final DecimalFormat df = new DecimalFormat("#.#######");
		
		//mCG
		double result=(double)this.getmCG_WATSON()/(double)this.getCGDepth_WATSON();										
		sb.append("\t").append(this.getmCG_WATSON()).append("\t").append(this.getCGDepth_WATSON()).append("\t").append(df.format(result));
		result=(double)this.getmCG_CRICK()/(double)this.getCGDepth_CRICK();										
		sb.append("\t").append(this.getmCG_CRICK()).append("\t").append(this.getCGDepth_CRICK()).append("\t").append(df.format(result));
		result=(double)(this.getmCG_WATSON()+this.getmCG_CRICK())/(double)(this.getCGDepth_WATSON()+this.getCGDepth_CRICK());
		sb.append("\t").append(this.getmCG_WATSON()+this.getmCG_CRICK()).append("\t").append(this.getCGDepth_WATSON()+this.getCGDepth_CRICK()).append("\t").append(df.format(result));

		//mCHG
		result=(double)this.getmCHG_WATSON()/(double)this.getCHGDepth_WATSON();										
		sb.append("\t").append(this.getmCHG_WATSON()).append("\t").append(this.getCHGDepth_WATSON()).append("\t").append(df.format(result));
		result=(double)this.getmCHG_CRICK()/(double)this.getCHGDepth_CRICK();										
		sb.append("\t").append(this.getmCHG_CRICK()).append("\t").append(this.getCHGDepth_CRICK()).append("\t").append(df.format(result));
		result=(double)(this.getmCHG_WATSON()+this.getmCHG_CRICK())/(double)(this.getCHGDepth_WATSON()+this.getCHGDepth_CRICK());
		sb.append("\t").append(this.getmCHG_WATSON()+this.getmCHG_CRICK()).append("\t").append(this.getCHGDepth_WATSON()+this.getCHGDepth_CRICK()).append("\t").append(df.format(result));
		
				
		//mCHH
		result=(double)this.getmCHH_WATSON()/(double)this.getCHHDepth_WATSON();										
		sb.append("\t").append(this.getmCHH_WATSON()).append("\t").append(this.getCHHDepth_WATSON()).append("\t").append(df.format(result));
		result=(double)this.getmCHH_CRICK()/(double)this.getCHHDepth_CRICK();										
		sb.append("\t").append(this.getmCHH_CRICK()).append("\t").append(this.getCHHDepth_CRICK()).append("\t").append(df.format(result));
		result=(double)(this.getmCHH_WATSON()+this.getmCHH_CRICK())/(double)(this.getCHHDepth_WATSON()+this.getCHHDepth_CRICK());
		sb.append("\t").append(this.getmCHH_WATSON()+this.getmCHH_CRICK()).append("\t").append(this.getCHHDepth_WATSON()+this.getCHHDepth_CRICK()).append("\t").append(df.format(result));
		
		return sb.toString();
	}
	
	public static RegionMethylation unmarshall(String line) {
		final String[] tokens = line.split("\t");
		
		return new RegionMethylation(
				tokens[0],
				Integer.parseInt(tokens[1]),
				Integer.parseInt(tokens[2]),
				Integer.parseInt(tokens[10]),
				Integer.parseInt(tokens[11]),
				Integer.parseInt(tokens[19]),
				Integer.parseInt(tokens[20]),
				
				Integer.parseInt(tokens[4]),
				Integer.parseInt(tokens[5]),
				Integer.parseInt(tokens[13]),
				Integer.parseInt(tokens[14]),
				Integer.parseInt(tokens[22]),
				Integer.parseInt(tokens[23])
		);
	}
}
