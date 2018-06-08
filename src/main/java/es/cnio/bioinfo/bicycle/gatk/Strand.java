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

import java.util.regex.Pattern;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

public enum Strand {

	WATSON(".*CT.*WATSON.*"), CRICK(".*CT.*CRICK.*");

	private Pattern fileRegExp;

	private Strand(String suffix) {
		this.fileRegExp = Pattern.compile(suffix);
	}

	public Pattern getFileRegExp() {
		return this.fileRegExp;
	}

	public Context getContext(ReferenceContext ref, long pos) {
		if (this == WATSON) {
			if (ref.getBase() == 'C') {
				if (ref.getWindow().getStop() - pos == 2) {

					byte[] downstream = ref.getBases();
					// correct index if pos <= 2
					int correction = pos == 2? -1:pos == 1? -2: 0;
					if (downstream[3+correction] == 'G') {
						return Context.CG;
					} else if (downstream[4+correction] == 'G') {
						return Context.CHG;
					} else {
						return Context.CHH;
					}
				} else {
					System.err.println("null context in " + ref.getLocus().getContig() + ":" + pos + " strand Watson");
					return null;
				}
			} else {
				throw new IllegalArgumentException("Can't compute Watson context in a reference which is not C");
			}
		} else if (this == CRICK) {
			if (ref.getBase() == 'G') {
			  if (pos <= 2) {
			    System.err.println("Cannot compute context for " + ref.getLocus().getContig() + ":" + pos + " strand Crick");
          return null;
			  }
				if (pos - ref.getWindow().getStart() == 2) {

					byte[] bases = ref.getBases();
					byte[] upstream = new byte[2];
					System.arraycopy(bases, 0, upstream, 0, 2);

					if (upstream[1] == 'C') {
						return Context.CG;
					} else if (upstream[0] == 'C') {
						return Context.CHG;
					} else {
						return Context.CHH;
					}
				} else {
					System.err.println("null context in " + ref.getLocus().getContig() + ":" + pos + " strand Crick");
					return null;
				}
			} else {
				throw new IllegalArgumentException("Can't compute Watson context in a reference which is not C");
			}
		} else {
			throw new IllegalArgumentException("this strand " + this + " does not support context calculation");
		}
	}

	public char getCytosineBase() {
		if (this == WATSON) {
			return 'C';
		}
		if (this == CRICK) {
			return 'G';
		}
		throw new IllegalArgumentException("this strand " + this + " does not support getCytonsine");
	}

	public char getGuanineBase() {
		if (this == WATSON) {
			return 'G';
		}
		if (this == CRICK) {
			return 'C';
		}
		throw new IllegalArgumentException("this strand " + this + " does not support getGuanine");
	}

	public char getThymineBase() {
		if (this == WATSON) {
			return 'T';
		}
		if (this == CRICK) {
			return 'A';
		}
		throw new IllegalArgumentException("this strand " + this + " does not support getThymine");
	}

	public long downstream(long pos) {
		if (this == WATSON) {
			return pos + 1;
		}
		if (this == CRICK) {
			return pos - 1;
		}
		throw new IllegalArgumentException("this strand " + this + " does not support downstream");
	}

	public boolean isNegative() {
		if (this == WATSON) {
			return false;
		}
		if (this == CRICK) {
			return true;
		}
		throw new IllegalArgumentException("this strand " + this + " does not support isNegative");
	}


}
