# bicycle
bicycle (bisulfite-based methylcytosine caller) is a next-generation sequencing bioinformatics pipeline able to perform a full DNA methylation level analysis. More info at the [bicycle](http://sing.ei.uvigo.es/bicycle) project page.

## How it works?
**bicycle** (_**bi**sulfite-based methyl**cy**tosine **c**al**le**r_) is a next-generation sequencing bioinformatic pipeline aimed to analyze whole genome bisulfite sequencing data. It can process data from directional (Lister) and non-directional (Cokus) bisulfite sequencing protocols, and from single-end and paired-end sequencing, and performs methylation calls for cytosines in CG and non-CG contexts (CHG and CHH).

**bicycle** uses as input the bisulfite sequencing files from the different samples (FASTQ format) and a reference genome (FASTA format). It then performs: **generation and indexing of Watson and Crick bisulfited versions of the reference genome**, **in-silico bisulfitation of sequenced reads**, **read alignment**, **error estimation in bisulfite conversion**, **identification of clonal and ambiguous reads**, **cytosine methylation detection in CG and non-CG contexts, with non-CG to CG context correction when appropriated**, **calculates methylation ratios, beta scores and weighted mean of cytosine methylation status**, and performs **genomic annotation of methylated regions**, and **differential methylation for cytosines (DMC) and genomic regions (DMR)**.

## Requirements
Bicycle requirements:
- Operating System: Linux/OSX
- Java 1.8+
- Bowtie aligner 0.12.7+
- samtools 0.1.8+

## Coordinators
- Florentino Fdez-Riverola ([SING Research Group](http://sing.ei.uvigo.es/). Informatics Department. University of Vigo)
- David Pisano ([Bioinformatics Unit](http://www.cnio.es/es/grupos/plantillas/presentacion.asp?grupo=50004300). CNIO: Spanish National Cancer Research Centre)


## Authors
  * Daniel Glez-Peña (SING Research Group. Informatics Department. University of Vigo)
  * Osvaldo Graña (Bioinformatics Unit. CNIO: Spanish National Cancer Research Centre)
  * Hugo López-Fernández (SING Research Group. Informatics Department. University of Vigo)

## License
  * GNU Lesser General Public License v3.0

