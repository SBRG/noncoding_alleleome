# E. coli Non-Coding Alleleome

This repository contains code and data for generating and analyzing the E. coli non-coding alleleome and associated figures. This alleleome contains 1,169 reference non-coding regions from strain K-12 MG1655 that are well-annotated with non-coding functional sites such as transcription factor binding sites, core promoters, and transcriptional attenuators from [RegulonDB](http://regulondb.ccg.unam.mx/). These regions are mapped across 2,350 E. coli strains downloaded from [BV-BRC](https://www.bv-brc.org/) (formerly known as PATRIC).

## Abstract

Microbial genome sequences are rapidly accumulating, enabling large-scale studies of sequence variation. Existing studies primarily focus on coding regions to study patterns of amino acid substitutions in proteins. However, non-coding regulatory regions also play a distinct role in determining physiologic responses. To investigate intergenic sequence variation on a large-scale, we identified non-coding regulatory region alleles across 2,350 fully-sequenced Escherichia coli strains. This “alleleome” consists of 117,781 unique alleles for 1,169 reference regulatory regions (transcribing 1,975 genes) at single base-pair resolution. We find a high degree of conservation in non-coding sequences; overall 64% of nucleotide positions are invariant, and variant positions vary in a median of just 0.6% of strains. Additionally, non-coding alleles are sufficient to recover E. coli phylogroups. We find that core promoter elements and transcription factor binding sites are conserved, especially those located upstream of essential or highly-expressed genes. However, variability in conservation of transcription factor binding sites is significant both within and across regulons. Finally, we contrast mutations acquired during adaptive laboratory evolution with wild-type variation, finding that the former preferentially alter positions that the latter conserves. Overall, this analysis elucidates the wealth of information found in E. coli non-coding sequence variation and expands pangenomic studies to non-coding regions at single-nucleotide resolution.

## Required Software

The following software are required to successfully run this workflow and analysis:

- Python 3.10+
- Python packages (can be installed with `pip <package name>`):
	- biopython
	- logomaker
	- matplotlib
	- numpy
	- pandas
	- pymodulon
	- scipy
	- sklearn 
	- statsmodels
	- seaborn
- Command line programs
	- BLAST+ (https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
	- CD-HIT (https://github.com/weizhongli/cdhit)
	- MUSCLE (https://github.com/rcedgar/muscle)
	- ClermonTyping (https://github.com/A-BN/ClermonTyping#dependencies--installation)

## Overview

The processing and analysis workflow are organized into 2 Jupyter notebooks: `0__build_alleleome` and `1__analyze_alleleome`. Only the latter need be run to analyze the pre-processed and analyzed alleleome reported in this publication. Otherwise, the former notebook may be used to re-generate an alleleome with more genomes/reference regions, or for a different organism. See these notebooks for detailed instructions on execution and analysis.
