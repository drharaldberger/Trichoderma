# Trichoderma
Tools for analysing transcriptome and ChIP data

The input for these scripts are coverage traces from RNA-seq experiments. Typically these can be calculated from sorted and indexed .bam files as results from read mappings. These traces can be produced using e.g. bedtools like:
bedtools genomecov -d -scale <normalization scale> -strand <+/->  -ibam input.bam  > output_trace.sgr

This produced a file with the coverage per base pair of the chromosomes. If strand specific RNAseq has been performed two .sgr files can be calculated what allows better separation between gene specific transcripts.

Multiple traces can be summed up using sum_up_sgr.R for better estimations of Transcript Start or End Points (TSP or TEP), what may be important especially for low abundant transcripts.
Usage of sum_up_sgr.R:

Rscript sum_up_sgr.R /input/folder/with_sgr_files output.rda

The sgr files to process need to be in one folder and have the ending .sgr!
The output will be an R-data file, which is much smaller then the sgr file.

In the next step the TSPs or TEPs are estimated based on the idea that the highest slopes of transcript coverage traces upstream of the 5' ATG may indicate possible TSP or downstream of the 3' STOP codon may indicate TEPs.

Rscript TSP_from_transcript.R input.rda Annotation_file.gff3 Remove_in_gene_name Threshold output.bed
Rscript TEP_from_transcript.R input.rda Annotation_file.gff3 Remove_in_gene_name Threshold output.bed

The input file is the output from sum_up_sgr.R, the annotation file needs to have 'exon' features. Column 9 has only the gene name and something like ID= or Parent= what can be removed using Remove_in_gene_name.
A Threshold can be provided, by using 0 it will be calculated internally.

output is a bed file with a score (0-1000), strand and name info.

