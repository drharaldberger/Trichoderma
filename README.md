# Trichoderma
<h3>Theory:</h3>
These scripts are based on the ideas used in mountainClimber (mountainClimber Identifies Alternative Transcription Start and Polyadenylation Sites in RNA-Seq) by Ashley A Cass and Xinshu Xiao 2019
<a href="https://github.com/gxiaolab/mountainClimber">github.com/gxiaolab/mountainClimber</a> <br>
but has been modified for gene dense genomes as common in fungi.
Basically the stonges upstream slopes or downstream slopes are considered as TSP or TEP, respectively.<br><br>
<h3>Walkthrough:</h3>
The input for these scripts are coverage traces from RNA-seq experiments. Typically these can be calculated from sorted and indexed .bam files as results from read mappings. These traces can be produced using e.g. bedtools like:

<tt>bedtools genomecov -d -scale <normalization scale> -strand <+/->  -ibam input.bam  > output_trace.sgr</tt>

This produced a file with the coverage per base pair of the chromosomes. If strand specific RNAseq has been performed two .sgr files can be calculated what allows better separation between gene specific transcripts.

Multiple traces can be summed up using sum_up_sgr.R for better estimations of Transcript Start or End Points (TSP or TEP), what may be important especially for low abundant transcripts.
Usage of sum_up_sgr.R:

<tt>Rscript sum_up_sgr.R /input/folder/with_sgr_files output.rda</tt>

The sgr files to process need to be in one folder and have the ending .sgr!
The output will be an R-data file, which is much smaller then the sgr file.

In the next step the TSPs or TEPs are estimated based on the idea that the highest slopes of transcript coverage traces upstream of the 5' ATG may indicate possible TSP or downstream of the 3' STOP codon may indicate TEPs.

<tt>Rscript TSP_from_transcript.R input.rda Annotation_file.gff3 Remove_in_gene_name Threshold output.bed</tt>
<br> or <br>
<tt>Rscript TEP_from_transcript.R input.rda Annotation_file.gff3 Remove_in_gene_name Threshold output.bed</tt>

The input file is the output from sum_up_sgr.R, the annotation file needs to have 'exon' features. Column 9 has only the gene name and something like ID= or Parent= what can be removed using Remove_in_gene_name. 
A line in the gff file should look like:

<tt>unitig_0_consensus	maker	exon	7295	9560	.	+	.	Parent=mRNA_TrA0001W</tt>

A Threshold can be provided, by using 0 it will be calculated internally.

Requirements for these scripts are the packages GenomicRanges and stringr which need to be in libpath(). This is sometimes a mess in Windows.

output is a bed file with a score (0-1000), strand and name info.

