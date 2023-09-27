In this repository you can find scripts, which were used while creating and using small RNA annotation, published in the article `ITAS: integrated transcript annotation for small RNA` (https://doi.org/10.3390/ncrna8030030). 


Files HUMAN_DATABASES.R, MOUSE_DATABASE.R, RATS_DATABASE.R, CELEGANS_DATABASE.R, DROSOPHILA_DATABASE.R contain commands which were performed to process raw database files, to create consensus database, to search for conflicts in used database for each organism. File exel_to_final_gtf.R contains the script for creating presented gtf-files based on the tables, obtained on the previous step.
File LITERATURE_ANALYSIS.R was used to analyze case studies.


Archived databases gtf-files can be found in the directory https://github.com/EpiEpiMSU/ITAS 


Alignment
=========


#alignment to reference genome by hisat2 


hisat2 -x /path/to/reference -U ./path/to/fastq -S /path/to/sam/output -k 1 --no-spliced-alignment --no-softclip --summary-file /path/to/summary/output


#or by bowtie


bowtie -x /path/to/reference -q ./path/to/fastq -S /path/to/sam/output --large-index -v 1 -m 100 -k 1 --best --strata


#for our databases we used the following versions of reference genomes: 


hg38 for human, mm39 for mouse, rn7 for rat, ce11 for C.elegans and dm6 for drosophila. 


Calculating gene counts
=======================


#gene counts were calculated based on our snRNA databases by function featureCounts from Rsubread library


featureCounts(
files=sample_file, 
annot.ext=database_file, 
isGTFAnnotationFile = TRUE, 
GTF.featureType = "exon", 
GTF.attrType = "transcript_id", 
minFragLength=15)