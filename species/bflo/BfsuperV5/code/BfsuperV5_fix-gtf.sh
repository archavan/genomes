#! /bin/bash

# The Bflo superTranscriptome put together by N. Holland's lab at Oxford was downloaded from https://ora.ox.ac.uk/objects/uuid:51c9ff65-63c4-4010-a267-4b795cc40abc. When I cellranger mkref with fasta and gtf files from this assembly, I got the error that "property transcript_id was not found in the gft in column  9", which in fact is present. Turns out the problem occurs because this gtf is not formatted properly. There needs to be a space between each entry in column 9, but this gtf doesn't have that. See one line below for example. 
# FE_G19231       superTranscript exon    1       510     .       +       .       gene_id "FETU36640";transcript_id "FETU36640";gene_Name "lncRNA";Note "lncRNA predited by CPC2"

# I am using sed's substitute function to edit this gtf to include a space after each occurance of ";" in all rows. We shouldn't expect any semicolons in any column other than 9, so this should only alter the 9th column.   

cd /home/arc78/scratch60/macrophages/bflo/genome/BfsuperV5/files

# replace ; with '; '. Use/g to enable global change. By default sed s only replaces the first occurance on each line. 
sed 's/;/; /g' BfsuperV5-4_FIX.gtf > BfsuperV5-4_FIX_spacing-fixed.gtf


