# 2020-08-07
# This script fixes the GTF file for mitochondrial genome of Branchiostoma floridae. 
# Source of GTF: downloaded as GFF file from NCBI Genbank (https://www.ncbi.nlm.nih.gov/nuccore/AF098298.1). The GFF file was then converted to GTF using gffread v0.12.1, with --force-exons option, which turns the smallest child item of each gene into an exon. 
# Problem with GTF: 24 lines out of 87 lines do not have a gene_id feature in column name. 
# Solution: Add a gene_id feature for these lines, use transcript_id as the gene_id

### packages ==================================================================
library(tidyverse)

### read file =================================================================
gtf <- read.table("species/bflo/genome/genbank/GCA_000003815.1_Version_2_with-mito/AF098298.1_bflo_mito/AF098298.1_bflo_mitochondrion_complete-genome.gtf", 
                  header = FALSE, sep = '\t', quote = "", stringsAsFactors = FALSE)

### fix the file ==============================================================
# fine lines that only have transcript_id and rewrite them to include gene_id
gtf$V9 <- gsub("(transcript_id)( \"rna-AF098298.1:[0-9]*\\.\\.[0-9]*\";)($)", 
               replacement = "\\1\\2 gene_id \\2 \\3", 
               gtf$V9)

### write file ================================================================
write.table(gtf, 
            file = "species/bflo/genome/genbank/GCA_000003815.1_Version_2_with-mito/AF098298.1_bflo_mito/AF098298.1_bflo_mitochondrion_complete-genome_fixed.gtf", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

### end =======================================================================