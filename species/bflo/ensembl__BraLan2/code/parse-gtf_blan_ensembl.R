# parse gtf file to get protein-ids corresponding to all gene-ids.
# started: 2020-07-10
# description: Proteins FASTA file was used for gene annotations done with eggNOG-mapper, therefore the eggNOG-mapper output files is organized by protein IDs, while the sc-RNA-seq data is organized by gene-ids. Need to make a table with gene-id protein-id associations in order to more easily use the annotations from eggNOG-mapper. 

### packages ==================================================================
library(tidyverse)

### read data =================================================================
gtf <- read.table("species/bflo/genome/ensembl_BraLan2/Branchiostoma_lanceolatum.BraLan2.46.gtf", 
                  header = FALSE, sep = '\t', quote = "", stringsAsFactors = FALSE)

egg <- read.table("species/bflo/genome/ensembl_BraLan2/eggnog/job_MM_ueu6z4wr_annotations.tsv", 
                  header = FALSE, sep = '\t', quote = "", stringsAsFactors = FALSE)

### extract gene and protein ids ==============================================
# extract rows with 'protein_id' attribute
pid <- gtf[grep('protein_id', gtf$V9), ]

# extract gene_id and protein_id
pid$gene_id <- gsub(pattern = '(.*)(gene_id ")(BL\\d*)(")(.*)', 
                    replacement = '\\3', 
                    x = pid$V9)

pid$protein_id <- gsub(pattern = '(.*)(protein_id ")(BL.*)(")(.*)',
                       replacement = '\\3',
                       x = pid$V9)

# clean up
pid <- pid %>% 
  select(gene_id, protein_id) %>%
  unique()

### add eggnog annotation =====================================================
pid$eggnog_protein_name <- egg$V6[match(x = pid$protein_id, table = egg$V1)]

### write =====================================================================
write.csv(pid, 
          "species/bflo/genome/ensembl_BraLan2/blan_ensembl_gene-id_protein-id.csv", 
          row.names = FALSE)

### end =======================================================================
