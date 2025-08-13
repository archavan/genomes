# parse gtf file to create a dataframe with:
# 1. gene id
# 2. protein id
# 3. genomic scaffold (to get MT genes)
# 4. gene name (from eggnog OG groups)
# started: 2020-08-14

### packages ==================================================================
library(tidyverse)

### read data =================================================================
gtf <- read.table("species/bflo/genome/genbank/GCA_000003815.1_Version_2_with-mito/appended/GCA_000003815.1_Version_2_genomic_with-mito.gtf", 
                  header = FALSE, sep = '\t', quote = "", stringsAsFactors = FALSE)

egg <- read.table("species/bflo/genome/orthologs/eggnog/genbank/job_MM_r6ttp7pa_annotations.tsv", 
                  header = FALSE, sep = '\t', quote = "", stringsAsFactors = FALSE) # this does not include the mitochondrial genes, but that's fine. The mitochodnrial genes are only about 10 and their names are in the gtf. 

egg$V6[egg$V6 == ""] <- NA # replace no gene names with NA

### extract attributes from gtf and add them as columns =======================
# gene_id
gtf$gene_id <- gsub(pattern = '(.*)(gene_id\\s+")([^;"]+)(";.*)',
                    replacement = '\\3',
                    gtf$V9)

# gene_name
# not all rows have a gene_name attribute. when no match is found gsub return the entire input string. We'd rather have an empty string or NA. To get around, use grepl to get indices where a match is found and replace only for those rows.
gtf$gene_name <- ifelse(
  grepl(pattern = '(.*)(gene_name\\s+")([^;"]+)(";.*)', gtf$V9), 
  yes = gsub(pattern = '(.*)(gene_name\\s+")([^;"]+)(";.*)',
             replacement = '\\3',
             gtf$V9), 
  no = NA
)

# protein_id
gtf$protein_id <- ifelse(
  grepl(pattern = '(.*)(protein_id\\s+")([^;"]+)(";.*)', gtf$V9), 
  yes = gsub(pattern = '(.*)(protein_id\\s+")([^;"]+)(";.*)',
             replacement = '\\3',
             gtf$V9), 
  no = NA
)

colnames(gtf)[1] <- "location"
### create separate dfs to organized attributes ===============================
## pull out gene_id and gene_name with location and get unique rows -----------
gid.gname <- gtf[, c("location", "gene_id", "gene_name")] %>% unique()

## pull out gene_id and protein id and get unique rows ------------------------
gid.pid <- gtf[, c("gene_id", "protein_id")]
gid.pid <- gid.pid[complete.cases(gid.pid), ] %>% unique()

## add attributes to gid.gname from gid.pid and from egg ----------------------
# add protein_ids
gid.gname$protein_id <- gid.pid$protein_id[match(x = gid.gname$gene_id,
                                                 table = gid.pid$gene_id)]

# gene name from eggnog output
gid.gname$eggnog_name <- egg$V6[match(x = gid.gname$protein_id, table = egg$V1)]

# add eggnog name values to gene name column. Mito genes already have names, avoid those rows by filtering on NA
i <- which(is.na(gid.gname$gene_name) & !is.na(gid.gname$eggnog_name))
gid.gname$gene_name[i] <- gid.gname$eggnog_name[i]

# where gene name is NA, replace with gene id
j <- which(is.na(gid.gname$gene_name))
gid.gname$gene_name[j] <- gid.gname$gene_id[j]

gid.gname$eggnog_name <- NULL # not needed anymore

# add gene description from eggnog output
gid.gname$gene_description <- egg$V22[match(x = gid.gname$protein_id,
                                            table = egg$V1)]

# add GO terms from eggnog output
gid.gname$GO_terms <- egg$V7[match(x = gid.gname$protein_id,
                                   table = egg$V1)]

### write =====================================================================
write.csv(gid.gname, 
          "species/bflo/genome/genbank/GCA_000003815.1_Version_2_with-mito/appended/gene-info.csv", 
          row.names = FALSE)
### end =======================================================================
