#' ---
#' title: "Fix GTF file to run with cellranger mkref"
#' author: "Arun Chavan"
#' output: github_document
#' ---
#'
#' Started: 2022-04-29
#'

# background ==================================================================
# We already converted the KH2013 gff3 file to gtf using agat. Here we read in
# the agat-converted gtf file and simplify the attributes column to avoid issues
# with cellranger mkref. We will only keep the gene id and transcript id
# attributes and get rid of everything else.

# setup =======================================================================
library(tidyverse)
library(bioseq)
library(glue)

# functions ===================================================================
extract_attribute <- function(gxf, x) {
  query <- glue::glue(".*({x} \")([^\\s]+?)\".*")
  ifelse(grepl(x, gxf$attribute), 
         gsub(query, "\\2", gxf$attribute), 
         "") 
} 

# data ========================================================================
gxf_headers <- c("seqname", 
                 "source", 
                 "feature", 
                 "start", 
                 "end", 
                 "score", 
                 "strand", 
                 "frame", 
                 "attribute")

gtf <- read_tsv("species/crob/genome/KH13/data/ghost-db/KH.KHGene.2013_agat-converted.gtf",
                comment = "#",
                col_names = gxf_headers,
                col_types = cols())

# extract gene_id and transcript_ids ==========================================
gtf <- gtf %>% 
  mutate(gene_id       = extract_attribute(., "gene_id"),
         transcript_id = extract_attribute(., "transcript_id")) 

# replace attribute column with a simplified one with only gene_id and transcript_id
gtf <- gtf %>%
  mutate(attribute = ifelse(transcript_id == "",
                            glue("gene_id \"{gene_id}\";"),
                            glue("gene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";"))) %>% 
  select(-gene_id, -transcript_id)

# write file ==================================================================
write_tsv(gtf, 
          "species/crob/genome/KH13/data/ghost-db/KH.KHGene.2013_agat-converted_only-gid-tid.gtf", 
          col_names = FALSE)
# end =========================================================================
