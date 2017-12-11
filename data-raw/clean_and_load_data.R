# script to generate raw data for vignettes
# Chrom 12 only for space
library(tidyverse)
# Load mutect1 variant call_stats
load_call_stats <- function(variant_file) {
  col_spec <- readr::cols_only(contig='c',position='i',ref_allele='c',alt_allele='c',judgement='c',tumor_f='d')
  readr::read_tsv(variant_file,comment='#',col_types=col_spec) %>%
    dplyr::rename('chrom'=contig,'loc'=position,'ref'=ref_allele,'alt'=alt_allele,'vaf'=tumor_f) %>%
    dplyr::filter(chrom=="12")
}

# Load annotations from cleaned oncotator file
load_annotations <- function(annotation_file){
  fr <- readr::read_tsv(annotation_file,col_types='cccccicccccccccccccc',comment='#',col_names=F)
  fr <- fr %>% dplyr::rename(gene=X1,chrom=X5,loc=X6,
                             variant_classification=X9,ref=X11,alt=X13) %>%
    select(gene,chrom,loc,variant_classification,ref,alt)
  fr %>% dplyr::filter(chrom=="12")
}

# Load variant scores from mapexr output Only chromosome 12 for space
load_scores <- function(score_file){
  read_tsv(score_file,col_types='cidc') %>% dplyr::filter(chrom=="12")
}

variants <- load_call_stats('./data-raw/PH12_29x1a1.call_stats.chr12.txt')
annotations <- load_annotations('./data-raw/PH12_29x1a1.oncotated.cleaned.txt')
scores <- load_scores('./data-raw/PH12_29x1a1.mapexr.scored.tsv')

# merge data
merged <- variants %>% left_join(annotations,by=c('chrom','loc','ref','alt')) %>%
  left_join(scores,by=c('chrom','loc'))

devtools::use_data(variants,annotations,scores,overwrite=T)
