# score reads
.score_reads <- function(read_fr){
  read_len <- max(read_fr$qlen) # assume the reads are same length as the longest query length
  out <- read_fr %>% dplyr::mutate(read_score=dplyr::case_when(.$var_chrom != .$chrom ~ 0, # off target
                                                               .$align_start <= (.$var_start + (read_len - .$qlen)) & # on or near target
                                                                 .$align_end >= (.$var_start - (read_len - .$qlen)) ~ 1,
                                                               TRUE ~ 0),
                                   reason=dplyr::case_when(.$var_chrom != .$chrom & grepl('m',.$chrom) ~ 'mouse',
                                                           .$var_chrom != .$chrom & !grepl('m',.$chrom) ~ 'off_target',
                                                           .$align_start <= (.$var_start + (read_len - .$qlen)) & # on or near target
                                                             .$align_end >= (.$var_start - (read_len - .$qlen)) ~ 'on_target',
                                                           TRUE ~ 'off_target'))
  out
}

.score_variants <- function(scored_read_fr){
  out <- scored_read_fr %>% dplyr::group_by(variant) %>%
    dplyr::summarise(variant_score=mean(read_score),
                     reason=names(which.max(table(reason)))) %>%
    tidyr::separate(variant,c('chrom','loc','stop'),remove=T) %>%
    dplyr::select(chrom,loc,variant_score,reason)
  out
}
