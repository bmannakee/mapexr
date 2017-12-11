#' function to load blast results in format #6 i.e. -outfmt 6
#' @param blast_output /path/to/blast/output
#' @return tibble with named columns

.load_blast <- function(blast_output){
  col_names <- c('query_id','chrom','percent_id','qlen','mismatches','gaps',
                 'qstart','qend','align_start','align_end','eval','bitscore')
  col_spec <- readr::cols('query_id'='c','chrom'='c','percent_id'='d','qlen'='i',
                   'mismatches'='c','gaps'='c','qstart'='c','qend'='c',
                   'align_start'='i','align_end'='i','eval'='d',
                   'bitscore'='c')
  fr <- readr::read_tsv(blast_output,col_names=col_names,col_types=col_spec)
  fr <- fr %>% dplyr::select(query_id,chrom,percent_id,qlen,align_start,align_end,eval,bitscore) %>%
    tidyr::separate(query_id,into=c('variant','read_id'),sep='_',remove=T) %>%
    tidyr::separate(variant,into=c('var_chrom','var_start','var_end'),convert=T,remove=F)
  # Blast often returns more than one best hit per read. they come back ordered best to worst. get the best one
  fr <- fr %>% dplyr::group_by(variant,read_id) %>% dplyr::filter(bitscore==max(bitscore))
  fr <- fr %>% dplyr::ungroup()
  fr
}
