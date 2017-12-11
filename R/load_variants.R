#' Load a variant call file.
#' Allowed types are either maflite, vcf (version >=4), call_stats (MuTect1 format),
#' or snp (Varscan2 format) files
#'
#' The function is generally internal, but is exported to help the user in debugging input
#'
#' @param variant_file /path/to/variant/file
#' @param file_type either maf, vcf, call_stats, or snp
#' @return string of the form chrom:startposition-startposition_ref_alt
#'     used to identify variants when extracting reads from BAM files.
#'
#' @examples
#' # An example of converting oncotator output to maflite and calling load_variants
#' \dontrun{
#'
#'
#' }
.load_variants <- function(variant_file,file_type='maf'){
  ftypes <- c('vcf','call_stats','snp','maf')
  ftype <- match(file_type,ftypes)

  if (! file.exists(variant_file)){
    stop('Variant file does not exist - Please check the path')
  }
  if (ftype==1){
    vars <- .get_variant_locs_vcf(variant_file)
  }
  if (ftype==2){
    vars <- .get_variant_locs_cs(variant_file)
  }
  if (ftype==3){
    vars <- .get_variant_locs_snp(variant_file)
  }
  if (ftype==4) {
    vars <- .get_variant_locs_maf(variant_file)
  }
  vars
}

# Return variant locations from maflite
.get_variant_locs_maf <- function(variant_file){
  cols <- c('chr','start','end','ref_allele','alt_allele')
  format <- c('ccccc')
  readr::read_tsv(variant_file,comment='#',col_names=cols,col_types=format) %>%
    dplyr::mutate(loc_string=paste0(chr,':',start,'-',start,'_',ref_allele,'_',alt_allele)) %>%
    .pull(loc_string)
}
# Return variant locations from vcf
.get_variant_locs_vcf <- function(variant_file){
  cols <- c('chr','pos','id','ref','alt','qual','filter','info','format','tumor','normal')
  format <- c('ccccccccccc')
  readr::read_tsv(variant_file,comment="#",col_names=cols,col_types=format) %>%
    dplyr::filter(filter=='PASS') %>%
    dplyr::mutate(loc_string=paste0(chr,':',pos,'-',pos,'_',ref,'_',alt)) %>%
    .pull(loc_string)
}

# Return variant locations from Mutect1 call_stats file
.get_variant_locs_cs <- function(variant_file){
  col_spec <- readr::cols_only(contig='c',position='c',ref_allele='c',alt_allele='c',judgement='c')
  readr::read_tsv(variant_file,comment='#',col_types=col_spec) %>%
    dplyr::filter(judgement=='KEEP') %>%
    dplyr::mutate(loc_string=paste0(contig,':',position,'-',position,'_',ref_allele,'_',alt_allele)) %>%
    .pull(loc_string)
}

# Return variant locations from Varscan .snp file
.get_variant_locs_snp <- function(variant_file){
  # Varscan somatic leaves GL chromosome variants in. They are removed here.
  message('High confidence calls are generated with the default varscan algorithm, ',
          'tumor VAF > 15% normal VAF < 5%, p-val < 0.03')
  col_spec <- readr::cols_only(chrom='c',position='c',ref='c',var='c',
                               somatic_status='c',normal_var_freq='c',tumor_var_freq='c',
                               somatic_p_value='d')
  readr::read_tsv(variant_file,comment='#',col_types=col_spec) %>%
    dplyr::filter(! stringr::str_detect(chrom,'GL')) %>%
    dplyr::mutate(normal_var_freq=readr::parse_number(normal_var_freq),
                  tumor_var_freq=readr::parse_number(tumor_var_freq)) %>%
    dplyr::filter(somatic_status=="Somatic" & normal_var_freq < 5 & tumor_var_freq > 15 & somatic_p_value < 0.03) %>%
    dplyr::mutate(loc_string=paste0(chrom,':',position,'-',position,'_',ref,'_',var)) %>%
    .pull(loc_string)
}
