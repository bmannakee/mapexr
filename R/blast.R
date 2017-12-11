#' Run local blast on prepared sequenc.
#' @param variant_locs list of variant locations
#' @param blast_out
#' @param bam /path/to/tumor/bam
#' @param blast_path /path/to/blast/executable
#' @param dp_path /path/to/combined/genome/blastdb
#' @param blast_threads number of thread to run blast with

#' @return list containing best hit blast scores


.run_blast <- function(variant_locs,blast_out,bam_path,bam_idx_path,blast_path,db_path,blast_threads,min_mapq){
  blast_in <- tempfile(pattern='blast_in',tmpdir=tempdir())
  blast_input <- lapply(variant_locs,.generate_blast_input,bam=bam_path,bam_idx=bam_idx_path,min_mapq=min_mapq)
  if (is.null(blast_out)) {
    blast_out <- tempfile(pattern='blast_out',tmpdir=tempdir())
  }
  invisible(lapply(blast_input,write,blast_in,append=T,ncolumns=1))
  cmd <- paste0(blast_path, ' -task megablast', ' -query ', blast_in, ' -num_threads ', blast_threads,
                ' -db ', db_path, ' -max_target_seqs 1',
                ' -outfmt 6', ' -out ', blast_out)
  system(cmd)
  fr <- .load_blast(blast_out)
  fr
}
