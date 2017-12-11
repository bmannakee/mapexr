
.generate_blast_input <- function(variant_loc,bam,bam_idx,min_mapq){
  is_insertion <- F # flag to use to determine how to extract variant reads
  is_deletion <- F
  tmp <- stringr::str_split(variant_loc,'_',simplify=T)
  # for indels we want to pile and extract reads at the variant site + 1
  # but we need to keep the NAME of the variant the same so that it will match
  # up with the chromosome and position we have in the vcfs. So we have two locations.
  # One for the pileup, and one to correctly name the reads for the variants we are considering.
  variant_loc <- tmp[[1]]
  loc <- tmp[[1]]
  ref_allele <- tmp[[2]]
  alt_allele <- tmp[[3]]
  if (nchar(ref_allele) > nchar(alt_allele)){
    # we have a deleletion. skip to next position
    is_deletion <- T
  }
  if (nchar(alt_allele) > nchar(ref_allele)){
    # we have an insertion. We set a flag, and call differently
    is_insertion <- T
  }

  pos <- GenomicRanges::GRanges(loc)
  mybamflag <- Rsamtools::scanBamFlag(isNotPassingQualityControls=NA,
                           isDuplicate=FALSE,
                           hasUnmappedMate=NA)
  bamparams <- Rsamtools::ScanBamParam(flag=mybamflag,mapqFilter = min_mapq, simpleCigar=F, tag=c('RG','HC'),which=pos,what='seq')
  if (is_insertion){
    # Insertions and deletions are hard because the insertion can be the same as the actual next sequence
    # The only way to be sure you are getting the insertion reads is to look for I or D in
    # the cigar string for insertions or deletions respectively
    gal <- GenomicAlignments::readGAlignments(bam,index=bam_idx,param=bamparams,use.names=T)
    alt_gal <- gal[which(stringr::str_detect(GenomicAlignments::cigar(gal),'I'))]
    alt_read_names <- names(alt_gal)
    alt_reads <- GenomicRanges::mcols(alt_gal)$seq
    # In at least some samples the read groups are not unique. This gaurantees unique names
    alt_read_names <- paste(alt_read_names,seq(1,length(alt_read_names)),sep=':')
  } else if (is_deletion){
    gal <- GenomicAlignments::readGAlignments(bam,index=bam_idx,param=bamparams,use.names=T)
    alt_gal <- gal[which(stringr::str_detect(GenomicAlignments::cigar(gal),'D'))]
    alt_read_names <- names(alt_gal)
    alt_reads <- GenomicRanges::mcols(alt_gal)$seq
    # In at least some samples the read groups are not unique. This gaurantees unique names
    alt_read_names <- paste(alt_read_names,seq(1,length(alt_read_names)),sep=':')
  } else{
    gal <- GenomicAlignments::stackStringsFromBam(bam,bam_idx,what='seq',param=bamparams,use.names=T)
    #seqlevels(gal) <- seqlevelsInUse(pos)
    # get the piled nucleotides at pos.
    nucs <- as.data.frame(gal)$x # '-' for deletion, the second letter of alt_allele for insertion
    # limit gal to only those rows with a variant read
    alt_gal <- gal[which(nucs==alt_allele)]
    alt_read_names <- names(alt_gal)
    alt_reads <- GenomicRanges::mcols(alt_gal)$seq
    # In at least some samples the read groups are not unique. This gaurantees unique names
    alt_read_names <- paste(alt_read_names,seq(1,length(alt_read_names)),sep=':')
  }

  if (length(alt_reads)==0){
    # sometimes, particularly for Indels, MuTect2's local realignment
    # seems to create them, when they are not visible in the original BAM file.
    # There is nothing we can do here.
    warning(paste0('No variant reads for ',variant_loc))
    return(c('\n'))
  }
  alt_read_names_line <- paste0('>',variant_loc,'_',alt_read_names) # adding pos ensures uniqueness if read covers more than one variant
  blast_input <- paste(alt_read_names_line,alt_reads,sep='\n')
  # blast_out_file <- paste0(blast_output_path,'/',sample_id,'_blast_input.txt')
  # read_list_out_file <- paste0(blast_output_path,'/',sample_id,'_read_list.txt')
  # lapply(paste(alt_read_names_line,sep='\n'),write,read_list_out_file,append=T,ncolumns=1)
  # lapply(blast_input,write,blast_out_file,append=T,ncolumns=1)
  # loc
  blast_input
}

