#' Oncotator annotations for SNVs on chromosome 12 of a PDX sample
#' 
#' @format A data frame with 13055 rows and 6 variables
#' \describe{
#'   \item{gene}{Gene names (HUGO symbol)}
#'   \item{chrom}{The chromosome the variant is on}
#'   \item{loc}{The genomic coordinates of the SNV}
#'   \item{variant_classification}{The oncotator annotation for the type of mutation (i.e. Missense_Mutation)}
#'   \item{ref}{The reference allele}
#'   \item{alt}{The alternate allele}
#' }
"annotations"

#' MuTect1 SNV calls for Chromosome 12 of a PDX sample
#' 
#' @format A data frame with 78505 rows and 6 variables
#' \describe{
#'   \item{chrom}{The chromosome the variant is on}
#'   \item{loc}{The genomic coordinates of the SNV}
#'   \item{ref}{The reference allele}
#'   \item{alt}{The alternate allele}
#'   \item{vaf}{Frequency of the alternate allele}
#'   \item{judgement}{The judgement of MuTect1 on the variant (KEEP/REJECT)}
#' }
"variants"

#' MAPEX scores for the passed variants in \code{variants}
#' 
#' @format MAPEX scores
#' \describe{
#'   \item{chrom}{The chromosome the variant is on}
#'   \item{loc}{The genomic coordinates of the SNV}
#'   \item{variant_score}{The variant score (between 0.0 and 1.0)}
#'   \item{reason}{The most common classification of reads supporting the variant (mouse,on_target,off_target)}
#'   }
"scores"