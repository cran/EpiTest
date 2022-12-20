
#' Phenotypes for the American maize NAM dataset
#'
#' \code{Pheno.list} is a list of three data.frames corresponding to three
#' populations of the American maize NAM experiment.
#' Each data.frame contains the phenotypic values of 191 to 196 recombinant
#' inbred lines (RILs) evaluated for five different traits, along with
#' two extra columns:
#' Genotype (providing the RIL ID) and Weight (providing the weight associated
#' with the RIL, see function EpiTest.fit for details).
#'
#' @format A list of three data.frames
#'
"Pheno.list"

#' Allele ancestries for the American maize NAM dataset
#'
#' \code{Ancestry.list} is a list of numeric matrices corresponding to three populations
#' of the American maize NAM experiment. Each matrix has 0/1 entries representing
#' the ancestry of alleles (0: homozygous for B73 alleles and 1: homozygous
#' for alternative parent alleles).
#' Each matrix includes 191 to 196 recombinant inbred lines characterized at 1,106 markers.
#' @format A list of three numeric matrices
#'
"Ancestry.list"

#' Phenotypes for parental lines of the American maize NAM dataset
#'
#' \code{Parents} is a data.frame with four rows. The first row corresponds to the
#' parental  inbred line common to all bi-parental populations of the American maize NAM experiment
#' (i.e., B73), whereas the three other rows correspond to alternative inbred parental lines used to
#' generate other bi-parental populations. For each parental inbred line, the phenotypic values
#' for five traits are provided, along with two extra columns:
#' Genotype (providing the line ID) and Weight (providing the weight associated
#' with the individual, see function EpiTest.fit for details).
#' @format A data.frame
#'
"Parents"
