#' @import utils
## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))



#' Chi-Squared Mixtures Distribution Function
#'
#' @param quant a quantile
#' @param s number of fixed effects to be tested
#' @param q  number of random effects to be tested
#' @param lower.tail logical. if TRUE (default), probabilities are \eqn{P[X<=x]}, otherwise \eqn{P[X>x]}.
#' @description The approximate null distribution of a likelihood ratio for 2 nested mixed models,
#' where both fixed and random effects are tested simultaneously, is a very specific mixture
#' of chi-square distributions. It depends on both the number of random effects
#' and the number of fixed effects to be tested simultaneously. Note that this function is a copy
#' of the \code{pchisqmix} function of the \code{TcGSA} package.
#' @return a probability.
#' @importFrom stats pchisq
pchisqmix <- function (quant, s, q, lower.tail = TRUE){
  if (q > 0) {
    mixprobs <- numeric(q + 1)
    for (k in s:(q + s)) {
      mixprobs[k - s + 1] <- choose(q, k - s) * 2^(-q)
    }
    mix <- mixprobs
  }
  else {
    mix = 1
  }
  res <- numeric(length(mix))
  for (k in (s:(q + s))) {
    res[k - s + 1] <- mix[k - s + 1] * stats::pchisq(quant,
                                                     df = k, lower.tail = lower.tail)
  }
  return(sum(res)/sum(mix))
}


#' Compute covariance matrices
#'
#' @description This function computes the covariance matrices associated with
#' the variance components of the EpiTest model: the segregation variance,
#' the (segregation \code{x} segregation) variance and the error variance

#' @param Ancestry a matrix with 0/1 entries representing the allele ancestries
#' of a bi-parental population at a set of markers  (0: homozygous for parent A
#' alleles and 1: homozygous for parent B alleles)
#'
#' @return A list of three covariance matrices
#'
#' @export
#'
#' @examples
#' ## One bi-parental population, no weighting
#' data(Ancestry.list)
#' Ancestry <- Ancestry.list[[1]]
#' VarList <- VarList.comp(Ancestry = Ancestry)
#' purrr::map(VarList,~.x[1:5,1:5])
VarList.comp <- function(Ancestry){
  pi_i <- rowMeans(Ancestry)
  Delta <- tcrossprod(Ancestry-pi_i)/ncol(Ancestry)
  Delta2 <- Delta^2
  res <- list(sigma2_S = Delta,
              sigma2_SxS = Delta2,
              sigma2_E = diag(nrow(Delta)))
  return(res)
}


#' Inference and Test
#'
#' @param Ancestry a numeric matrix with 0/1 entries representing the allele ancestries
#' of a bi-parental population at a set of markers (0: homozygous for parent A alleles and
#' 1: homozygous for parent B alleles)
#' @param Pheno a data.frame with a \code{Genotype} character column indicating the names
#' of individuals, a numeric column for each trait to be analyzed, and possibly
#' a \code{Weight} numeric column indicating the number of observation for each
#' individual. Missing values are discarded if present for a trait
#' @param ParentName a character vector of length 2 providing the parent names where the
#' first element should be named P0 (parent whose allele ancestries are coded 0) and
#' the second element should be named P1 (parent whose allele ancestries are coded 1)
#' @param Parents (optional) a data.frame with a \code{Genotype} character column indicating
#' the names of parents, a \code{Family} character column indicating the names of
#' families/populations, a numeric column for each trait to be analyzed, and possibly
#' a \code{Weight} numeric column indicating the number of observation for each parent
#' (see Details)
#' @param Trait a string indicating the name of the phenotypic trait
#' @param Weight a boolean indicating whether weights indicating the number of
#' observations should be used for the inference. If TRUE, a column \code{Weight}
#' should be found in Pheno and Parents (see Details)
#'
#' @return A list of four items: Beta, Sigma2, Test_fixed, Test_random
#'
#' @description This function fits the EpiTest model on a bi-parental population
#' dataset, then outputs the estimates of the fixed and random effects, along
#' with their associated test results.
#' @details Depending on cases, the data may be organized in two different ways.
#'
#' In Case 1, the \code{Ancestry} matrix and the \code{Pheno} data.frame respectively provide
#' the whole ancestry and phenotypic data, including the parental data.
#'
#' In Case 2, the phenotypic data are split into two parts, the parental data being included
#' in a separate \code{Parents} data.frame.
#'
#' In both cases, one may provide a set of weights associated to each individual by including
#' a \code{Weight} column in the \code{Pheno} and \code{Parents} data.frames. These weights
#' correspond to the number of observations that where used to compute the parental BLUPS/mean values.
#' These weights must be provided as a \code{Weight} column in the \code{Parents} and \code{Pheno}
#' data.frame.
#'
#' The fitted EpiTest model includes three variance components: the segregation variance,
#' the (segregation \code{x} segregation) variance and the error variance.
#' For each genetic variance a likelihood ratio test is performed.
#'
#' The fitted EpiTest model includes three fixed parameters: the intercept,
#' the linear regression coefficient (beta) on parent proportions and the quadratic
#' regression coefficient (delta) on parent proportions  that only involves epistatic effects.
#' Each fixed parameter is tested using a Wald test.
#'
#' Additionally, epistasis can be tested by testing that both the (segregation \code{x} segregation) variance component
#' and the quadratic mean component (delta) are null through a likelihood ratio test. Note that in this
#' case the LRT is based on the full (i.e. unrestricted) likelihoods.
#'
#' The function outputs a list of five items: Beta is the vector of estimated
#' fixed effects, Sigma2 is the vector of estimated variances, Test_full, Test_fixed,
#' Text_random are 3 data.table that provide the results of the tests for no epistasis (jointly tested on fixed and variance
#' components), and for the fixed and
#' random effects, respectively.
#' @import MM4LMM
#' @import purrr
#' @import dplyr
#' @export
#'
#' @examples
#' ## One bi-parental population, no weighting and no parental phenotypes
#' data(Pheno.list)
#' data(Ancestry.list)
#' Ancestry <- Ancestry.list[[1]]
#' Pheno <- Pheno.list[[1]]
#' ParentName <- c(P0 = 'B73',P1 = 'B97')
#' ETest.1 <- EpiTest.fit(Ancestry = Ancestry,
#'                        Pheno = Pheno,
#'                        ParentName = ParentName,
#'                        Trait = "GDD_DTA")
#'
#' ## One bi-parental population, with weights and parental phenotypes
#' data(Parents)
#' ETest.2 <- EpiTest.fit(Ancestry = Ancestry,
#'                        Pheno = Pheno,
#'                        ParentName = ParentName,
#'                        Trait = "GDD_DTA",
#'                        Parents = Parents,
#'                        Weight=TRUE)
#'
#' ## Full NAM analysis, with weights and parental phenotypes
#' Parent.list <- Parents$Genotype[-1]
#' names(Parent.list) <- Parents$Family[-1]
#' ETest.nam <- purrr::imap(Parent.list, ~ EpiTest.fit(Ancestry = Ancestry.list[[.y]],
#'                                                     Pheno = Pheno.list[[.y]],
#'                                                     ParentName=c(P0 = 'B73',P1 = .x),
#'                                                     Parents = Parents,
#'                                                     Trait = 'GDD_DTA',
#'                                                     Weight = TRUE))
EpiTest.fit <- function(Ancestry,Pheno,Trait,ParentName,Parents=NULL,Weight=FALSE){

  ## Checks
  #Ancestry
  if(!is.matrix(Ancestry)|!is.numeric(Ancestry)){
    stop("Ancestry should be a numeric matrix")
  }
  if(any(is.na(Ancestry))){
    stop("Missing values are not allowed in the Ancestry matrix, consider imputation")
  }
  if(max(Ancestry,na.rm = T)>1){
    stop(paste("Ancestry matrix should not have values over 1, the maximum value found here is:",max(Ancestry,na.rm = T)))
  }
  if(min(Ancestry,na.rm = T)<0){
    stop(paste("Ancestry matrix should not have values below 0, the minimum value found here is:",min(Ancestry,na.rm = T)))
  }
  if(length(which(!(Ancestry %in% c(0,1))))>0){
    message("Some Ancestry entries are different from 0 or 1...")
  }
  #Trait
  if(!is.character(Trait)){
    stop("Trait should be a string")
  }
  #Pheno
  if(!is.data.frame(Pheno)){
    stop("Pheno should be a data.frame")
  }
  if(!("Genotype"%in%colnames(Pheno))){
    stop("Pheno should include a Genotype column")
  }
  if(!(Trait%in%colnames(Pheno))){
    stop(paste("Column",Trait,"not found in data.frame Pheno"))
  }
  if(any(is.na(Pheno))){
    message(paste("Missing values were detected in Pheno for",Trait,"and will be discarded from the analysis"))
  }
  #ParentName
  if((!is.character(ParentName))|(length(ParentName)!=2)|(!length(names(ParentName))==2)|!all(names(ParentName)%in%c("P0","P1"))){
    stop("ParentName should be a two components character vector with names P0 and P1")
  }
  if(!(ParentName[["P0"]]%in%c(Parents$Genotype,Pheno$Genotype))){
    message(paste("No phenotypic value found for parent",ParentName[1]))
  }
  if(!(ParentName[["P1"]]%in%c(Parents$Genotype,Pheno$Genotype))){
    message(paste("No phenotypic value found for parent",ParentName[2]))
  }
  #Parents
  if(!is.null(Parents)){
    if(!is.data.frame(Parents)){
      stop("Parents should be a data.frame")
    }
    if(any(is.na(Parents[[Trait]]))){
      message(paste("Missing values were detected in Parents for",Trait,"and will be discarded from the analysis"))
    }
  }
  #Weight
  if(!is.logical(Weight)){
    stop("Weight should be logical")
  }
  if(Weight&!is.null(Parents)&(!("Weight"%in%colnames(Parents))|!("Weight"%in%colnames(Pheno)))){
    stop("A Weight column should be included in Pheno and Parents if Weight is TRUE")
  }
  if(Weight&is.null(Parents)&(!("Weight"%in%colnames(Pheno)))){
    stop("A Weight column should be included in Pheno if Weight is TRUE")
  }
  if(Weight&!is.null(Parents)&(!is.numeric(Parents$Weight)|!is.numeric(Pheno$Weight))){
    stop("The Weight column in Pheno and Parents should be numeric")
  }
  if(Weight&is.null(Parents)&(!is.numeric(Pheno$Weight))){
    stop("The Weight column in Pheno should be numeric")
  }

  ## Order Ancestry and Pheno files
  Pheno <- Pheno %>%
    mutate_if(is.factor,as.character)
  Pheno <- Pheno[!is.na(Pheno[[Trait]]),]
  IndInter <- intersect(rownames(Ancestry),Pheno$Genotype)
  if(length(IndInter)<2){
    stop("Check the correspondance between the row names of Ancestry and the Genotype column of Pheno")
  }
  Pheno <- Pheno[match(IndInter,Pheno$Genotype),]
  Ancestry <- Ancestry[match(IndInter,rownames(Ancestry)),]
  if(!is.null(Parents)){
    Parents <- Parents %>%
      mutate_if(is.factor,as.character)
    Parents <- Parents[!is.na(Parents[[Trait]][match(ParentName,Parents$Genotype)]),]
    if(nrow(Parents)==0){
      Parents <- NULL
      message(paste("Missing values for both parents, Parents data.frame ignored for ",Trait))
    }
  }

  ## Build ParentProp and VarList from allele ancestries
  if(is.null(Parents)){
    Ancestry_all <- Ancestry
    if(Weight){
      Weights <- Pheno$Weight
    }
  }else{
    if(all(ParentName%in%Parents$Genotype)){
      Ancestry_all <- rbind(P0 = 0,P1 = 1,Ancestry)
      if(Weight){
        Weights <- c(Parents$Weight[match(ParentName[c("P0","P1")],Parents$Genotype)],Pheno$Weight)
      }
    }else if(ParentName[["P0"]]%in%Parents$Genotype){
      Ancestry_all <- rbind(P0 = 0,Ancestry)
      if(Weight){
        Weights <- c(Parents$Weight[match(ParentName[c("P0")],Parents$Genotype)],Pheno$Weight)
      }
    }else if(ParentName[["P1"]]%in%Parents$Genotype){
      Ancestry_all <- rbind(P1 = 1,Ancestry)
      if(Weight){
        Weights <- c(Parents$Weight[match(ParentName[c("P1")],Parents$Genotype)],Pheno$Weight)
      }
    }
  }
  ParentProp <- rowMeans(Ancestry_all) %>%
    data.frame(mu = 1,beta = .,delta = .^2) %>%
    as.matrix()
  VarList <- VarList.comp(Ancestry_all)
  if(Weight){
    VarList$sigma2_E <- diag(1/Weights)
  }
  VarList <- map(VarList, ~ .x/mean(diag(.x)))

  ## Get the phenotypes
  if(is.null(Parents)){
    Y <- Pheno[[Trait]]
  } else {
    Y <- c(Parents[[Trait]][Parents$Genotype==ParentName[["P0"]]],
           Parents[[Trait]][Parents$Genotype==ParentName[["P1"]]],
           Pheno[[Trait]])
  }

  ## Fit the 3 models and test variance components
  VarListComp <- list(M1 = 3,M2 = c(1,3),M3 = 1:3)
  Models <- purrr::map(VarListComp, ~ MM4LMM::MMEst(Y=Y,
                                                    Cofactor = ParentProp,
                                                    VarList = VarList[.x],
                                                    Verbose = FALSE))
  Test_random <- map_dbl(Models, ~ .$NullModel$`LogLik (Reml)`) %>%
    {c(sigma2_S=-2*(.[1]-.[2]), sigma2_SxS=-2*(.[2]-.[3]))} %>%
    replace(.,which(.<0),0) %>%
    data.frame(LR_stat=.,pval=0.5*stats::pchisq(.,1,lower.tail = F) + 0.5*(.==0),df=1)
  rownames(Test_random) <- names(VarList)[1:2]

  ## Test fixed effects with all variance components
  Test_fixed <- MM4LMM::AnovaTest(Models$M3,Type = "TypeI")[[1]] %>% as.data.frame
  Beta <- Models$M3$NullModel$Beta
  Sigma2 <- Models$M3$NullModel$Sigma

  ## Test absence of epistasis in both mean and variance
  NullModel.Reml <- MM4LMM::MMEst(Y=Y,
                                  Cofactor = ParentProp[,1:2],
                                  VarList = VarList[c(1,3)],
                                  MaxIter = 200,
                                  Verbose = FALSE)
  NullModel.ML <- MM4LMM::MMEst(Y=Y,
                                Cofactor = ParentProp[,1:2],
                                VarList = VarList[c(1,3)],
                                Init = NullModel.Reml$NullModel$Sigma2,
                                Method = 'ML',MaxIter = 1,
                                Verbose = FALSE)
  AltModel.Reml <- MM4LMM::MMEst(Y=Y,
                                 Cofactor = ParentProp,
                                 VarList = VarList[1:3],
                                 MaxIter = 200,
                                 Verbose = FALSE)
  AltModel.ML <- MM4LMM::MMEst(Y=Y,
                               Cofactor = ParentProp,
                               VarList = VarList[1:3],
                               Init = AltModel.Reml$NullModel$Sigma2,
                               Method = 'ML', MaxIter = 1,
                               Verbose = FALSE)

  LRT.ep <- -2*(NullModel.ML$NullModel$`LogLik (ML)` - AltModel.ML$NullModel$`LogLik (ML)`)
  Pval.ep <- pchisqmix(quant = LRT.ep,s = 1,q = 1,lower.tail = FALSE)
  Test.ep <- data.frame(LR_stat=LRT.ep,pval=Pval.ep)

  ## Return output
  res <- list(Beta = Beta,
              Sigma2 = Sigma2,
              Test_full = Test.ep,
              Test_fixed = Test_fixed,
              Test_random = Test_random)
  return(res)
}

