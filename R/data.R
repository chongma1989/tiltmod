#' Marfan
#' 
#' Gene expression data includes the treatment group and 4132 gene expressions.
#' 
#' @format A data frame for 101 samples with 4133 variables: \code{treatment}, 
#'   \code{X1},...,\code{X4132}. \code{treatment} contains 41 samples from 
#'   the control group and 60 samples from the Marfan group.  
#'   
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2174953/}.   
"Marfan"

#' pickrell
#' 
#' The RNA-Seq profiles were made of cell lines derived from lymphoblastoid 
#' cells from 69 different Yoruba individuals from Ibadan, Nigeria. Pickrell
#' data consists of 40 females and 29 males for 17310 gene counts data, which
#' are well annotated and being at least 1 count-per-million (cpm) in at least
#' 20 samples. The raw RNA-Seq data for pickrell is available in R package 
#' \emph{tweeDEseqCountData}.
#' 
#' @format A DGEList S4 class, contains the gene count data, sample information, 
#'   and gene annotation data.  
#'   
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3089435/}.
"pickrell"