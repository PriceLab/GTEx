#----------------------------------------------------------------------------------------------------
#' @import methods
#' @import TrenaProject
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#'
#' @title GTEx-class
#'
#' @name GTEx-class
#' @rdname GTEx-class
#' @aliases GTEx
#' @exportClass GTEx
#'

.GTEx <- setClass("GTEx",
                  representation=representation(
                     state="environment",
                     matrix="matrix",
                     metadata="data.frame",
                     quiet="logical"
                     ))

#----------------------------------------------------------------------------------------------------
setGeneric('createSubMatrix',  signature='obj', function(obj, tissueName) standardGeneric('createSubMatrix'))
setGeneric('normalize',  signature='obj', function(obj, matrix, strategy) standardGeneric('normalize'))
#----------------------------------------------------------------------------------------------------
#' Define an object of class GTEx
#'
#' @description
#' Expression, variant and covariate data for the genes of interest (perhaps unbounded) for pre-term birth studies
#'
#' @rdname GTEx-class
#'
#' @export
#'
#' @return An object of the GTEx class
#'

GTEx <- function(quiet=TRUE)

{
   state <- new.env(parent=emptyenv())
   printf("loading full geneSymbol GTEx v8 matrix....")

   filename <- "~/github/GTEx/misc/raw/mtx.gtex.v8.geneSymbols.25337x17382.RData"
   mtx <- get(load(filename))

   filename <- "~/github/GTEx/misc/raw/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
   tbl.md <- read.table(filename, sep="\t", nrow=-1, header=TRUE, as.is=TRUE, quote="")
   rownames(tbl.md) <- gsub("-", ".", tbl.md$SAMPID, fixed=TRUE)

   .GTEx(state=state, matrix=mtx, metadata=tbl.md, quiet=quiet)

} # GTEx, the constructor
#----------------------------------------------------------------------------------------------------
#' create a tissue-specific expression matrix
#'
#' Given a string which matches all or part of a GTEx tissue name, return the matrix sith geneSymbol
#' names.
#'
#' @param obj An instance of the GTEx class
#' @param tissueName A character string, like "lung"
#'
#' @return An invisible expression matrix
#'
#' @export
#'
#' @aliases createSubMatrix
#' @rdname createSubMatrix

setMethod('createSubMatrix', 'GTEx',

    function(obj, tissueName){

       tbl.md <- obj@metadata
       mtx <- obj@matrix

       matching.tissue.names <- grep(tissueName, tbl.md$SMTS, ignore.case=TRUE, value=TRUE)
       if(!obj@quiet)
          print(table(matching.tissue.names))
       tissue.rows <- grep(tissueName, tbl.md$SMTS, ignore.case=TRUE)
       tissue.samples <- rownames(tbl.md)[tissue.rows]

       tissue.samples.with.expression <- intersect(tissue.samples, colnames(obj@matrix))
       length(intersect(tissue.samples, colnames(mtx)))  # [1] 578
       mtx.sub <- mtx[, tissue.samples.with.expression]
       dim(mtx.sub)  # [1] 56200   578
       invisible(mtx.sub)
       }) # createSubMatrix

#----------------------------------------------------------------------------------------------------
#' normalize
#'
#' transform counts using the specified strategy
#'
#' @param obj An instance of the GTEx class
#' @param matrix an expression matrix
#' @param strategy A character string: "asinh",
#'
#' @return An invisible expression matrix
#'
#' @export
#'
#' @aliases normalize
#' @rdname normalize

setMethod('normalize', 'GTEx',

    function(obj, matrix, strategy){

       stopifnot(strategy %in% c("asinh", "cory"))

       if(strategy == "asinh")
         mtx.out <- asinh(matrix)

       if(strategy == "cory"){
          minValue <- min(matrix[matrix > 0])
          if(minValue == 0)
             minValue <- .Machine$double.eps
          mtx.1 <- matrix + minValue
          mtx.2 <- log10(mtx.1)
          mtx.3 <- t(scale(t(mtx.2)))
          #means <- apply(mtx.2, 2, mean)
          mtx.3[is.na(mtx.3)] <- 0
          mtx.out <- mtx.3
          browser()
          xyz <- 99
          } # cory's strategy

       invisible(mtx.out)

       }) # normalize

#----------------------------------------------------------------------------------------------------
