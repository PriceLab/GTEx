library(GTEx)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("gtex"))
    gtex <- GTEx()
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_extractTissueSpecificExpressionMatrix()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   checkTrue("GTEx" %in% is(gtex))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_extractTissueSpecificExpressionMatrix <- function()
{
   printf("--- test_extractTissueSpecificExpressionMatrix")

   mtx.uterus <- createSubMatrix(gtex, "uterus")
   checkEquals(dim(mtx.uterus), c(25337, 142))

   mtx.lung   <- createSubMatrix(gtex, "lung")
   dim(mtx.lung)
   checkEquals(dim(mtx.lung), c(25337, 578))

   mtx.prostate   <- createSubMatrix(gtex, "prostate")
   checkEquals(dim(mtx.prostate), c(25337, 245))

} # test_extractSpecificExpressionMatrix
#------------------------------------------------------------------------------------------------------------------------
test_normalize <- function()
{
   printf("--- test_normalize")

   mtx.prostate   <- createSubMatrix(gtex, "prostate")
   mtx.asinh <- normalize(gtex, mtx.prostate, "asinh")
   checkTrue(all(mtx.asinh >= 0))
   checkEqualsNumeric(max(mtx.asinh), 12.18, tol=0.1)
   checkEquals(dim(mtx.prostate), dim(mtx.asinh))

   mtx.cory  <- normalize(gtex, mtx.prostate, "cory")
   checkTrue(all(mtx.cory >= -12))
   checkEqualsNumeric(max(mtx.cory), 15.6, tol=0.1)
   checkEquals(length(which(is.na(mtx.cory))), 0)
   checkEquals(dim(mtx.prostate), dim(mtx.cory))

} # test_normalize
#------------------------------------------------------------------------------------------------------------------------
# create and save matrices to be used by other packages and scripts
createMatrices <- function()
{
   mtx.counts   <- createSubMatrix(gtex, "lung")
   mtx.lung <- normalize(gtex, mtx.counts, "asinh")
   checkTrue(all(mtx.lung >= 0))
   checkEqualsNumeric(max(mtx.lung), 12.18, tol=0.1)
   save(mtx.lung, file="~/tmp/mtx.lung.RData")

} # createMatrices
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
