library(GTEx)
library(rnaSeqNormalizer)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("gtex"))
    gtex <- GTEx()
#------------------------------------------------------------------------------------------------------------------------
metadata.filename <- "~/github/GTEx/misc/raw/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
stopifnot(file.exists(metadata.filename))
tbl.md <- read.table(metadata.filename, sep="\t", nrow=-1, header=TRUE, as.is=TRUE, quote="")
dim(tbl.md)


   #----------------------------------------------------------
   # used to use the SMTSD column, but no longer (8 nov 2019)
   #----------------------------------------------------------

# tbl.xtab <- as.data.frame(table(tbl.md$SMTSD))
# colnames(tbl.xtab) <- c("tissue", "count")
# tbl.xtab <- tbl.xtab[order(tbl.xtab$count, decreasing=TRUE),]
# dim(tbl.xtab)
# dim(subset(tbl.xtab, count >= 100))   # 50

   #-------------------------------------------
   # simpler names come from the SMTS column
   #-------------------------------------------

tbl.xtab <- as.data.frame(table(tbl.md$SMTS))
colnames(tbl.xtab) <- c("tissue", "count")
tbl.xtab <- tbl.xtab[order(tbl.xtab$count, decreasing=TRUE),]
dim(tbl.xtab)
subset(tbl.xtab, count <= 100)   # loses only 3: bladder (21), cervix uteri (19), fallopian tube (9)

tissues <- as.character(subset(tbl.xtab, count >= 100)$tissue)

for(tissue in tissues){
   mtx <- createSubMatrix(gtex, tissue)
   printf("%20s: %d %d", tissue, nrow(mtx), ncol(mtx))
   if(ncol(mtx) >= 100){
      normalizer <- rnaSeqNormalizer.gtex(mtx, algorithm="vst", duplicate.selection.statistic="median")
      mtx.vst.median <- getNormalizedMatrix(normalizer)
      tissue.name <- gsub(" ", ".", tissue)
      tissue.name <- tolower(tissue.name)
      save(mtx.vst.median, file=sprintf("%s-vst-median.RData", tissue.name))
      } # at least 100 samples
   } # for tissue





