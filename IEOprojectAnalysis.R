## ----setup, echo=FALSE, cache=FALSE-------------------------------------------
library(knitr) ## kable()
library(kableExtra) ## kable_styling(), save_kable()
library(here) ## here()
library(usethis) ## use_directory()

knitr::opts_chunk$set(
  collapse=TRUE,
  comment="",
  fig.align="center",
  fig.wide=TRUE,
  cache=FALSE
)

## this option avoid use_directory() being verbose
options(usethis.quiet=TRUE)

## create these paths at build time if they do not exist
use_directory(file.path("doc"))
use_directory(file.path("inst", "doc"))

## fetch the package root directory
path2pkg <- here()


## ----eval=FALSE---------------------------------------------------------------
#  devtools::build_vignettes()

## ----eval=FALSE---------------------------------------------------------------
#  devtools::document()

## ----message=FALSE------------------------------------------------------------
library(SummarizedExperiment)
library(here)

se <- readRDS(here("GSE79209.rds"))
se


## -----------------------------------------------------------------------------
head(rowData(se))

## -----------------------------------------------------------------------------
dim(colData(se))
head(colData(se), n=3)

## -----------------------------------------------------------------------------
length(unique(se$geo_accession))
table(lengths(split(colnames(se), se$geo_accession)))

## ----message=FALSE------------------------------------------------------------
library(edgeR)

dge <- DGEList(counts=assays(se)$counts, genes=rowData(se))
head(dge$samples, 3)
dim(dge)

## -----------------------------------------------------------------------------
assays(se)$logCPM <- cpm(dge, log=TRUE)
assays(se)$logCPM[1:5, 1:5]

## -----------------------------------------------------------------------------
# Delete this at the end

table(se$characteristics_ch1.1) # age or table(se$age.ch1)
table(se$characteristics_ch1.2) # sex or table(se$Sex.ch1)
table(se$characteristics_ch1.3) # smoking status (current - former) or table(se$smoking.status.ch1)
table(se$characteristics_ch1.4) # pack-years (how many cigarettes you have smoked in your lifetime 1 pack = 20 ciggs) or table(se$pack.years.ch1)
table(se$characteristics_ch1.5) # copd status or table(se$copd.status.ch1)
table(se$characteristics_ch1.6) # max histology - hyperplasia, metaplasia, mild dysplasia, moderate dysplasia, normal and severe dysplasia or table(se$max.histology.ch1)
table(se$characteristics_ch1.7) # dysplasia status - dysplasia, NA, normal or table(se$dysplasia.status.ch1)
table(se$age.ch1) # mean gc content or table(se$mean.gc.content.ch1)

# Delete this at the end

## -----------------------------------------------------------------------------
table(se$characteristics_ch1.4)
table(se$characteristics_ch1.5)

## -----------------------------------------------------------------------------
table(se$characteristics_ch1.6)

## -----------------------------------------------------------------------------
se$copd_status <- gsub("copd status: ", "", se$characteristics_ch1.5)
se$copd_status <- factor(se$copd_status, levels = c("COPD", "No COPD"))

se$histology <- gsub("max histology: ", "", se$characteristics_ch1.6)
se$histology <- factor(se$histology, levels(se$histology) <- c("Hyperplasia", "Metaplasia", "Mild dysplasia", "Moderate dysplasia", "Normal", "Severe dysplasia"))

## -----------------------------------------------------------------------------
se$characteristics_ch1.4 <- gsub("pack years: ", "", se$characteristics_ch1.4)
class(se$characteristics_ch1.4)
se$characteristics_ch1.4 <- as.numeric(se$characteristics_ch1.4)

quartiles <- quantile(se$characteristics_ch1.4, probs = c(0.25, 0.5, 0.75))


## -----------------------------------------------------------------------------

se$pack_years <- cut(se$characteristics_ch1.4,
                          breaks = c(-Inf, quartiles[1], quartiles[2], quartiles[3], Inf),
                          labels = c("Low", "Medium", "High", "Very High"),
                          include.lowest = TRUE)


## ----pheno, echo=FALSE, message=FALSE-----------------------------------------
tmpdf <- data.frame("Identifer"=colnames(se),
                    "Pack Years"=se$pack_years,
                    Histology=se$histology,
                    "COPD Status"=se$copd_status,
                    check.names=FALSE)
ktab <- kable(tmpdf, caption="Phenotypic variables.")
kable_styling(ktab, position="center")

## ----libsizes, echo=FALSE, height=8, width=8, out.width="600px", fig.cap="Library sizes in increasing order."----
par(mar=c(7, 5, 2, 2))
ord <- order(dge$sample$lib.size/1e6)
ordmreads <- dge$sample$lib.size[ord]/1e6
names(ordmreads) <- colnames(se)[ord]
bp <- barplot(ordmreads, las=1, ylab="Millions of reads",
              xlab="", col=c("blue", "green", "red", "purple", "orange", "pink")[factor(se$histology[ord])], las=2, ylim = c(0, 50))
legend("topleft", c("Hyperplasia", "Metaplasia", "Mild dysplasia", "Moderate dysplasia", "Normal", "Severe dysplasia"), fill=c("blue", "green", "red", "purple", "orange", "pink"), inset=0.01, cex=0.85)


## ----distRawExp, echo=FALSE, fig.height=5, fig.width=5, out.width="600px", fig.cap="Non-parametric density distribution of expression profiles per sample.", message=FALSE----
library(geneplotter)
par(mar=c(4, 5, 1, 1))
lst <- as.list(as.data.frame(assays(se)$logCPM))
multidensity(lst, xlab="log 2 CPM", legend=NULL,
             main="", las=1)

## ----exprdist, echo=FALSE, out.width="600px", fig.cap="Distribution of average expression level per gene."----
avgexp <- rowMeans(assays(se)$logCPM)
hist(avgexp, xlab="log2 CPM", main="", las=1)

## -----------------------------------------------------------------------------
mask <- filterByExpr(dge, group=se$samplegroup)
se.filt <- se[mask, ]
dim(se.filt)
dge.filt <- dge[mask, ]
dim(dge.filt)

## -----------------------------------------------------------------------------
dge.filt <- calcNormFactors(dge.filt)


## -----------------------------------------------------------------------------
assays(se.filt)$logCPM <- cpm(dge.filt, log=TRUE,
                              normalized.lib.sizes=TRUE)

## ----maPlots, fig.height=18, fig.width=10, dpi=100, echo=FALSE, fig.cap="MA-plots of filtered and normalized expression values."----
par(mfrow=c(5, 4), mar=c(4, 5, 3, 1))
for (i in 1:ncol(se.filt)) {
  A <- rowMeans(assays(se.filt)$logCPM)
  M <- assays(se.filt)$logCPM[, i] - A
  smoothScatter(A, M, main=colnames(se.filt)[i], las=1)
  abline(h=0, col="blue", lwd=2)
  lo <- lowess(M ~ A)
  lines(lo$x, lo$y, col="red", lwd=2)
}

## -----------------------------------------------------------------------------
table(se.filt$cell_line, se.filt$treatment)

