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

se <- readRDS(file.path(system.file("extdata",
                                    package="IEOproject"),
                        "GSE79209.rds"))
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
mask <- rowMeans(assays(se)$logCPM) > 1
  se.filt <- se[mask, ]
  dge.filt <- dge[mask, ]
  dim(se.filt)
  
cpmcutoff <- round(10/min(dge$sample$lib.size/1e6), digits=1)
cpmcutoff

nsamplescutoff <- 5
nsamplescutoff
mask <- avgexp >= 4 & avgexp <= 10
se.filt <- se[mask, ]
dge.filt <- dge[mask, ]
dim(se.filt)

## ----lowexpgenes, echo=FALSE, out.width="600px", fig.cap="Distribution of lowly-expressed genes."----
par(mar=c(4, 5, 1, 1))
h <- hist(avgexp, xlab=expression("Expression level (" * log[2] * "CPM)"),
            main="", las=1, col="grey", cex.axis=1.2, cex.lab=1.5)
x <- cut(rowMeans(assays(se.filt)$logCPM), breaks=h$breaks)
lines(h$mids, table(x), type="h", lwd=10, lend=1, col="red")
legend("topright", c("All genes", "Filtered genes"), fill=c("grey", "red"))

## -----------------------------------------------------------------------------
dge.filt <- calcNormFactors(dge.filt)

## -----------------------------------------------------------------------------
assays(se.filt)$logCPM <- cpm(dge.filt, log=TRUE,
                              normalized.lib.sizes=TRUE)

## ----maplots, fig.height=18, fig.width=10, dpi=100, echo=FALSE, fig.cap="MA-plots of filtered and normalized expression values."----
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
table(se.filt$histology, se.filt$library_strategy)

## -----------------------------------------------------------------------------
table(se.filt$pack_years, se.filt$histology)

## ----sampleClustering, fig.height=8, fig.width=10, dpi=100, echo=FALSE, fig.cap="Figure S6: Hierarchical clustering of the samples. Labels correspond to pack-years and sample identifer, while colors indicate histology groups."----
logCPM <- cpm(dge.filt, log=TRUE, prior.count=3)
d <- as.dist(1-cor(logCPM, method="spearman"))
sampleClustering <- hclust(d)
batch <- as.integer(se.filt$histology)
sampleDendrogram <- as.dendrogram(sampleClustering, hang=0.1)
names(batch) <- colnames(se.filt)
outcome <- paste(se.filt$pack_years, colnames(se), sep="\n")
names(outcome) <- colnames(se.filt)
sampleDendrogram <- dendrapply(sampleDendrogram,
                               function(x, batch, labels) {
                                 if (is.leaf(x)) {
                                   attr(x, "nodePar") <- list(lab.col=as.vector(batch[attr(x, "label")]))
                                   attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                 }
                                 x
                               }, batch, outcome)
plot(sampleDendrogram, main="Hierarchical clustering of samples",
     cex=1)
legend("topright", levels(se.filt$histology),
       fill=seq_len(nlevels(se.filt$histology)))


## ----mdsPlot, fig.height=5, fig.width=8, dpi=100, echo=FALSE, fig.cap="Figure S7: Multidimensional scaling plot of the samples. Labels correspond to treatment and colors indicate sample group."----
outcome <- se.filt$pack_years
names(outcome) <- colnames(se.filt)
plotMDS(dge.filt, labels=outcome, col=batch)
legend("bottomleft", levels(se.filt$histology),
       fill=seq_len(nlevels(se.filt$histology)), inset=0.05)

## -----------------------------------------------------------------------------
copd <- rowMeans(logCPM[, se.filt$copd_status=="COPD"])
no_copd <- rowMeans(logCPM[, se.filt$copd_status=="No COPD"])

## -----------------------------------------------------------------------------
par(mfrow=c(1, 2))
plot(copd, no_copd, xlab="COPD", ylab="No COPD", pch=".", cex=4, las=1)
plot((copd+no_copd)/2, copd-no_copd, xlab="(COPD + No COPD)/2", ylab="COPD - No COPD", pch=".", cex=4, las=1)

## -----------------------------------------------------------------------------
se.filt$dysplasia[se.filt$histology == "Normal"] <- "No Dysplasia"
se.filt$dysplasia[se.filt$histology == "Hyperplasia"] <- "No Dysplasia"
se.filt$dysplasia[se.filt$histology == "Metaplasia"] <- "No Dysplasia"
se.filt$dysplasia[se.filt$histology == "Mild dysplasia"] <- "Dysplasia"
se.filt$dysplasia[se.filt$histology == "Moderate dysplasia"] <- "Dysplasia"
se.filt$dysplasia[se.filt$histology == "Severe dysplasia"] <- "Dysplasia"

## -----------------------------------------------------------------------------

#se.filt.dysplasia <- se.filt[, !is.na(se.filt$dysplasia)]
se.filt.dysplasia <- se.filt

no_dysplasia <- rowMeans(logCPM[, se.filt.dysplasia$dysplasia=="No Dysplasia"])
dysplasia <- rowMeans(logCPM[, se.filt.dysplasia$dysplasia=="Dysplasia"])


## -----------------------------------------------------------------------------
par(mfrow=c(1, 2))
plot(no_dysplasia, dysplasia, xlab="No Dysplasia", ylab="Dysplasia", pch=".", cex=4, las=1)
plot((no_dysplasia+dysplasia)/2, no_dysplasia-dysplasia, xlab="(No Dysplasia + Dysplasia)/2", ylab="No Dysplasia - Dysplasia", pch=".", cex=4, las=1)

## -----------------------------------------------------------------------------
mod <- model.matrix(~ se.filt$copd_status, colData(se.filt))
mod0 <- model.matrix(~ 1, colData(se.filt))

## ----message=FALSE------------------------------------------------------------
library(sva)
library(limma)

pv <- f.pvalue(assays(se.filt)$logCPM, mod, mod0)
sum(p.adjust(pv, method="fdr") < 0.05)

## ----pdistCOPD, echo=FALSE, out.width="600px", fig.cap="Distribution of raw p-values for an F-test on every gene between samples with COPD and samples without COPD."----
hist(pv, main="", las=1)

## ----message=FALSE------------------------------------------------------------
mask <- p.adjust(pv, method="fdr") < 0.5
DEgenesEGs <- names(pv)[mask]
DEgenesSyms <- mcols(se.filt)[DEgenesEGs, "symbol"]
DEgenesPvalue <- pv[mask]
DEgenesDesc <- mcols(se.filt)[DEgenesEGs, "description"]
DEgenesDesc <- sub(" \\[.+\\]", "", DEgenesDesc)
DEgenesTab <- data.frame(EntrezID=DEgenesEGs,
                         Symbol=DEgenesSyms,
                         Description=DEgenesDesc,
                         "P value"=DEgenesPvalue,
                         stringsAsFactors=FALSE, check.names=FALSE)
DEgenesTab <- DEgenesTab[order(DEgenesTab[["P value"]]), ] ## order by p-value
rownames(DEgenesTab) <- 1:nrow(DEgenesTab)

## ----echo=FALSE---------------------------------------------------------------
## generate full table in a CSV file and store it in the 'doc' directory
## twice, once in 'doc' to enable quickly look up during vignette editing
## and building with 'devtools::build_vignettes()' and a second time in
## 'inst/doc' to make these files available at install.
fnameCSV <- "DEgenesCOPD.csv"
fpathCSV <- file.path(path2pkg, "doc", fnameCSV)
write.csv(DEgenesTab, fpathCSV, row.names=FALSE)
fpathCSV <- file.path(path2pkg, "inst", "doc", fnameCSV)
write.csv(DEgenesTab, fpathCSV, row.names=FALSE)

## generate full table in HTML and store it into the 'doc' directory
## twice, just as we did with the CSV file. note that because the
## table caption is not translated from Markdown, but directly copied
## into HTML, we need to avoid using the '<' symbol, as in FDR < 10%,
## and put its HTML code instead (&lt;)
ktab <- kable(DEgenesTab, "html", escape=FALSE, row.names=TRUE,
              caption=sprintf("Differentially expressed genes. Differentially expressed genes between samples with COPD and samples without COPD, with FDR &lt; 5%% (CSV <a href=\"%s\" download>file</a>).",
                              fnameCSV))
ktab <- kable_styling(ktab,
                      bootstrap_options=c("stripped", "hover", "responsive"),
                      fixed_thead=TRUE)
fnameHTML <- "DEgenesCOPD.html"
fpathHTML <- file.path(path2pkg, "doc", fnameHTML)
save_kable(ktab, file=fpathHTML, self_contained=TRUE)
fpathHTML <- file.path(path2pkg, "inst", "doc", fnameHTML)
save_kable(ktab, file=fpathHTML, self_contained=TRUE)

## ----DEgenesCOPD, echo=FALSE--------------------------------------------------
ktab <- kable(DEgenesTab[1:10, ], "html", escape=FALSE, row.names=TRUE, 
              caption=sprintf("Differentially expressed genes. Top-10 differentially expressed genes with lowest p-value between samples with COPD and without COPD with FDR &lt; 5%%. To see the full list of DE genes, follow this <a href=\"%s\" target=\"_blank\">link</a> or download this CSV <a href=\"%s\" download>file</a>.",
                              fnameHTML, fnameCSV))
kable_styling(ktab, position="center")

## -----------------------------------------------------------------------------
mod <- model.matrix(~ se.filt.dysplasia$dysplasia, colData(se.filt.dysplasia))
mod0 <- model.matrix(~ 1, colData(se.filt.dysplasia))

## ----message=FALSE------------------------------------------------------------
library(sva)
library(limma)

pv <- f.pvalue(assays(se.filt.dysplasia)$logCPM, mod, mod0)
sum(p.adjust(pv, method="fdr") < 0.05)

## ----pdistDysplasia, echo=FALSE, out.width="600px", fig.cap="Distribution of raw p-values for an F-test on every gene between samples with Dysplasia and samples with Hyperplasia/Normal samples."----
sum(p.adjust(pv, method="fdr") < 0.05)
hist(pv, main="", las=1)

## ----message=FALSE------------------------------------------------------------
mask <- p.adjust(pv, method="fdr") < 0.05
DEgenesEGs <- names(pv)[mask]
DEgenesSyms <- mcols(se.filt.dysplasia)[DEgenesEGs, "symbol"]
DEgenesPvalue <- pv[mask]
DEgenesDesc <- mcols(se.filt.dysplasia)[DEgenesEGs, "description"]
DEgenesDesc <- sub(" \\[.+\\]", "", DEgenesDesc)
DEgenesTab <- data.frame(EntrezID=DEgenesEGs,
                         Symbol=DEgenesSyms,
                         Description=DEgenesDesc,
                         "P value"=DEgenesPvalue,
                         stringsAsFactors=FALSE, check.names=FALSE)

DEgenesTab <- DEgenesTab[order(DEgenesTab[["P value"]]), ] ## order by p-value
rownames(DEgenesTab) <- 1:nrow(DEgenesTab)

## ----echo=FALSE---------------------------------------------------------------
## generate full table in a CSV file and store it in the 'doc' directory
## twice, once in 'doc' to enable quickly look up during vignette editing
## and building with 'devtools::build_vignettes()' and a second time in
## 'inst/doc' to make these files available at install.
fnameCSV <- "DEgenesDysplasia.csv"
fpathCSV <- file.path(path2pkg, "doc", fnameCSV)
write.csv(DEgenesTab, fpathCSV, row.names=FALSE)
fpathCSV <- file.path(path2pkg, "inst", "doc", fnameCSV)
write.csv(DEgenesTab, fpathCSV, row.names=FALSE)

## generate full table in HTML and store it into the 'doc' directory
## twice, just as we did with the CSV file. note that because the
## table caption is not translated from Markdown, but directly copied
## into HTML, we need to avoid using the '<' symbol, as in FDR < 10%,
## and put its HTML code instead (&lt;)
ktab <- kable(DEgenesTab, "html", escape=FALSE, row.names=TRUE,
              caption=sprintf("Differentially expressed genes. Differentially expressed genes between samples with Dysplasia and samples with no Dysplasia, with FDR &lt; 5%% (CSV <a href=\"%s\" download>file</a>).",
                              fnameCSV))
ktab <- kable_styling(ktab,
                      bootstrap_options=c("stripped", "hover", "responsive"),
                      fixed_thead=TRUE)
fnameHTML <- "DEgenesDysplasia.html"
fpathHTML <- file.path(path2pkg, "doc", fnameHTML)
save_kable(ktab, file=fpathHTML, self_contained=TRUE)
fpathHTML <- file.path(path2pkg, "inst", "doc", fnameHTML)
save_kable(ktab, file=fpathHTML, self_contained=TRUE)

## ----DEgenesDysplasia, echo=FALSE---------------------------------------------
ktab <- kable(DEgenesTab[1:10, ], "html", escape=FALSE, row.names=TRUE, 
              caption=sprintf("Differentially expressed genes. Top-10 differentially expressed genes with lowest p-value between samples with Dysplasia and samples with No Dysplasia, with FDR &lt; 5%%. To see the full list of DE genes, follow this <a href=\"%s\" target=\"_blank\">link</a> or download this CSV <a href=\"%s\" download>file</a>.",
                              fnameHTML, fnameCSV))
kable_styling(ktab, position="center")

## -----------------------------------------------------------------------------
sessionInfo()

