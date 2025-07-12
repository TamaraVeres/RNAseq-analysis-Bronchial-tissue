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


## ----message=FALSE------------------------------------------------------------
library(SummarizedExperiment)

se <- readRDS(here("GSE79209.rds"))
se


## -----------------------------------------------------------------------------
head(rowData(se))

## -----------------------------------------------------------------------------
dim(colData(se))
head(colData(se), n=3)

## -----------------------------------------------------------------------------
colnames(colData(se))

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
interesting_se <- data.frame(row.names=se$title)
technical_se <- data.frame(row.names=se$title)
for(column in colnames(colData(se))){
  se[[column]] <- factor(se[[column]])
  if(length(levels(se[[column]]))>1 && length(levels(se[[column]]))<81){
    interesting_se[[column]] <- se[[column]]
  } else{
    technical_se[[column]] <- se[[column]]
  }
}

## -----------------------------------------------------------------------------
head(interesting_se)

## -----------------------------------------------------------------------------
se$sex <- gsub("Sex: ", "", se$characteristics_ch1.2)
se$sex <- factor(se$sex, levels = c("Male", "Female"))

se$smoker <- gsub("smoking status: ", "", se$characteristics_ch1.3)
se$smoker <- factor(se$smoker, levels = c("Current smoker", "Former smoker"))

se$copd_status <- gsub("copd status: ", "", se$characteristics_ch1.5)
se$copd_status <- factor(se$copd_status, levels = c("COPD", "No COPD"))

se$histology <- gsub("max histology: ", "", se$characteristics_ch1.6)
se$histology <- factor(se$histology, levels(se$histology) <- c("Hyperplasia", "Metaplasia", "Mild dysplasia", "Moderate dysplasia", "Normal", "Severe dysplasia"))

se$dysplasia <- gsub("Subject L.*, ", "", se$title)
se$dysplasia <- factor(se$dysplasia)

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


## -----------------------------------------------------------------------------

se$age <- gsub("age: ", "", se$characteristics_ch1.1)
class(se$age)
se$age <- as.numeric(se$age)

quartiles <- quantile(se$age, probs = c(0.25, 0.5, 0.75))
quartiles


## -----------------------------------------------------------------------------
se$age_intervals <- cut(se$age,
                          breaks = c(-Inf, quartiles[1], quartiles[2], quartiles[3], Inf),
                          labels = c("age <= 58", "58 < age <= 64", "64 < age <= 69", "69 < age"),
                          include.lowest = TRUE)


## ----pheno, echo=FALSE, message=FALSE-----------------------------------------
tmpdf <- data.frame(Identifer=colnames(se),
                    Age=se$age_intervals,
                    Sex=se$sex,
                    "Smokeing status"=se$smoker,
                    "Pack Years"=se$pack_years,
                    "COPD Status"=se$copd_status,
                    Histology=se$histology,
                    Dysplasia=se$dysplasia,
                    check.names=FALSE)
ktab <- kable(tmpdf, caption="Phenotypical variables.")
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
mask_initial <- rowMeans(assays(se)$logCPM) > 1
se_initial_filt <- se[mask_initial, ]
dge_initial_filt <- dge[mask_initial, ]
dim(se_initial_filt)

cpmcutoff <- round(10/min(dge$samples$lib.size/1e6), digits=1)
cpmcutoff

if(any(is.na(se$histology))) {
  stop("Missing values found in the histology column.")
}

pml_status <- se$histology %in% c("Mild dysplasia", "Moderate dysplasia", "Severe dysplasia")
nsamplescutoff <- sum(pml_status)
nsamplescutoff

pml_indices <- which(pml_status)

cpm_pml <- cpm(dge_initial_filt)[, pml_indices]

mask_final <- rowSums(cpm_pml > cpmcutoff) >= nsamplescutoff
sum(mask_final)
se.filt <- se_initial_filt[mask_final, ]
dge.filt <- dge_initial_filt[mask_final, ]
dim(se.filt)

## ----lowexpgenes, echo=FALSE, out.width="600px", fig.cap="Distribution of lowly-expressed genes."----
par(mar=c(4, 5, 1, 1))
avgexp <- rowMeans(assays(se)$logCPM)
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

