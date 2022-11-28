BiocManager::install("maftools")
library(maftools)
library(base)
setwd("...")

c<-read.table("TCGA.COAD.maf", sep="\t", quote = "",header=TRUE)#downloaded from dbGap as available in February 2018 (VarScan, public)
r<-read.table("TCGA.READ.maf", sep="\t", quote = "",header=TRUE)#downloaded from dbGap as available in February 2018 (VarScan, public)
colnames(c)
colnames(r)
m<-rbind(c,r)

g<-read.table("TCGA.COAD.metafile.tsv", sep="\t", header=TRUE)#this is the metafile with the TCGA clinical attributes and the T cell infiltration metrics

m$Tumor_Sample_Barcode<-strtrim(m$Tumor_Sample_Barcode, 16)

library(dplyr)
library(rms)

t<-subset(m, Tumor_Sample_Barcode %in% g$Tumor_Sample_Barcode)

describe(t$Variant_Classification)#useful description

h<-g$Tumor_Sample_Barcode[g$infiltration=="VeryHigh"]

l<-g$Tumor_Sample_Barcode[g$infiltration=="NoVeryHigh"]

cl<-subset(t, Tumor_Sample_Barcode %in% l)

ch<-subset(t, Tumor_Sample_Barcode %in% h)

library(rms)

Low<-read.maf(
  maf=cl,
  clinicalData = NULL,
  removeDuplicatedVariants = TRUE,
  useAll = TRUE,
  gisticAllLesionsFile = NULL,
  gisticAmpGenesFile = NULL,
  gisticDelGenesFile = NULL,
  gisticScoresFile = NULL,
  cnLevel = "all",
  cnTable = NULL,
  isTCGA = TRUE,
  vc_nonSyn=c(),
  verbose = TRUE
)

High<-read.maf(
  maf=ch,
  clinicalData = NULL,
  removeDuplicatedVariants = TRUE,
  useAll = TRUE,
  gisticAllLesionsFile = NULL,
  gisticAmpGenesFile = NULL,
  gisticDelGenesFile = NULL,
  gisticScoresFile = NULL,
  cnLevel = "all",
  cnTable = NULL,
  isTCGA = TRUE,
  vc_nonSyn=c(),
  verbose = TRUE
)

coOncoplot(m1 = Low, m2 = High, m1Name = 'Not very high', m2Name = 'Very high') 
quartz.save("Oncoplot_comp.PDF", type="PDF")
dev.off()

oncoplot(Low)
quartz.save("Oncoplot_Low.PDF", type="PDF")
dev.off()

oncoplot(High)
quartz.save("Oncoplot_High.PDF", type="PDF")
dev.off()

BiocManager::install("BSgenome")
BiocManager::install("NMF")
library("BiocManager")
library("BSgenome")
library("NMF")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BSgenome::available.genomes()
Low #to see what is reference genome, in this case it is GRCh38
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
library("BSgenome.Hsapiens.NCBI.GRCh38")
coadL.tnm <- trinucleotideMatrix(maf = Low, ref_genome = 'BSgenome.Hsapiens.NCBI.GRCh38',
                                prefix = 'chr', add = FALSE, useSyn = TRUE)

coadH.tnm <- trinucleotideMatrix(maf = High, ref_genome = 'BSgenome.Hsapiens.NCBI.GRCh38',
                                prefix = 'chr', add = FALSE, useSyn = TRUE)
coadL.sign <- extractSignatures(mat = coadL.tnm, plotBestFitRes = FALSE, n = 3, pConstant = 0.01)
coadH.sign <- extractSignatures(mat = coadH.tnm, plotBestFitRes = FALSE, n = 3, pConstant = 0.01)
coad.v3.cosmL = compareSignatures(nmfRes = coadL.sign, sig_db = "SBS")
coad.v3.cosmH = compareSignatures(nmfRes = coadH.sign, sig_db = "SBS")
maftools::plotSignatures(nmfRes = coadL.sign, title_size = 1.2, sig_db = "SBS")
quartz.save("SBS_signature_COREAD_NotVeryHighInf.PDF", type="PDF")
maftools::plotSignatures(nmfRes = coadH.sign, title_size = 1.2, sig_db = "SBS")
quartz.save("SBS_signature_COREAD_VeryHighInf.PDF", type="PDF")
pt.vs.rt <- mafCompare(m1 = Low, m2 = High, m1Name = 'Not very high', m2Name = 'Very high', minMut = 5)
print(pt.vs.rt)
class(pt.vs.rt)
head(pt.vs.rt$results$Hugo_Symbol)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.005, geneFontSize = 0.3,
           titleSize = 0.8,
           lineWidth = 1)
quartz.save("ForrestPlot_mutations_InfiltrationGroups.PDF", type="PDF")
dev.off()

subsetMaf(maf=High, genes=c("BRAF"), mafObj = FALSE)
braf<-subsetMaf(maf=High, genes=c("BRAF"), mafObj = FALSE)

braf$Tumor_Sample_Barcode
list<-strtrim(braf$Tumor_Sample_Barcode, 16)

list %in% m$Tumor_Sample_Barcode
colnames(m)
(braf$Tumor_Sample_Barcode %in% strtrim(m$Tumor_Sample_Barcode, 12))


