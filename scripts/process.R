.libPaths("/pasteur/zeus/projets/p01/uFlu/reassortment_project/resources/R_libs/4.1.0")

##########################################################################################################
##### This script has been written by Mary O'Neill in August 2021 for the uFlu project.
##### It should be executed by a shell script for each sample.
##########################################################################################################

library(SingleCellExperiment)
library(scran)
library(scater)
library(dplyr)
library(tidyr)
library(data.table)
library(DropletUtils)
library(SoupX)
library(cowplot)
library(uwot)

PROJECT="NGS11"
ALIGNDIR="/pasteur/zeus/projets/p01/uFlu/reassortment_project/02_STARsolo_outputs"
TASKDIR="/pasteur/zeus/projets/p01/uFlu/reassortment_project/03_custom_analyses"
REFDIR="/pasteur/zeus/projets/p01/uFlu/reassortment_project/resources/Ref_data"

args <- commandArgs(TRUE)
SAMPLE=args[1]

set.seed(1000)

#make theme for figures
theme_mary <- function() {
  font <- "Arial"
  
  theme_light() + 
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=14),
          #axis.text.x=element_text(angle=45, hjust=1),
          legend.position="bottom",
          legend.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA))
}

#load the data selected by the empty drops method
eD50 <- read10xCounts(sprintf('%s/%s/%s/%s.Solo.out/Gene/emptyDrops_50counts', ALIGNDIR, PROJECT, SAMPLE, SAMPLE))
eD75 <- read10xCounts(sprintf('%s/%s/%s/%s.Solo.out/Gene/emptyDrops_75counts', ALIGNDIR, PROJECT, SAMPLE, SAMPLE))

#read in whole dataset
sce=read10xCounts(sprintf('%s/%s/%s/%s.Solo.out/Gene/raw', ALIGNDIR, PROJECT, SAMPLE, SAMPLE)) #full dataset, very big

#filter out barcodes with no counts
cnts = colSums(counts(sce))
sce = sce[,cnts>0] 
saveRDS(sce, file=sprintf('%s/%s/%s/%s_allBarcodes.RDS', TASKDIR, PROJECT, SAMPLE, SAMPLE))

#plot barcode ranks
sce$isCell50 <- ifelse(sce$Barcode %in% eD50$Barcode, "YES", "NO")
sce$isCell75 <- ifelse(sce$Barcode %in% eD75$Barcode, "YES", "NO")
sceranks <- barcodeRanks(counts(sce))
sceranks <- data.frame(sceranks)
sceranks$isCell50 <- sce$isCell50
sceranks$isCell75 <- sce$isCell75
png(sprintf('%s/%s/%s/n50/%s_rankBarcodePlot_emptyDrops50.tiff', TASKDIR, PROJECT, SAMPLE, SAMPLE))
ggplot(data.frame(sceranks)) +
	ggtitle(paste(SAMPLE, 'n50', sep=" ")) +
	geom_point(aes(x=rank, y=total, color=isCell50), alpha=0.5) +
	scale_y_log10() +
	scale_x_log10() +
	xlab("Barcode Rank (log10)") +
	ylab("Total Counts (log10)") +
	theme_mary()	
dev.off()

png(sprintf('%s/%s/%s/n75/%s_rankBarcodePlot_emptyDrops75.tiff', TASKDIR, PROJECT, SAMPLE, SAMPLE))
ggplot(data.frame(sceranks)) +
	ggtitle(paste(SAMPLE, 'n75', sep=" ")) +
	geom_point(aes(x=rank, y=total, color=isCell75), alpha=0.5) +
	scale_y_log10() +
	scale_x_log10() +
	xlab("Barcode Rank (log10)") +
	ylab("Total Counts (log10)") +
	theme_mary()	
dev.off()

#append QC metrics
is.H1N1 <- grepl("^H1N1_", rowData(sce)$Symbol)
is.H3N2 <- grepl("^H3N2_", rowData(sce)$Symbol)

eD50 <- addPerCellQC(eD50, subsets=list(H1N1=is.H1N1, H3N2=is.H3N2))
eD75 <- addPerCellQC(eD75, subsets=list(H1N1=is.H1N1, H3N2=is.H3N2))

#normalize count matrices
eD50 <- logNormCounts(eD50)
eD75 <- logNormCounts(eD75)

#flag barcodes with a mix of H1N1 and H3N2 by simple threshold 20/80%
eD50$call <- ifelse(eD50$subsets_H1N1_percent > 80, "H1N1", ifelse(eD50$subsets_H1N1_percent < 20, "H3N2", "flag"))
eD75$call <- ifelse(eD75$subsets_H1N1_percent > 80, "H1N1", ifelse(eD75$subsets_H1N1_percent < 20, "H3N2", "flag"))

#interpret the barcodes
A=read.table(sprintf("%s/A_whitelist.txt", REFDIR), header=F, stringsAsFactors=F)
B=read.table(sprintf("%s/B_whitelist.txt", REFDIR), header=F, stringsAsFactors=F)
C=read.table(sprintf("%s/C_whitelist.txt", REFDIR), header=F, stringsAsFactors=F)
D=read.table(sprintf("%s/D_whitelist.txt", REFDIR), header=F, stringsAsFactors=F)

A$ind <- paste0("A", 1:96)
B$ind <- paste0("B", 1:96)
C$ind <- paste0("C", 1:96)
D$ind <- paste0("D", 1:96)

#convert the single cell objects to data frames
merge_cnts_with_meta <- function(sce) {
	#extract meta info
	df <- data.frame(colData(sce))
	
	#decode the barcodes
	df <- separate(df, Barcode, into=c("A", "B", "C", "D"), remove=F)
	df$A <- A[match(df$A, A$V1), 'ind']
	df$B <- B[match(df$B, B$V1), 'ind']
	df$C <- C[match(df$C, C$V1), 'ind']
	df$D <- D[match(df$D, D$V1), 'ind']
	df$Index <- paste(df$A, df$B, df$C, df$D, sep="_")
	rownames(df) <- df$Index

	#merge counts into dataframe
	cnts <- data.frame(t(as.matrix(counts(sce))))
	names(cnts) <- paste(names(cnts), "raw", sep="_")
	df <- cbind(df, cnts)
	lcnts <- data.frame(t(as.matrix(logcounts(sce))))
	names(lcnts) <- paste(names(lcnts), "norm", sep="_")
	df <- cbind(df, lcnts)
	
	#calculate the umap coordinates
  	set.seed(1000)
  	um <- umap(df[,grepl('_norm', colnames(df))])
  	um <- data.frame(um)
  	df <- cbind(df, um)
	
	return(df)
}

df50 <- merge_cnts_with_meta(eD50)
df75 <- merge_cnts_with_meta(eD75)

#compare segment ratio *** RETIRED in favor of a more sophisticated approach
#call_segment <- function(df) {
#  df$PB2 <- ifelse(df$H1N1_PB2_norm == df$H3N2_PB2_norm, "NA", 
#                         ifelse(df$H1N1_PB2_norm > df$H3N2_PB2_norm, "H1N1", "H3N2"))
#  df$PB1 <- ifelse(df$H1N1_PB1_norm == df$H3N2_PB1_norm, "NA", 
#                         ifelse(df$H1N1_PB1_norm > df$H3N2_PB1_norm, "H1N1", "H3N2"))
#  df$PA <- ifelse(df$H1N1_PA_norm == df$H3N2_PA_norm, "NA", 
#                        ifelse(df$H1N1_PA_norm > df$H3N2_PA_norm, "H1N1", "H3N2"))
#  df$HA <- ifelse(df$H1N1_HA_norm == df$H3N2_HA_norm, "NA", 
#                        ifelse(df$H1N1_HA_norm > df$H3N2_HA_norm, "H1N1", "H3N2"))
#  df$NP <- ifelse(df$H1N1_NP_norm == df$H3N2_NP_norm, "NA", 
#                        ifelse(df$H1N1_NP_norm > df$H3N2_NP_norm, "H1N1", "H3N2"))
#  df$NA6 <- ifelse(df$H1N1_NA_norm == df$H3N2_NA_norm, "NA", 
#                         ifelse(df$H1N1_NA_norm > df$H3N2_NA_norm, "H1N1", "H3N2"))
#  df$M <- ifelse(df$H1N1_M_norm == df$H3N2_M_norm, "NA", 
#                       ifelse(df$H1N1_M_norm > df$H3N2_M_norm, "H1N1", "H3N2"))
#  df$NS <- ifelse(df$H1N1_NS_norm == df$H3N2_NS_norm, "NA", 
#                        ifelse(df$H1N1_NS_norm > df$H3N2_NS_norm, "H1N1", "H3N2"))
#  df$H1N1_segs <- rowSums(df[,c("PB2", "PB1", "PA", "HA", "NP", "NA6", "M", "NS")]=="H1N1")
#  df$H3N2_segs <- rowSums(df[,c("PB2", "PB1", "PA", "HA", "NP", "NA6", "M", "NS")]=="H3N2")
#  df$Indistinguishable_segs <- rowSums(df[,c("PB2", "PB1", "PA", "HA", "NP", "NA6", "M", "NS")]=="NA")
#  return(df)
#}

#eD50df <- call_segment(eD50df)
#eD75df <- call_segment(eD75df)

#Prep data in for SoupX - estimate the contamination rate
#library(SoupX)
#estimate_soup <- function(sceALL, sceSub, th) { #***RETIRED - merged and expanded with next function
#	tmpTOD <- counts(sceALL)
#	tmpTOC <- counts(sceSub)
#	colnames(tmpTOC) <- sceSub$Barcode
#	colnames(tmpTOD) <- sceALL$Barcode
#	SC <- SoupChannel(tod=tmpTOD, toc=tmpTOC, metaData = NULL, calcSoupProfile = FALSE)
#	SC <- estimateSoup(SC, soupRange=c(0,th), keepDroplets=FALSE)
#
#	H1N1.genes <- c("H1N1_PB2", "H1N1_PB1", "H1N1_PA", "H1N1_NP", "H1N1_M", "H1N1_NS", 
#				"H1N1_HA", "H1N1_NA")
#	H3N2.genes <- c("H3N2_PB2", "H3N2_PB1", "H3N2_PA", "H3N2_NP", "H3N2_M", "H3N2_NS", 
#				"H3N2_HA", "H3N2_NA")
#
#	#estimate contaminant fraction
#	toUseSoup = data.frame(colData(sceSub)[,c('Barcode', 'call')])
#	rownames(toUseSoup) <- toUseSoup$Barcode
#	toUseSoup$H1N1 <- ifelse(toUseSoup$call=="H1N1", TRUE, FALSE)
#	toUseSoup$H3N2 <- ifelse(toUseSoup$call=="H3N2", TRUE, FALSE)
#	test <- as.matrix(toUseSoup[,c("H1N1", "H3N2")])
#
#	SC = calculateContaminationFraction(SC, list(H1N1 = H3N2.genes, H3N2 = H1N1.genes), test)
#	print(unique(SC$metaData$rho))
	
	#***RETIRED in favor of obtaining p-values 
	#out = adjustCounts(SC)
	
	#sxcnts <- data.frame(t(as.matrix(out)))
	#names(sxcnts) <- paste(names(sxcnts), "sxAdj_raw", sep="_")
	
	#sxN = SingleCellExperiment(assays=list(counts=out))
	#sxN = logNormCounts(sxN)
	#SXcntsN <- data.frame(t(as.matrix(logcounts(sxN))))
	#names(SXcntsN) <- paste(names(SXcntsN), "sxAdj_norm", sep="_")

	#cnts <- cbind(sxcnts, SXcntsN)
	#return(cnts)
#}

#eDdf <- cbind(eDdf, eDSX)

#call_segment_sxAdj <- function(df) { #***RETIRED in favor of p-values
#	df$PB2_sxAdj <- ifelse(df$H1N1_PB2_sxAdj_norm == df$H3N2_PB2_sxAdj_norm, "NA", 
#						ifelse(df$H1N1_PB2_sxAdj_norm > df$H3N2_PB2_sxAdj_norm, "H1N1", "H3N2"))
#	df$PB1_sxAdj <- ifelse(df$H1N1_PB1_sxAdj_norm == df$H3N2_PB1_sxAdj_norm, "NA", 
#						ifelse(df$H1N1_PB1_sxAdj_norm > df$H3N2_PB1_sxAdj_norm, "H1N1", "H3N2"))
#	df$PA_sxAdj <- ifelse(df$H1N1_PA_sxAdj_norm == df$H3N2_PA_sxAdj_norm, "NA", 
#						ifelse(df$H1N1_PA_sxAdj_norm > df$H3N2_PA_sxAdj_norm, "H1N1", "H3N2"))
#	df$HA_sxAdj <- ifelse(df$H1N1_HA_sxAdj_norm == df$H3N2_HA_sxAdj_norm, "NA", 
#						ifelse(df$H1N1_HA_sxAdj_norm > df$H3N2_HA_sxAdj_norm, "H1N1", "H3N2"))
#	df$NP_sxAdj <- ifelse(df$H1N1_NP_sxAdj_norm == df$H3N2_NP_sxAdj_norm, "NA", 
#						ifelse(df$H1N1_NP_sxAdj_norm > df$H3N2_NP_sxAdj_norm, "H1N1", "H3N2"))
#	df$NA6_sxAdj <- ifelse(df$H1N1_NA_sxAdj_norm == df$H3N2_NA_sxAdj_norm, "NA", 
#						ifelse(df$H1N1_NA_sxAdj_norm > df$H3N2_NA_sxAdj_norm, "H1N1", "H3N2"))
#	df$M_sxAdj <- ifelse(df$H1N1_M_sxAdj_norm == df$H3N2_M_sxAdj_norm, "NA", 
#						ifelse(df$H1N1_M_sxAdj_norm > df$H3N2_M_sxAdj_norm, "H1N1", "H3N2"))
#	df$NS_sxAdj <- ifelse(df$H1N1_NS_sxAdj_norm == df$H3N2_NS_sxAdj_norm, "NA", 
#						ifelse(df$H1N1_NS_sxAdj_norm > df$H3N2_NS_sxAdj_norm, "H1N1", "H3N2"))
#	df$H1N1_segs_sxAdj <- rowSums(df[,c("PB2_sxAdj", "PB1_sxAdj", "PA_sxAdj", "HA_sxAdj", "NP_sxAdj", "NA6_sxAdj", "M_sxAdj", "NS_sxAdj")]=="H1N1")
#	df$H3N2_segs_sxAdj <- rowSums(df[,c("PB2_sxAdj", "PB1_sxAdj", "PA_sxAdj", "HA_sxAdj", "NP_sxAdj", "NA6_sxAdj", "M_sxAdj", "NS_sxAdj")]=="H3N2")
#	df$Indistinguishable_segs_sxAdj <- rowSums(df[,c("PB2_sxAdj", "PB1_sxAdj", "PA_sxAdj", "HA_sxAdj", "NP_sxAdj", "NA6_sxAdj", "M_sxAdj", "NS_sxAdj")]=="NA")
#	return(df)
#}

#eDdf <- call_segment_sxAdj(eDdf)

#Modify the function to return the pvalues rather than TRUE/FALSE
estimateNonExpressingCellsMod <- function (sc, nonExpressedGeneList, clusters = NULL, maximumContamination = 1,
                                           pCut = 0.05)
{ 
  if (!is(sc, "SoupChannel")) 
    stop("sc is not a valid SoupChannel object")
  if (is.null(clusters)) {
    if ("clusters" %in% colnames(sc$metaData)) {
      clusters = setNames(as.character(sc$metaData$clusters),
                          rownames(sc$metaData))
    }
  }
  if (is.null(clusters) || (length(clusters) == 1 && clusters ==
                            FALSE)) {
    message("No clusters found or supplied, using every cell as its own cluster.")
    clusters = setNames(rownames(sc$metaData), rownames(sc$metaData))
  }
  if (!all(colnames(sc$toc) %in% names(clusters))) 
    stop("Invalid cluster specification.  clusters must be a named vector with all column names in the table of counts appearing.")
  if (!is.list(nonExpressedGeneList)) 
    stop("nonExpressedGeneList must be a list of sets of genes.  e.g. list(HB = c('HBB','HBA2'))")
  tgtGns = unique(unlist(nonExpressedGeneList))
  dat = sc$toc[tgtGns, , drop = FALSE]
  cnts = do.call(rbind, lapply(nonExpressedGeneList, function(e) colSums(dat[e, 
                                                                             , drop = FALSE])))
  exp = outer(sc$soupProfile[tgtGns, "est"], sc$metaData$nUMIs *
                maximumContamination)
  colnames(exp) = colnames(cnts)
  rownames(exp) = tgtGns
  exp = do.call(rbind, lapply(nonExpressedGeneList, function(e) colSums(exp[e, 
                                                                            , drop = FALSE])))
  s = split(names(clusters), clusters)
  clustExp = ppois(cnts - 1, exp, lower.tail = FALSE)
  print(dim(clustExp))
  clustExp = do.call(rbind, lapply(s, function(e) apply(clustExp[,  
                                                                 e, drop = FALSE], 1, min)))
  #clustExp = clustExp >= pCut
  clustExp = clustExp[match(clusters, rownames(clustExp)),
                      , drop = FALSE]
  rownames(clustExp) = names(clusters)
  #if (sum(clustExp) == 0) {
  #  warning("No non-expressing cells identified.  Consider setting clusters=FALSE, increasing maximumContamination and/or pCut.")
  #}
  #if (sum(clustExp) > 0 && sum(clustExp) < 100) {
  #  warning("Fewer than 100 non-expressing cells identified.  The estimation of the contamination fraction may be inaccurate.  Consider setting clusters=FALSE, #increasing maximumContamination and/or pCut.")
  #}
  return(clustExp)

}

execute_soupX_pvals <- function(sceALL, sceSub, umi) {
    tmpTOD <- counts(sceALL)
    tmpTOC <- counts(sceSub)
    colnames(tmpTOC) <- sceSub$Barcode
    colnames(tmpTOD) <- sceALL$Barcode
    SC <- SoupChannel(tod=tmpTOD, toc=tmpTOC, metaData = NULL, calcSoupProfile = FALSE)
    SC <- estimateSoup(SC, soupRange=c(0,umi), keepDroplets=FALSE)
    H1N1.genes <- c("H1N1_PB2", "H1N1_PB1", "H1N1_PA", "H1N1_NP", "H1N1_M", "H1N1_NS", 
				"H1N1_HA", "H1N1_NA")
    H3N2.genes <- c("H3N2_PB2", "H3N2_PB1", "H3N2_PA", "H3N2_NP", "H3N2_M", "H3N2_NS", 
				"H3N2_HA", "H3N2_NA")

    #estimate contaminant fraction
    toUseSoup = data.frame(colData(sceSub)[,c('Barcode', 'call')])
    rownames(toUseSoup) <- toUseSoup$Barcode
    toUseSoup$H1N1 <- ifelse(toUseSoup$call=="H1N1", TRUE, FALSE)
    toUseSoup$H3N2 <- ifelse(toUseSoup$call=="H3N2", TRUE, FALSE)
    test <- as.matrix(toUseSoup[,c("H1N1", "H3N2")])

    SC = calculateContaminationFraction(SC, list(H1N1 = H3N2.genes, H3N2 = H1N1.genes), test)
    print(unique(SC$metaData$rho))

    mC = unique(SC$metaData$rho)
    postBool <- data.frame(estimateNonExpressingCells(SC, nonExpressedGeneList=list(H1N1 = H3N2.genes, H3N2 = H1N1.genes), maximumContamination = mC, FDR=0.01))
    postBool$H1N1only <- ifelse(postBool$H1N1 & !postBool$H3N2, TRUE, FALSE)
    postBool$H3N2only <- ifelse(postBool$H3N2 & !postBool$H1N1, TRUE, FALSE)
    postBool <- as.matrix(postBool[,c("H1N1only", "H3N2only")])
    SC = calculateContaminationFraction(SC, list(H1N1only = H3N2.genes, H3N2only = H1N1.genes), postBool)
    print(unique(SC$metaData$rho))

    pvalsLibList=list()
    for ( f in seq(1, 50, 1)) {
            pvals = estimateNonExpressingCellsMod(SC, nonExpressedGeneList = list(H1N1 = H1N1.genes, H3N2=H3N2.genes, 
        H1N1_PB2 = "H1N1_PB2", 
        H1N1_PB1 = "H1N1_PB1", 
        H1N1_PA = "H1N1_PA", 
        H1N1_NP = "H1N1_NP", 
        H1N1_M = "H1N1_M", 
        H1N1_NS = "H1N1_NS",
        H1N1_HA = "H1N1_HA", 
        H1N1_NA = "H1N1_NA",
        H3N2_PB2 = "H3N2_PB2", 
        H3N2_PB1 = "H3N2_PB1", 
        H3N2_PA = "H3N2_PA", 
        H3N2_NP = "H3N2_NP", 
        H3N2_M = "H3N2_M", 
        H3N2_NS = "H3N2_NS",
        H3N2_HA = "H3N2_HA", 
        H3N2_NA = "H3N2_NA"), maximumContamination = f/100)
        pvalsLibList[[f]]=pvals
    }
    
    return(pvalsLibList)
}

tt50 <- execute_soupX_pvals(sce, eD50, 50)
tt75 <- execute_soupX_pvals(sce, eD75, 75)

#save the pvalue lists for future use
saveRDS(tt50, file=sprintf('%s/%s/%s/n50/%s_50counts_pvalsList_1-50.RDS', TASKDIR, PROJECT, SAMPLE, SAMPLE))
saveRDS(tt75, file=sprintf('%s/%s/%s/n75/%s_75counts_pvalsList_1-50.RDS', TASKDIR, PROJECT, SAMPLE, SAMPLE))

#Cutom function to extract and manipulate the data
return_DF <- function(meta, pvalsList, cont, SAMPLE, th, COUNTS) { 
  
  #extract the pvals for the given contamination % and adjust the p-values
  df = data.frame(apply(pvalsList[[cont]], 2, function(x) p.adjust(x, method="fdr"))) #switched to FDR *less stringent than the bonferroni 
  dfx = data.frame(apply(df, 2, function(x) x < (th/100))) #use an FDR of 1% or 5%
  dfx$H1N1_TRUE_segs <- rowSums(dfx[,grepl("^H1N1_", colnames(dfx))]) #count number of H1N1 segs
  dfx$H3N2_TRUE_segs <- rowSums(dfx[,grepl("^H3N2_", colnames(dfx))]) #count number of H3N2 segs
  
  #for each segment, determine if there is statistical evidence for expression of the H1N1, H3N2, both, or neither 
  dfx$PB2_con <- factor(paste(dfx$H1N1_PB2, dfx$H3N2_PB2, sep="_"))
  levels(dfx$PB2_con) <- list(NEITHER="FALSE_FALSE", H3N2="FALSE_TRUE", H1N1="TRUE_FALSE", BOTH="TRUE_TRUE")
  
  dfx$PB1_con <- factor(paste(dfx$H1N1_PB1, dfx$H3N2_PB1, sep="_"))
  levels(dfx$PB1_con) <- list(NEITHER="FALSE_FALSE", H3N2="FALSE_TRUE", H1N1="TRUE_FALSE", BOTH="TRUE_TRUE")
  
  dfx$PA_con <- factor(paste(dfx$H1N1_PA, dfx$H3N2_PA, sep="_"))
  levels(dfx$PA_con) <- list(NEITHER="FALSE_FALSE", H3N2="FALSE_TRUE", H1N1="TRUE_FALSE", BOTH="TRUE_TRUE")
  
  dfx$HA_con <- factor(paste(dfx$H1N1_HA, dfx$H3N2_HA, sep="_"))
  levels(dfx$HA_con) <- list(NEITHER="FALSE_FALSE", H3N2="FALSE_TRUE", H1N1="TRUE_FALSE", BOTH="TRUE_TRUE")
  
  dfx$NP_con <- factor(paste(dfx$H1N1_NP, dfx$H3N2_NP, sep="_"))
  levels(dfx$NP_con) <- list(NEITHER="FALSE_FALSE", H3N2="FALSE_TRUE", H1N1="TRUE_FALSE", BOTH="TRUE_TRUE")
  
  dfx$NA_con <- factor(paste(dfx$H1N1_NA, dfx$H3N2_NA, sep="_"))
  levels(dfx$NA_con) <- list(NEITHER="FALSE_FALSE", H3N2="FALSE_TRUE", H1N1="TRUE_FALSE", BOTH="TRUE_TRUE")
  
  dfx$M_con <- factor(paste(dfx$H1N1_M, dfx$H3N2_M, sep="_"))
  levels(dfx$M_con) <- list(NEITHER="FALSE_FALSE", H3N2="FALSE_TRUE", H1N1="TRUE_FALSE", BOTH="TRUE_TRUE")
  
  dfx$NS_con <- factor(paste(dfx$H1N1_NS, dfx$H3N2_NS, sep="_"))
  levels(dfx$NS_con) <- list(NEITHER="FALSE_FALSE", H3N2="FALSE_TRUE", H1N1="TRUE_FALSE", BOTH="TRUE_TRUE")
  
  #
  dfx$dup <- apply(dfx[,grepl("_con", colnames(dfx))], 1, function(x) "BOTH" %in% x)
  dfx$dup_no <- apply(dfx[,grepl("_con", colnames(dfx))], 1, function(x) sum(x =="BOTH"))
 
  print("Calculating doublet rate: ")
  print(sum(dfx$dup)/length(dfx[,1]))

  dfx$virus <- ifelse(dfx$H1N1 & dfx$H3N2, 'BOTH',
                      ifelse(dfx$H1N1 & !(dfx$H3N2), 'H1N1',
                             ifelse(dfx$H3N2 & !(dfx$H1N1), 'H3N2', 'NEITHER')))
  
  print("Of the singlets, gene sets for the viral strains:")
  table(dfx[!dfx$dup, 'virus'])

  tmp <- data.frame(pvalsList[[cont]])
  names(tmp) <- paste(names(tmp), 'pval', sep='_')
  tmp$H1N1_segs_fdr <- dfx$H1N1_TRUE_segs
  tmp$H3N2_segs_fdr <- dfx$H3N2_TRUE_segs
  tmp$dup_segs_fdr <- dfx$dup_no
  tmp$virus <- dfx$virus
  
  tmp <- cbind(tmp, dfx[,grepl("_con", colnames(dfx))])  
  
  newmeta <- cbind(meta, tmp)
  
  plotDF <- newmeta[,c("Index", "dup_segs_fdr", colnames(newmeta)[grep("_norm", colnames(newmeta))])]
  plotDF <- melt(plotDF, id.vars=c("Index", "dup_segs_fdr"))
  plotDF <- tidyr::separate(plotDF, variable, into=c("strain", "segment", "type"), sep="_") 
  plotDF$type <- NULL
  plotDF <- dcast(plotDF, Index + dup_segs_fdr + segment ~ strain)
  plotDF$segment <- factor(plotDF$segment, levels=c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))

  dupPlot <- ggplot(plotDF) +
	  facet_wrap(~segment, ncol=2) +
	  geom_point(aes(x=H1N1, y=H3N2, color=dup_segs_fdr), alpha=0.7, size=0.5) +
    	  scale_color_gradientn(colors=c("#CCCCCC", "#403891ff", "#6b4596ff", "#90548bff", "#b8627dff", "#de7065ff", "#f68f46ff", "#f9b641ff", "#efe350ff")) +
    	  theme_mary()

  png(sprintf('%s/%s/%s/n%s/%s_%scounts_contamination%s_fdr%s_doublets.png', TASKDIR, PROJECT, SAMPLE, COUNTS, SAMPLE, COUNTS, cont, th), units='in', width=5, height=10, res=300)
  plot(dupPlot)
  dev.off()
  
  newmeta$consensus <- paste(newmeta$PB2_con, newmeta$PB1_con, newmeta$PA_con, newmeta$HA_con, newmeta$NP_con, newmeta$NA_con, newmeta$M_con, newmeta$NS_con, sep="_")
  
  newmeta$droplet <- ifelse(newmeta$dup_segs_fdr > 0, "DOUBLET", 
                  ifelse(newmeta$H1N1_segs_fdr == 0 & newmeta$H3N2_segs_fdr == 8, "complete_H3N2",
                  ifelse(newmeta$H1N1_segs_fdr == 8 & newmeta$H3N2_segs_fdr == 0, "complete_H1N1",
                  ifelse(newmeta$H1N1_segs_fdr == 0, "partial_H3N2",
                  ifelse(newmeta$H3N2_segs_fdr == 0, "partial_H1N1", "REASSORTMENT")))))
  table(newmeta$droplet)
  
  umapPlot <- ggplot(newmeta, aes(X1, X2)) +
    geom_point(aes(color=droplet), alpha=0.7,size=0.5) +
    scale_color_manual(values=c("#440154FF", "#FDE725FF", "#CCCCCC", "#29AF7FFF", "#B8DE29FF", "#9900CC")) +
    #scale_color_gradientn(colors=c("#CCCCCC", "#403891ff", "#6b4596ff", "#90548bff", "#b8627dff", "#de7065ff", "#f68f46ff", "#f9b641ff", "#efe350ff")) +
    theme_mary()
  
  png(sprintf('%s/%s/%s/n%s/%s_%scounts_contamination%s_fdr%s_umap.png', TASKDIR, PROJECT, SAMPLE, COUNTS, SAMPLE, COUNTS, cont, th), units='in', width=5, height=6, res=300)
  plot(umapPlot)
  dev.off()
  
  #RA <- filter(newmeta, droplet == "REASSORTMENT") 
  #RA <- RA[,c(grepl("_con", colnames(RA)))]
  #RAm <- RA
  #RAm$Index <- rownames(RAm)
  #RAm <- melt(RAm, id.vars=c("Index"))
  
  #mat <- matrix(as.numeric(c("H1N1" = 1, "H3N2" = 2, "NEITHER" = 0)[as.matrix(RA)]), ncol=8)
  #colnames(mat) <- colnames(RA)
  
  #png(sprintf('%s/210721/FDR%s/%s_contamination%s_fdr%s_reassortmentHeatmap.png', TASKDIR, th, nameString, cont, th), units='in', width=5, height=5, res=300)
  #draw(Heatmap(t(mat),
  #        cluster_rows=F,
  #        clustering_method_columns = "ward.D2",
  #        col=viridis(3)))
  #dev.off()

  newmeta$pos_segs_fdr <- newmeta$H1N1_segs_fdr + newmeta$H3N2_segs_fdr - newmeta$dup_segs_fdr
  
  saveRDS(newmeta, file=sprintf('%s/%s/%s/n%s/%s_%scounts_contamination%s_fdr%s.RDS', TASKDIR, PROJECT, SAMPLE, COUNTS, SAMPLE, COUNTS, cont, th))
  write.table(newmeta, file=sprintf('%s/%s/%s/n%s/%s_%scounts_contamination%s_fdr%s_data.txt', TASKDIR, PROJECT, SAMPLE, COUNTS, SAMPLE, COUNTS, cont, th), sep="\t", quote=F, col.names=T, row.names=F, eol="\n")
  
  reassortments <- filter(newmeta, droplet == "REASSORTMENT") %>% count(consensus) %>% data.frame() %>% arrange(desc(n)) 
  reassortments <- separate(reassortments, consensus, into=c("seg1", "seg2", "seg3", "seg4", "seg5", "seg6", "seg7", "seg8"), remove=F)
  reassortments$missing <- rowSums(reassortments[,c(grepl("seg", colnames(reassortments)))]=="NEITHER")
  
  write.table(reassortments, file=sprintf('%s/%s/%s/n%s/%s_%scounts_contamination%s_fdr%s_reassortments.txt', TASKDIR, PROJECT, SAMPLE, COUNTS, SAMPLE, COUNTS, cont, th), sep="\t", quote=F, col.names=T, row.names=F, eol="\n")
  
  #return(newmeta)
}

return_DF(df50, tt50, 5, SAMPLE, 1, 50)
return_DF(df75, tt75, 5, SAMPLE, 1, 75)
return_DF(df50, tt50, 10, SAMPLE, 1, 50)
return_DF(df75, tt75, 10, SAMPLE, 1, 75)
return_DF(df50, tt50, 15, SAMPLE, 1, 50)
return_DF(df75, tt75, 15, SAMPLE, 1, 75)
return_DF(df50, tt50, 20, SAMPLE, 1, 50)
return_DF(df75, tt75, 20, SAMPLE, 1, 75)

return_DF(df50, tt50, 5, SAMPLE, 5, 50)
return_DF(df75, tt75, 5, SAMPLE, 5, 75)
return_DF(df50, tt50, 10, SAMPLE, 5, 50)
return_DF(df75, tt75, 10, SAMPLE, 5, 75)
return_DF(df50, tt50, 15, SAMPLE, 5, 50)
return_DF(df75, tt75, 15, SAMPLE, 5, 75)
return_DF(df50, tt50, 20, SAMPLE, 5, 50)
return_DF(df75, tt75, 20, SAMPLE, 5, 75)
