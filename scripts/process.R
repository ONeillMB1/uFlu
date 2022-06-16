.libPaths("/pasteur/zeus/projets/p01/uFlu/reassortment_project/resources/R_libs/4.1.0")

##########################################################################################################
##### This script has been developed by Mary O'Neill in 2021-2022 for the uFlu project.
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

ALIGNDIR="/pasteur/zeus/projets/p01/uFlu/reassortment_project/02_STARsolo_outputs"
TASKDIR="/pasteur/zeus/projets/p01/uFlu/reassortment_project/03_custom_analyses"
REFDIR="/pasteur/zeus/projets/p01/uFlu/reassortment_project/resources/Ref_data"

args <- commandArgs(TRUE)
SAMPLE=args[1]
PROJECT=args[2]

set.seed(1000)

#make theme for figures
theme_mary <- function() {
  font <- "Arial"
  
  theme_light() + 
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          text = element_text(size=10),
          #axis.text.x=element_text(angle=45, hjust=1),
          legend.position="bottom",
          legend.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA))
}

#load the data selected by the empty drops method
eD50 <- read10xCounts(sprintf('%s/%s/%s/%s.Solo.out/Gene/emptyDrops_50counts', ALIGNDIR, PROJECT, SAMPLE, SAMPLE))

#read in whole dataset
sce=read10xCounts(sprintf('%s/%s/%s/%s.Solo.out/Gene/raw', ALIGNDIR, PROJECT, SAMPLE, SAMPLE)) #full dataset, very big

#filter out barcodes with no counts
cnts = colSums(counts(sce))
sce = sce[,cnts>0] 
saveRDS(sce, file=sprintf('%s/%s/%s/%s_allBarcodes.RDS', TASKDIR, PROJECT, SAMPLE, SAMPLE))

#plot barcode ranks
sce$isCell50 <- ifelse(sce$Barcode %in% eD50$Barcode, "YES", "NO")
sceranks <- barcodeRanks(counts(sce))
sceranks <- data.frame(sceranks)
sceranks$isCell50 <- sce$isCell50
png(sprintf('%s/%s/%s/%s_rankBarcodePlot_emptyDrops50.tiff', TASKDIR, PROJECT, SAMPLE, SAMPLE))
ggplot(data.frame(sceranks)) +
	ggtitle(paste(SAMPLE, 'n50', sep=" ")) +
	geom_point(aes(x=rank, y=total, color=isCell50), alpha=0.5) +
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

#normalize count matrices
eD50 <- logNormCounts(eD50)

#flag barcodes with a mix of H1N1 and H3N2 by simple threshold 20/80%
eD50$call <- ifelse(eD50$subsets_H1N1_percent > 80, "H1N1", ifelse(eD50$subsets_H1N1_percent < 20, "H3N2", "flag"))

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

#Modify the SoupX function to return the pvalues rather than TRUE/FALSE
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
    for ( f in c(seq(1, 50, 1),unique(SC$metaData$rho)*100)) {
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
        pvalsLibList[[paste0(f)]]=pvals
    }
    
    return(pvalsLibList)
}

tt50 <- execute_soupX_pvals(sce, eD50, 50)

#save the pvalue lists for future use
saveRDS(tt50, file=sprintf('%s/%s/%s/%s_50counts_pvalsList_1-50.RDS', TASKDIR, PROJECT, SAMPLE, SAMPLE))

#Newly added analysis
identical(df50$Barcode, rownames(tt50[[5]])) #check that they are in the same order

return_DF <- function(meta, pvalsList, cont, th, umi) { 
  
  meta <- meta[meta$sum >= umi, ]
  meta$umiTH <- umi
  meta$cont <- cont
  meta$fdrTH <- th
  
  #extract the pvals for the given contamination % and adjust the p-values
  df = data.frame(apply(pvalsList[[cont]][meta$Barcode,], 2, function(x) p.adjust(x, method="fdr"))) #switched to FDR *less stringent than the bonferroni 
  dfx = data.frame(apply(df, 2, function(x) x < (th/100))) #use an FDR of 1% or 5%
  
  #count the number of segments for each strain
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
  
  tmp <- data.frame(pvalsList[[cont]][meta$Barcode,])
  names(tmp) <- paste(names(tmp), 'pval', sep='_')
  tmp$H1N1_segs_fdr <- dfx$H1N1_TRUE_segs
  tmp$H3N2_segs_fdr <- dfx$H3N2_TRUE_segs
  tmp$dup_segs_fdr <- dfx$dup_no
  tmp$virus <- dfx$virus
  
  tmp <- cbind(tmp, dfx[,grepl("_con", colnames(dfx))])  
  
  newmeta <- cbind(meta, tmp)
  
  newmeta$consensus <- paste(newmeta$PB2_con, newmeta$PB1_con, newmeta$PA_con, newmeta$HA_con, newmeta$NP_con, newmeta$NA_con, newmeta$M_con, newmeta$NS_con, sep="_")
  
  newmeta$droplet <- ifelse(newmeta$dup_segs_fdr > 0, "DOUBLET", 
                     ifelse(newmeta$H1N1_segs_fdr == 0 & newmeta$H3N2_segs_fdr == 8, "complete_H3N2",
                     ifelse(newmeta$H1N1_segs_fdr == 8 & newmeta$H3N2_segs_fdr == 0, "complete_H1N1",
                     ifelse(newmeta$H1N1_segs_fdr == 0, "partial_H3N2",
                     ifelse(newmeta$H3N2_segs_fdr == 0, "partial_H1N1", "REASSORTMENT")))))
  table(newmeta$droplet)
  
  newmeta$pos_segs_fdr <- newmeta$H1N1_segs_fdr + newmeta$H3N2_segs_fdr - newmeta$dup_segs_fdr
  
  #Calculate a few other things
  newmeta$missing <- rowSums(newmeta[,c(grepl("_con", colnames(newmeta)))]=="NEITHER")
  newmeta$duplicated <- rowSums(newmeta[,c(grepl("_con", colnames(newmeta)))]=="BOTH")
  newmeta$background <- ifelse(newmeta$dup_segs_fdr > 0, "DOUBLET", 
                               ifelse(newmeta$H3N2_segs_fdr < 5 & newmeta$H1N1_segs_fdr < 5, "INDISTINGUISHABLE",
                                      ifelse(newmeta$H1N1_segs_fdr > newmeta$H3N2_segs_fdr, "H1N1", "H3N2")))
  
  #ATTENTION: note that droplet2 was dropped from previous developments!
              
		      
  return(newmeta)
  
}


paramList = list()
for(th in c(1,5)) {
	for(umi in c(50, 75, 100, 125, 150)) {
		for(cont in c(seq(1,10,1), seq(15,50, 5), as.numeric(names(tt50)[[51]]))) {
			paramList[[paste0(umi, "UMIs", "_", cont, "cont", "_", th, "FDR")]] = return_DF(df50, tt50, cont, th, umi)
		}
  }
}

master <- data.frame(rbindlist(paramList, idcol="params"))
master$exp <- paste(PROJECT, SAMPLE, sep="_")

saveRDS(master, file=sprintf('%s/%s/%s/%s_master_paramEval.RDS', TASKDIR, PROJECT, SAMPLE, SAMPLE))


print("Finished pipeline successfully")
