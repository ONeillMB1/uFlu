.libPaths("/pasteur/zeus/projets/p01/uFlu/reassortment_project/resources/R_libs/4.1.0")

##########################################################################################################
##### This script has been written by Mary O'Neill in January 2022 for the uFlu project.
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

###
master = readRDS(sprintf('%s/%s/%s/%s_master_paramEval.RDS', TASKDIR, PROJECT, SAMPLE, SAMPLE))

table(master$background)

#Update background rules
master$background <- ifelse(master$dup_segs_fdr > 0, "DOUBLET", 
                        ifelse(master$H3N2_segs_fdr < 5 & master$H1N1_segs_fdr < 5, "INDISTINGUISHABLE",
                               ifelse(master$H1N1_segs_fdr > master$H3N2_segs_fdr, "H1N1", "H3N2")))
table(master$background)                           

#Drop the old droplet2 column that was based on the old background definintion
master$droplet2 <- NULL

#Updated
makeup <- master %>% group_by(umiTH, cont, fdrTH) %>% count(paste(droplet, background)) %>% data.frame() %>%
  group_by(umiTH, cont, fdrTH) %>%
  mutate(Percent = 100*n/sum(n)) %>%
  data.frame()

tab1 <- master %>% 
  group_by(exp, umiTH, cont, fdrTH) %>% 
  count(droplet=paste(droplet, background, sep="_")) %>% 
  data.frame() %>% 
  reshape2::dcast(exp+umiTH+cont+fdrTH~droplet)

tab1[is.na(tab1)] <- 0

tab1$cell_containing <- rowSums(tab1[,c(5:14)])
tab1$single_cell <- tab1$cell_containing - tab1$DOUBLET_DOUBLET
tab1$doublet_per <- tab1$DOUBLET/tab1$cell_containing*100

tab2 <- master %>% 
  group_by(umiTH, cont, fdrTH) %>% 
  count(background) %>% 
  data.frame() %>% 
  reshape2::dcast(umiTH+cont+fdrTH~background)

tab1$H1N1_background <- tab2$H1N1
tab1$H3N2_background <- tab2$H3N2
tab1$indistinguishable_background <- tab2$INDISTINGUISHABLE

tab3 <- master %>% 
  group_by(umiTH, cont, fdrTH) %>% 
  count(test=pos_segs_fdr==8 & droplet=="REASSORTMENT") %>% 
  data.frame() %>% 
  reshape2::dcast(umiTH+cont+fdrTH~test)


tab_8seg <- master %>%
  filter(pos_segs_fdr==8 & droplet != "DOUBLET") %>%
  group_by(umiTH, cont, fdrTH) %>%
  count(droplet) %>%
  data.frame() %>%
  reshape2::dcast(umiTH+cont+fdrTH ~ droplet)              

tab_8seg$RR <- tab_8seg$REASSORTMENT / rowSums(tab_8seg[,3:5]) *100 
tab_8seg$H1H3 <- tab_8seg$complete_H1N1 / tab_8seg$complete_H3N2


tab1$eight_seg_reassortment <- tab3$`TRUE`

tab1$eight_segment <- tab1$complete_H1N1_H1N1 + tab1$complete_H3N2_H3N2 + tab1$eight_seg_reassortment
tab1$reassortment_rate_8seg <- tab1$eight_seg_reassortment / tab1$eight_segment * 100
tab1$reassortment_rate_global <- rowSums(tab1[,grepl("REASSORTMENT", colnames(tab1))])/tab1$single_cell*100

tab1$ratio_pure_8seg <- tab1$complete_H1N1_H1N1 / tab1$complete_H3N2_H3N2 
tab1$ratio_background_pure <- (rowSums(tab1[,c("complete_H1N1_H1N1", "partial_H1N1_H1N1")]) / rowSums(tab1[,c("complete_H3N2_H3N2", "partial_H3N2_H3N2")]))
tab1$ratio_background_withRA <- tab1$H1N1_background / tab1$H3N2_background 

write.table(tab1, file=sprintf('%s/%s/%s/%s_summaryStats_updated2022.txt', TASKDIR, PROJECT, SAMPLE, SAMPLE), sep="\t", quote=F, col.names=T, row.names=F, eol="\n")


#Plot some summary stats
library(viridis)


#function needed to find rho
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

rho=as.numeric(levels(as.factor(tab1$cont)))[!is.wholenumber(as.numeric(levels(as.factor(tab1$cont))))]

RR <- ggplot(tab1, aes(x=cont, y=reassortment_rate_8seg, shape=as.factor(fdrTH), color=umiTH)) +
  facet_wrap(~exp, nrow=1, scales="free_y") +
  geom_point() +
  scale_color_viridis(end=0.9, option="B") +
  theme_mary() +
  theme(legend.position="none") +
  geom_vline(xintercept=rho, linetype="dashed", color="grey") +
  xlab("Maximum Contamination Fraction") +
  ylab("Reassortment Rate (% single-cells)") +
  theme(aspect.ratio=1)

doublet <- ggplot(tab1, aes(x=cont, y=doublet_per, shape=as.factor(fdrTH), color=umiTH)) +
  facet_wrap(~exp, nrow=1) +
  geom_point() +
  scale_color_viridis(end=0.9, option="B") +
  theme_mary() +
  theme(legend.position="none") +
  geom_vline(xintercept=rho, linetype="dashed", color="grey") +
  xlab("Maximum Contamination Fraction") +
  ylab("Doublet Rate (% cell-containing droplets)") +
  theme(aspect.ratio=1)

droplets <- ggplot(tab1, aes(x=cont, y=cell_containing, shape=as.factor(fdrTH), color=umiTH)) +
  facet_wrap(~exp, nrow=1, scales="free_y") +
  geom_point() +
  scale_color_viridis(end=0.9, option="B") +
  theme_mary() +
  theme(legend.position="none") +
  geom_vline(xintercept=rho, linetype="dashed", color="grey") +
  xlab("Maximum Contamination Fraction") +
  ylab("Cell Containing Droplets") +
  theme(aspect.ratio=1)

singlets <- ggplot(tab1, aes(x=cont, y=single_cell, shape=as.factor(fdrTH), color=umiTH)) +
  facet_wrap(~exp, nrow=1, scales="free_y") +
  geom_point() +
  scale_color_viridis(end=0.9, option="B") +
  theme_mary() +
  theme(legend.position="none") +
  geom_vline(xintercept=rho, linetype="dashed", color="grey") +
  xlab("Maximum Contamination Fraction") +
  ylab("Single Cells Identified") +
  theme(aspect.ratio=1)

rat <- ggplot(tab1, aes(x=cont, y=ratio_background_withRA, shape=as.factor(fdrTH), color=umiTH)) +
  facet_wrap(~exp, nrow=1) +
  geom_point() +
  scale_color_viridis(end=0.9, option="B") +
  theme_mary() +
  theme(legend.position="bottom") +
  geom_vline(xintercept=rho, linetype="dashed", color="grey") +
  geom_hline(yintercept=1, linetype="solid", color="grey") +
  xlab("Maximum Contamination Fraction") +
  ylab("Ratio of H1N1 / H3N2") +
  theme(aspect.ratio=1)

sumPlotHor <- plot_grid(droplets, doublet, singlets, RR, rat, nrow=1, align='hv', axis='tblr')
sumPlotVer <- plot_grid(droplets, doublet, singlets, RR, rat, ncol=1, align='hv', axis='tblr')

pdf(sprintf('%s/%s/%s/%s_paramEval_horizontal_updated2022.pdf', TASKDIR, PROJECT, SAMPLE, SAMPLE), height=4, width=15)
sumPlotHor
dev.off()

pdf(sprintf('%s/%s/%s/%s_paramEval_vertical_updated2022.pdf', TASKDIR, PROJECT, SAMPLE, SAMPLE), height=16, width=4)
sumPlotVer
dev.off()

#Correlations between 8 segments
RA <- filter(master, droplet == "REASSORTMENT" & pos_segs_fdr == 8)

RA <- RA %>% 
  group_by(umiTH, cont, fdrTH, exp, params) %>%
  count(consensus) %>% 
  data.frame() %>% 
  group_by(umiTH, cont, fdrTH, exp, params) %>%
  mutate(Percent = 100*n/sum(n)) %>% 
  arrange(desc(Percent)) %>%
  data.frame()

RA <- arrange(RA, fdrTH, umiTH, cont)

RA.short <- reshape2::dcast(RA[, c(1:3,6,8)], fdrTH+umiTH+cont~consensus)
RA.short[is.na(RA.short)] <- 0
rownames(RA.short) <- paste('fdr', RA.short$fdrTH, '_umi', RA.short$umiTH, '_cont', RA.short$cont, sep="")

x <- cor(t(as.matrix(RA.short[,4:length(RA.short)])))

x <- data.frame(x)
x$comp1 <- rownames(x)
xm <- reshape2::melt(x, id.vars=c("comp1"))
xm$comp1 <- factor(xm$comp1, levels=levels(xm$variable))

corTile <- ggplot(xm, aes(x=comp1, y=variable, fill=value)) +
  geom_tile() +
  theme_mary() +
  theme(axis.text.x=element_text(angle=90, hjust=1),
        aspect.ratio=1) +
  scale_x_discrete(labels=gsub("cont", "", sapply(strsplit(as.character(xm$comp1), "_"), `[[`, 3, simplify=TRUE))) +
  scale_y_discrete(labels=gsub("cont", "", sapply(strsplit(as.character(levels(xm$variable)), "_"), `[[`, 3, simplify=TRUE))) +
  scale_fill_gradientn(limits=c(0,1), colors=viridis(100)) +
  xlab("") +
  ylab("")

pdf(sprintf('%s/%s/%s/%s_corr8seg_updated2022.pdf', TASKDIR, PROJECT, SAMPLE, SAMPLE), height=28, width=30)
corTile
dev.off()

print("Finished pipeline successfully")

