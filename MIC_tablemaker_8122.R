library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(stringr)
library(scales)
library(hash)
library(epitools)
library(extrafont)
library(vcd)
library(rcompanion)

ecoffs <- read.csv('/Users/joshuacarter/Desktop/paper-cryptic-mic/ecoffs.csv')
drug_list <- c('RIF','RFB','INH','ETH','EMB','KAN','AMI','BDQ','CFZ','LEV','LZD','MXF','DLM') #
flag <- 0

#read in all non_int data to make Table S1
for (drug in drug_list) { 
  
  #read in input data file prepared in Jupyter to get mutation names
  dat <- read.csv(file = paste('/Users/joshuacarter/Desktop/cryptic-analysis/Interval_regression/Output_files/no_int/R/',drug,'/',drug,'_no_int_main_effects.csv',sep=''), header = TRUE)
  dat$DRUG = drug
  if (flag==0) {
    out <- dat
  } else {out <- rbind(out,dat)}
  flag <- 1
  
}

homoplasy <- read.csv(file='/Users/joshuacarter/Desktop/CRyPTIC_GPI_homoplasy.csv', header=TRUE)
out <- merge(out, homoplasy, on='MUTS', all.x = TRUE)
out[is.na(out$HOMOPLASY.COUNT),"HOMOPLASY.COUNT"] <- 0
out$HOMOPLASY2 <- ifelse(out$HOMOPLASY.COUNT>0, 'True', 'False')

#print total number of mutations tested
length(out[!(out$GENE %in% c('DRUG','LINEAGE','SITEID')),'MUTS'])

#print total number of significant mutations
length(out[out$COEF>0 & out$BH_PVAL<0.05 & !(out$GENE %in% c('DRUG','LINEAGE','SITEID')),'MUTS'])
#print number sig genes
length(unique(out[out$BH_PVAL<0.05 & !(out$GENE %in% c('DRUG','LINEAGE','SITEID')),'GENE']))
#print total genes tested
length(unique(out[!(out$GENE %in% c('DRUG','LINEAGE','SITEID')),'GENE']))
#list number of significant, MIC-elevating mutations by drug
out[out$BH_PVAL<0.05 & out$COEF>0 & !(out$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>%
  group_by(DRUG) %>%
  summarise(no_rows = length(MUTS))

#list number of significant, MIC-elevating mutations by drug
out[out$BH_PVAL<0.05 & (out$GENE %in% c('SITEID')),] %>%
  group_by(DRUG) %>%
  summarise(min = min(COEF), mean = median(COEF), max=max(COEF))

#list number of significant, MIC-elevating mutations by drug
out[out$BH_PVAL<0.05 & (out$GENE %in% c('LINEAGE')),] %>%
  group_by(DRUG) %>%
  summarise(min = min(COEF), mean = median(COEF), max=max(COEF))

#arrange dataframe by MICs
out[out$COEF>0 & out$BH_PVAL<0.05 & !(out$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>% arrange(desc(COEF))

#print max un-adjusted pval that passes FDR 0.05 by drug
out[out$BH_PVAL<0.05,] %>%
  group_by(DRUG) %>%
  summarise(no_rows = max(P))

#shift homoplasy to T/F
#out$HOMOPLASY <- out$HOMOPLASY2

#Stratum Fisher exact tests for odds of homoplasy and MIC-elevating status
mh_vect <- c()
for (drug in drug_list) {
  mh_vect <- append(mh_vect, c(length(out[out$DRUG==drug & out$HOMOPLASY.COUNT==0 & !(out$BH_PVAL<0.05 & out$COEF>0),'MUTS']),length(out[out$DRUG==drug & out$HOMOPLASY.COUNT>0 & !(out$BH_PVAL<0.05 & out$COEF>0),'MUTS']),length(out[out$DRUG==drug & out$HOMOPLASY.COUNT==0 & (out$BH_PVAL<0.05 & out$COEF>0),'MUTS']),length(out[out$DRUG==drug & out$HOMOPLASY.COUNT>0 & (out$BH_PVAL<0.05 & out$COEF>0),'MUTS'])), after = length(mh_vect))
}
homo <- array(mh_vect,
              dim = c(2,2,13),
              dimnames = list(
              Homplasy = c("False", "True"),
              MIC = c("SUS", "RES"),
              Penicillin.Level = drug_list))

woolf_test(homo)

for (i in c(1:13)) {
  print(drug_list[i])
  tmp <- fisher.test(homo[,,i])
  print(tmp$estimate)
  print(tmp$p.value)}

groupwiseCMH(homo,
             group   = 3,
             fisher  = TRUE,
             gtest   = FALSE,
             chisq   = FALSE,
             method  = "fdr",
             correct = "none",
             digits  = 3)

#test for homoplastic mutations having greater effects on MIC using Wilcoxon rank-sum
for (drug in drug_list) {
  test <- wilcox.test(out[out$HOMOPLASY.COUNT>0 & out$DRUG==drug,'COEF'],out[out$HOMOPLASY.COUNT==0 & out$DRUG==drug,'COEF'], alternative = 'greater',conf.int = TRUE)
  print(c(drug,test$p.value))
}

blah <- out[out$BH_PVAL<0.05 & !is.na(out$HOMOPLASY.COUNT) & !out$GENE %in% c('DRUG','LINEAGE','SITEID'),]
ggplot(data=blah, aes(x=blah$HOMOPLASY.COUNT, y=blah$COEF, color=blah$DRUG)) +
  geom_point() + scale_x_log10() + 
  theme_classic() + facet_grid(.~DRUG)


#Investigate the effects of promoter mutations on MIC compared to gene body mutations
out$PROM <- 0
out[out$POSITION<0,'PROM'] = 1
for (drug in drug_list) {
tmp <- lm(COEF ~ PROM + GENE, data=out[out$DRUG==drug & !(out$GENE %in% c('DRUG','LINEAGE','SITEID')),], weights = COUNTS)
print(drug)
if (summary(tmp)$coefficients[,4][2]<0.05) {
print(tmp$coefficients[2])
print(summary(tmp)$coefficients[,4][2])
}
}

#Stratum Fisher exact tests for odds of promoter and MIC-elevating status
mh_vect <- c()
for (drug in drug_list) {
  mh_vect <- append(mh_vect, c(length(out[out$DRUG==drug & out$PROM=='0' & !(out$BH_PVAL<0.05 & out$COEF>0),'MUTS']),length(out[out$DRUG==drug & out$PROM==1 & !(out$BH_PVAL<0.05 & out$COEF>0),'MUTS']),length(out[out$DRUG==drug & out$PROM==0 & (out$BH_PVAL<0.05 & out$COEF>0),'MUTS']),length(out[out$DRUG==drug & out$PROM==1 & (out$BH_PVAL<0.05 & out$COEF>0),'MUTS'])), after = length(mh_vect))
}
homo <- array(mh_vect,
              dim = c(2,2,13),
              dimnames = list(
                Homplasy = c("False", "True"),
                MIC = c("SUS", "RES"),
                Penicillin.Level = drug_list))

woolf_test(homo)

mantelhaen.test(homo)

for (i in c(1:13)) {
  print(drug_list[i])
  tmp <- fisher.test(homo[,,i])
  print(tmp$estimate)
  print(tmp$p.value)}

groupwiseCMH(homo,
             group   = 3,
             fisher  = TRUE,
             gtest   = FALSE,
             chisq   = FALSE,
             method  = "fdr",
             correct = "none",
             digits  = 3)

for (drug in drug_list) {
  if (length(out[out$BH_PVAL<0.05 & out$DRUG==drug & out$PROM==1 & !(out$GENE %in% c('DRUG','LINEAGE','SITEID')),'COEF'])>0) {
    test <- wilcox.test(out[out$BH_PVAL<0.05 & out$DRUG==drug & out$PROM==1 & !(out$GENE %in% c('DRUG','LINEAGE','SITEID')),'COEF'],out[out$BH_PVAL<0.05 & out$DRUG==drug & out$PROM==0 & !(out$GENE %in% c('DRUG','LINEAGE','SITEID')),'COEF'],conf.int = TRUE)
    print(c(drug,test$p.value))
  }
}


out$mut_aa <- substr(out$MUTS,nchar(out$MUTS),nchar(out$MUTS))

png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig2/prom_body_facet.png', width = 1028, height =768)
plot_dat <- out[out$ci_hi-out$ci_lo<4 & out$COUNTS>2 & ((((out$COUNTS>2 & out$COEF>1) | out$BH_PVAL<0.05) & out$GENE %in% c('embA','fabG1')) | ((out$COUNTS>2 & out$COEF>1) | out$BH_PVAL<0.05) & out$GENE %in% c('embB','inhA')) & out$COEF>0 & out$DRUG %in% c('INH','EMB') & !(out$GENE %in% c('DRUG','LINEAGE','SITEID')),]
plot_dat[plot_dat$GENE=='fabG1','PROM'] = 'PROM'
plot_dat <- plot_dat %>% arrange(DRUG, COEF)
row.names(plot_dat) <- NULL
plot_dat$X <- as.numeric(row.names(plot_dat))
hline.data <- data.frame(z = c(prom_ecoffs[5],prom_ecoffs[1]), DRUG = c("EMB","INH"))
ggplot(data=plot_dat,aes(x=X,y=COEF, color=factor(PROM,c('SNP','PROM')))) +
  theme_classic() + scale_x_continuous() + facet_grid(.~DRUG, scales="free_x", space="free_x") +
  scale_color_manual(values = c( 'lightgrey','#de7569')) +
  geom_pointrange(aes(ymin=ci_lo, ymax=ci_hi, alpha=ifelse(BH_PVAL<0.05,1,0.5), fatten=sqrt(COUNTS)), size=0.3) +
  #geom_point(position=position_jitterdodge(jitter.width = .5)) +
  #annotate("segment", y = prom_ecoffs[1], yend = prom_ecoffs[1], x = 0.6, xend = 1.5,color='tan', linetype='dashed', size=1.5) +
  #annotate("segment", y = prom_ecoffs[5], yend = prom_ecoffs[5], x = 1.6, xend = 2.5,color='tan', linetype='dashed', size=1.5) +
  #scale_x_discrete() + 
  coord_cartesian(clip = "off") +
  geom_hline(aes(yintercept = z) , color='tan', linetype='dashed', size=1.5, hline.data) + 
  geom_hline(aes(yintercept = 0) , color='black', size=1.5) + 
  geom_text(aes(label=ifelse(COUNTS>2,MUTS,''),hjust=0, vjust=-0.5,angle=90)) +
  theme(text=element_text(size=18, family="Arial Rounded MT Bold", color='black'), legend.position = 'none', axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank())
dev.off()

png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig2/prom_positions.png', width = 1028, height =768)
hline.data <- data.frame(a = c(0.3,0.3), DRUG = c("EMB","INH"))
ggplot(data=out[out$GENE!='ahpC' & out$DRUG %in% c('AMI','KAN','INH','ETH','EMB') & out$POSITION<0 & out$mut_aa %in% c('a','c','t','g','l'),],aes(x=POSITION,y=COEF, shape=ifelse(BH_PVAL<0.05,'True','False'), color=factor(mut_aa), size=5)) +
  theme_classic() + facet_grid(rows=vars(DRUG), space='free_y') + labs(color = "Mutant Nuc", shape='P<0.05', alpha='None') +
  scale_color_discrete(name = "Mutation", labels = c("a", "c", "g", "indel", "t")) +
  annotate("rect", xmin = c(-15), xmax = c(-5), ymin = -Inf, ymax = Inf, alpha = 0.3, fill = "tan") +
  geom_rect(inherit.aes = FALSE, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha=a), hline.data) +
  geom_hline(aes(yintercept = 0) , color='black', size=1) + 
  geom_point(aes(alpha=ifelse(BH_PVAL<0.05,1,1))) + 
  #geom_hline(aes(yintercept = z) , color='tan', linetype='dashed', size=1.5, hline.data) + 
  scale_shape_manual(values = c('True' = 16, 'False' = 1)) + scale_alpha(guide = 'none') +
  theme(text=element_text(size=18, family="Arial Rounded MT Bold", color='black'), panel.spacing.y=unit(0, "lines"))
dev.off()

#Mantel-Haenszel test for weighted odds of presence around -10 TATA element and MIC-elevating status
tmp <- out[out$GENE!='ahpC' & out$POSITION<0 & out$mut_aa %in% c('a','c','t','g','l'),]
mh_vect <- c()
for (drug in c('AMI','KAN','INH','ETH','EMB')) {
  mh_vect <- append(mh_vect, c(length(tmp[tmp$DRUG==drug & !(tmp$POSITION>-15 & tmp$POSITION< -5) & tmp$BH_PVAL>=0.05,'MUTS']),length(tmp[tmp$DRUG==drug & !(tmp$POSITION>-15 & tmp$POSITION< -5) & (tmp$BH_PVAL<0.05 & tmp$COEF>0),'MUTS']),length(tmp[tmp$DRUG==drug & (tmp$POSITION>=-15 & tmp$POSITION<= -5) & tmp$BH_PVAL>=0.05,'MUTS']),length(tmp[tmp$DRUG==drug & (tmp$POSITION>=-15 & tmp$POSITION<= -5) & (tmp$BH_PVAL<0.05 & tmp$COEF>0),'MUTS'])), after = length(mh_vect))
}
homo <- array(mh_vect,
              dim = c(2,2,5),
              dimnames = list(
                Significant = c("False", "True"),
                TATA = c("No", "Yes"),
                DRUG.Level = c('AMI','KAN','INH','ETH','EMB')))

woolf_test(homo)
mantelhaen.test(homo)

#add in ECOFFs to Table S1
ecoffs$ECOFF <- 0 
for (drug in drug_list) {
  ecoffs[ecoffs$DRUG==drug,'ECOFF'] <- log2(as.numeric(ecoffs[ecoffs$DRUG==drug,'ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG==drug,'COEF'][1]
}
out <- merge(out,ecoffs,on='DRUG')
#check how many significant MIC-elevating mutaitons are below respective ECOFF by drug
out[out$BH_PVAL<0.05 & out$COEF>0 & out$ci_lo<out$ECOFF & !(out$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>%
  group_by(DRUG) %>%
  summarise(no_rows = length(GENE))
#print total number of significant MIC-elevating mutations that are below respective ECOFF
length(out[out$BH_PVAL<0.05 & out$COEF>0 & out$ci_lo<out$ECOFF & !(out$GENE %in% c('DRUG','LINEAGE','SITEID')),'MUTS'])

#create dataframe with significant mutations after FDR correction
sig <- out[out$BH_PVAL<0.05,]

#create single resistant dataframe
res <- sig[sig$COEF>0,]
res <- res[res$GENE!='DRUG' & res$GENE!='SITEID' & res$GENE!='LINEAGE',]
sing_res <- res
#sing_res <- res[-grep(':',res$MUTS),]

#test for difference in promoter vs gene body
sing_res$PROM <- 0
sing_res[sing_res$POSITION<0,'PROM'] = 1
sing_res$PROM <- factor(sing_res$PROM)
wilcox.test(sing_res[sing_res$PROM==1,'COEF'],sing_res[sing_res$PROM==0,'COEF'],alternative='less')
for (drug in c('AMI','KAN','ETH','INH','EMB')) {
  test <- wilcox.test(sing_res[sing_res$PROM==1 & sing_res$DRUG==drug,'COEF'],sing_res[sing_res$PROM==0 & sing_res$DRUG==drug,'COEF'], alternative = 'less')
  print(c(drug,test$p.value))
}  

#add in ecoffs for promoter drugs and print Fig 2B
prom_ecoffs <- c()
for (drug in c('INH','ETH','AMI','KAN','EMB')) {
  prom_ecoffs <- append(prom_ecoffs,log2(as.numeric(ecoffs[ecoffs$DRUG==drug,'ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG==drug,'COEF'][1])
}

png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig2/promoter_v_body.png', width = 1028, height =768)
ggplot(sing_res[sing_res$DRUG %in% c('AMI','KAN','ETH','INH','EMB'),],aes(x=factor(DRUG,levels=c('INH','ETH','AMI','KAN','EMB')),y=COEF,weight=COUNTS, size=COUNTS, fill=factor(PROM,c(0,1)))) +
  geom_violin() + theme_classic() +
  scale_fill_manual(values = c( 'grey','#de7569')) +
  geom_point(position=position_jitterdodge()) +
  scale_x_discrete() +
  annotate("segment", y = prom_ecoffs[1], yend = prom_ecoffs[1], x = 0.6, xend = 1.4,color='tan', size=1) +
  annotate("segment", y = prom_ecoffs[2], yend = prom_ecoffs[2], x = 1.6, xend = 2.4,color='tan', size=1) +
  annotate("segment", y = prom_ecoffs[3], yend = prom_ecoffs[3], x = 2.6, xend = 3.4,color='tan', size=1) +
  annotate("segment", y = prom_ecoffs[4], yend = prom_ecoffs[4], x = 3.6, xend = 4.4,color='tan', size=1) +
  annotate("segment", y = prom_ecoffs[5], yend = prom_ecoffs[5], x = 4.6, xend = Inf,color='tan', size=1) +
  theme(text=element_text(size=24, family="Arial Rounded MT Bold", color='black'), axis.title.x=element_blank(), axis.title.y=element_blank())
dev.off()

#create wide plot for Figure 2A
sing_res$PROM <- 'SNP'
sing_res[sing_res$POSITION<0,'PROM'] = 'PROM'
sing_res$NON <- 0
sing_res[grep('!',sing_res$MUTS),'PROM'] = 'NONSENSE'
sing_res$INDEL <- 0
sing_res[grep('indel',sing_res$MUTS),'PROM'] = 'INDEL'
sing_res$PROM <- factor(sing_res$PROM)

png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig2/wide_effect_v_abundance.png', width = 1028, height =768)
ggplot(sing_res[sing_res$BH_PVAL<0.05,],aes(label=sing_res[sing_res$BH_PVAL<0.05,'MUTS'], x=factor(DRUG, levels=c('INH','ETH','RIF','RFB','EMB','LEV','MXF','AMI','KAN','LZD','DLM','BDQ','CFZ')),y=COEF,color=factor(PROM,c('PROM','SNP','NONSENSE','INDEL')), size=COUNTS)) +
  theme_classic() + xlab('DRUG') + ylab('Effect on log2MIC') + labs(color = "Mutation Type", size='Mutation count', shape='Homoplasy') +
  scale_x_discrete() +
  annotate("rect", xmin = c(1.5,3.5,5.5,7.5,9.5,11.5), xmax = c(2.5,4.5,6.5,8.5,10.5,12.5), ymin = -Inf, ymax = Inf, alpha = 0.3, fill = "grey") +
  geom_point(aes(shape=as.character(ifelse(sing_res[sing_res$BH_PVAL<0.05,'HOMOPLASY.COUNT']>0,'True','False'))), position=position_jitterdodge()) +
  annotate("segment",
           x = c(0.7, 1.7, 2.7, 3.7, 4.7, 5.7, 6.7, 7.7, 8.7, 9.7, 10.7, 11.7, 12.7),
           xend = c(1.3, 2.3, 3.3, 4.3, 5.3, 6.3, 7.3, 8.3, 9.3, 10.3, 11.3, 12.3, 13.3),
           y = c(ecoffs[ecoffs$DRUG=='INH','ECOFF'],ecoffs[ecoffs$DRUG=='ETH','ECOFF'],ecoffs[ecoffs$DRUG=='RIF','ECOFF'],ecoffs[ecoffs$DRUG=='RFB','ECOFF'],ecoffs[ecoffs$DRUG=='EMB','ECOFF'],ecoffs[ecoffs$DRUG=='LEV','ECOFF'],ecoffs[ecoffs$DRUG=='MXF','ECOFF'],ecoffs[ecoffs$DRUG=='AMI','ECOFF'],ecoffs[ecoffs$DRUG=='KAN','ECOFF'],ecoffs[ecoffs$DRUG=='LZD','ECOFF'],ecoffs[ecoffs$DRUG=='DLM','ECOFF'],ecoffs[ecoffs$DRUG=='BDQ','ECOFF'],ecoffs[ecoffs$DRUG=='CFZ','ECOFF']),
           yend = c(ecoffs[ecoffs$DRUG=='INH','ECOFF'],ecoffs[ecoffs$DRUG=='ETH','ECOFF'],ecoffs[ecoffs$DRUG=='RIF','ECOFF'],ecoffs[ecoffs$DRUG=='RFB','ECOFF'],ecoffs[ecoffs$DRUG=='EMB','ECOFF'],ecoffs[ecoffs$DRUG=='LEV','ECOFF'],ecoffs[ecoffs$DRUG=='MXF','ECOFF'],ecoffs[ecoffs$DRUG=='AMI','ECOFF'],ecoffs[ecoffs$DRUG=='KAN','ECOFF'],ecoffs[ecoffs$DRUG=='LZD','ECOFF'],ecoffs[ecoffs$DRUG=='DLM','ECOFF'],ecoffs[ecoffs$DRUG=='BDQ','ECOFF'],ecoffs[ecoffs$DRUG=='CFZ','ECOFF']),
           color='tan', size=1) + 
  scale_shape_manual(values = c('True' = 16, 'False' = 1))
  #geom_text_repel(size = 5, force=0.01, aes(label=ifelse(sing_res[sing_res$BH_PVAL<0.05,'COUNTS']>350,sing_res[sing_res$BH_PVAL<0.05,'MUTS'],''),hjust=0, vjust=0), color='black')
dev.off()

out$PROM <- 'SNP'
out[out$POSITION<0,'PROM'] = 'PROM'
out$NON <- 0
out[grep('!',out$MUTS),'PROM'] = 'NONSENSE'
out$INDEL <- 0
out[grep('indel',out$MUTS),'PROM'] = 'INDEL'
out$PROM <- factor(out$PROM)

png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig2/wide_effect_v_abundance.png', width = 1028, height =768)
ggplot(out[out$BH_PVAL<0.05 & out$COUNTS>=3,],aes(label=out[out$BH_PVAL<0.05 & out$COUNTS>=3,'MUTS'], x=factor(DRUG, levels=c('INH','ETH','RIF','RFB','EMB','LEV','MXF','AMI','KAN','LZD','DLM','BDQ','CFZ')),y=COEF,color=factor(PROM,c('PROM','SNP','NONSENSE','INDEL')), size=COUNTS+.1)) +
  theme_classic() + xlab('DRUG') + ylab('Effect on log2MIC') + labs(color = "Mutation Type", size='Mutation count', shape='Homoplasy') +
  scale_x_discrete() +
  geom_hline(yintercept = 0, color='black') +
  annotate("rect", xmin = c(1.5,3.5,5.5,7.5,9.5,11.5), xmax = c(2.5,4.5,6.5,8.5,10.5,12.5), ymin = -Inf, ymax = Inf, alpha = 0.3, fill = "grey") +
  geom_point(aes(shape=as.character(out[out$BH_PVAL<0.05 & out$COUNTS>=3,'HOMOPLASY2'])), position=position_jitterdodge()) +
  annotate("segment",
           x = c(0.7, 1.7, 2.7, 3.7, 4.7, 5.7, 6.7, 7.7, 8.7, 9.7, 10.7, 11.7, 12.7),
           xend = c(1.3, 2.3, 3.3, 4.3, 5.3, 6.3, 7.3, 8.3, 9.3, 10.3, 11.3, 12.3, 13.3),
           y = c(ecoffs[ecoffs$DRUG=='INH','ECOFF'],ecoffs[ecoffs$DRUG=='ETH','ECOFF'],ecoffs[ecoffs$DRUG=='RIF','ECOFF'],ecoffs[ecoffs$DRUG=='RFB','ECOFF'],ecoffs[ecoffs$DRUG=='EMB','ECOFF'],ecoffs[ecoffs$DRUG=='LEV','ECOFF'],ecoffs[ecoffs$DRUG=='MXF','ECOFF'],ecoffs[ecoffs$DRUG=='AMI','ECOFF'],ecoffs[ecoffs$DRUG=='KAN','ECOFF'],ecoffs[ecoffs$DRUG=='LZD','ECOFF'],ecoffs[ecoffs$DRUG=='DLM','ECOFF'],ecoffs[ecoffs$DRUG=='BDQ','ECOFF'],ecoffs[ecoffs$DRUG=='CFZ','ECOFF']),
           yend = c(ecoffs[ecoffs$DRUG=='INH','ECOFF'],ecoffs[ecoffs$DRUG=='ETH','ECOFF'],ecoffs[ecoffs$DRUG=='RIF','ECOFF'],ecoffs[ecoffs$DRUG=='RFB','ECOFF'],ecoffs[ecoffs$DRUG=='EMB','ECOFF'],ecoffs[ecoffs$DRUG=='LEV','ECOFF'],ecoffs[ecoffs$DRUG=='MXF','ECOFF'],ecoffs[ecoffs$DRUG=='AMI','ECOFF'],ecoffs[ecoffs$DRUG=='KAN','ECOFF'],ecoffs[ecoffs$DRUG=='LZD','ECOFF'],ecoffs[ecoffs$DRUG=='DLM','ECOFF'],ecoffs[ecoffs$DRUG=='BDQ','ECOFF'],ecoffs[ecoffs$DRUG=='CFZ','ECOFF']),
           color='tan', size=1) + 
  scale_shape_manual(values = c('True' = 16, 'False' = 1)) +
  scale_size_binned(trans='log10',range=c(0,6),breaks=c(4,10,50,300),labels=c("4","10","50","300"),guide="legend") +
  theme(text=element_text(size=18, family="Arial Rounded MT Bold", color='black'), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position=c(.95,.95), legend.title=element_text(size=12), legend.text=element_text(size=12), legend.background=element_blank())
dev.off()

#Figure 3
#non_int <- sig[-grep(':',sig$MUTS),]
non_int <- sig
non_int$RDRR <- FALSE
non_int[non_int$DRUG=='RIF' & non_int$GENE=='rpoB' & (non_int$POSITION>424 & non_int$POSITION<=452),'RDRR'] = 'RIF'
non_int[intersect(grep('indel',non_int$MUTS),grep('rpoB',non_int$MUTS)),'RDRR'] = 'RIF'
non_int[non_int$DRUG=='RIF' & non_int$GENE=='rpoC','RDRR'] = 'RIF'
non_int[non_int$DRUG=='RIF' & non_int$GENE=='Rv2752c','RDRR'] = 'RIF'

#test for difference in MICs between RDRR and non-RDRR mutations
wilcox.test(non_int[non_int$DRUG=='RIF' & non_int$COUNTS>2 & non_int$GENE=='rpoB' & non_int$RDRR=='RIF' & non_int$COEF>0,'COEF'],non_int[non_int$DRUG=='RIF'  & non_int$COUNTS>2 & non_int$GENE=='rpoB' & non_int$RDRR!='RIF' & non_int$COEF>0,'COEF'])


#RIF ONLY
png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig3/rif_genes_w_RDRR.png', width = 1028, height =768)
ggplot(non_int[(non_int$DRUG %in% c('RIF')) & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')),],aes(x=factor(GENE),y=COEF,fill=RDRR, size=COUNTS, weight=COUNTS)) +
  geom_violin() + theme_classic() +
  geom_point(position=position_jitterdodge(jitter.width = .15)) +
  geom_hline(yintercept = log2(as.numeric(ecoffs[ecoffs$DRUG=='RIF','ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG=='RIF','COEF'][1], color='red') +
  geom_hline(yintercept = 0) + xlab('Gene') + ylab('Effect on log2MIC')
dev.off()

#RIF + RFB
#unpaired wilcox
wilcox.test(non_int[non_int$DRUG=='RIF' & non_int$GENE=='rpoB','COEF'],non_int[non_int$DRUG=='RFB' & non_int$GENE=='rpoB','COEF'])
#paired wilcox
wilcox.test(non_int[non_int$DRUG=='RIF' & non_int$GENE=='rpoB' & (non_int$MUTS %in% non_int[non_int$DRUG=='RFB' & non_int$GENE=='rpoB','MUTS']) & order(non_int$MUTS),'COEF'], non_int[non_int$DRUG=='RFB' & non_int$GENE=='rpoB' & (non_int$MUTS %in% non_int[non_int$DRUG=='RIF' & non_int$GENE=='rpoB','MUTS']) & order(non_int$MUTS),'COEF'], paired=TRUE)
tmp = non_int[non_int$DRUG=='RIF' & non_int$GENE=='rpoB' & (non_int$MUTS %in% non_int[non_int$DRUG=='RFB' & non_int$GENE=='rpoB','MUTS']) & order(non_int$MUTS),'COEF']
ggplot(non_int$MUTS %in% tmp$MUTS, aes(x=)) +
  geom_point()


non_int[non_int$DRUG=='RFB','RDRR'] = 'FALSEF'
non_int[non_int$DRUG=='RFB' & (non_int$POSITION>424 & non_int$POSITION<=452),'RDRR'] = 'RIFB'
non_int[non_int$MUTS=='rpoB_1288_indel' & non_int$DRUG=='RFB','RDRR'] = 'RIFB'
non_int[non_int$DRUG=='RFB' & non_int$GENE=='rpoC','RDRR'] = 'RFB'
non_int[non_int$DRUG=='RFB' & non_int$GENE=='Rv2752c','RDRR'] = 'RFB'

png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig3/rif_rfb_genes_w_RDRR.png', width = 1028, height =768)
ggplot(non_int[(non_int$DRUG %in% c('RFB','RIF')) & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID','rpoA','rpoC')) & (non_int$BH_PVAL<0.05) & (!non_int$MUTS %in% c('rpoB_P45S','rpoC_F452L','rpoC_V483G','rpoC_D485Y','rpoC_I491T','rpoC_G519D','rpoC_G332S','rpoC_T812I','rpoC_G519D','rpoC_V483A','rpoC_I491V','rpoC_N698S','rpoC_P1040R','rpoC_G433S')),],aes(x=factor(GENE),y=COEF,fill=RDRR, size=COUNTS, weight=COUNTS)) +
  geom_violin(position=position_dodge(width = 1)) + theme_classic() + labs(shape='Homoplasy') +
  scale_fill_manual(values=c("#b3ffff", "#cccccc","#cccccc", "#b3ffff","#cccccc", "#b3ffff")) +
  geom_point(aes(group=RDRR, shape=as.character(non_int[(non_int$DRUG %in% c('RFB','RIF')) & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID','rpoA','rpoC')) & (non_int$BH_PVAL<0.05) & (!non_int$MUTS %in% c('rpoB_P45S','rpoC_F452L','rpoC_V483G','rpoC_D485Y','rpoC_I491T','rpoC_G519D','rpoC_G332S','rpoC_T812I','rpoC_G519D','rpoC_V483A','rpoC_I491V','rpoC_N698S','rpoC_P1040R','rpoC_G433S')),'HOMOPLASY2'])), position=position_jitterdodge(jitter.width = 0.15, dodge.width = 1)) +
  scale_shape_manual(values = c('True' = 16, 'False' = 1)) +
  scale_size_continuous(range = c(4,16)) +
  geom_hline(yintercept = log2(as.numeric(ecoffs[ecoffs$DRUG=='RIF','ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG=='RIF','COEF'][1], color='#72d2e8', size=1) +
  geom_hline(yintercept = log2(as.numeric(ecoffs[ecoffs$DRUG=='RFB','ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG=='RFB','COEF'][1], color="#333333", size=1) +
  geom_hline(yintercept = 0) + xlab('Gene') + ylab('Effect on log2MIC')
dev.off()

#create ordered drug comparison plots
drug1 <- 'LEV'
drug2 <- 'MXF'
tmp <- non_int[(non_int$DRUG %in% c(drug1,drug2)) & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')),]
tmp <- tmp[order(-tmp$COEF),]
counts_5_lab <- tmp[tmp$COUNTS>10 | tmp$HOMOPLASY==TRUE,"COUNTS"]
tmp$COUNTS2 <- tmp$COUNTS
tmp$COUNTS <- jitter(tmp$COUNTS)
tmp <- tmp[tmp$MUTS %in% tmp$MUTS[duplicated(tmp$MUTS)],]
tmp$COEF1 <- tmp[tmp$DRUG==drug1,'COEF']
tmp$COEF2 <- tmp[tmp$DRUG==drug2,'COEF']



ggplot(tmp[tmp$COUNTS2>10 | tmp$HOMOPLASY2==TRUE,], aes(x=factor(tmp[(tmp$COUNTS2>10 | tmp$HOMOPLASY2==TRUE) & order(-tmp$COEF),'MUTS'],levels=unique(tmp[(tmp$COUNTS2>10 | tmp$HOMOPLASY2==TRUE) & order(-tmp$COEF),'MUTS'])), y=COEF, color=DRUG, shape=HOMOPLASY2)) + 
  #scale_x_discrete(labels = counts_5_lab, limits=levels(tmp[tmp$COUNTS2>10 | tmp$HOMOPLASY==TRUE,'COUNTS'])) +
  geom_hline(yintercept=round(min(tmp[tmp$COUNTS2>10 | tmp$HOMOPLASY2==TRUE,'ci_lo'])):round(max(tmp[tmp$COUNTS2>10 | tmp$HOMOPLASY2==TRUE,'ci_hi'])), color='lightgrey') +
  geom_hline(yintercept=log2(as.numeric(ecoffs[ecoffs$DRUG==drug1,'ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG==drug1,'COEF'][1], color='#b3ffff', alpha=1, linetype='dashed') +
  geom_hline(yintercept=log2(as.numeric(ecoffs[ecoffs$DRUG==drug2,'ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG==drug2,'COEF'][1], color='#cccccc', alpha=1, linetype='dashed') +
  #geom_text(aes(7,log2(as.numeric(ecoffs[ecoffs$DRUG==drug,'ECOFF_CONC']))-log2(baseline),label = 'CRyPTIC ECOFF', vjust = -.5)) +
  geom_hline(yintercept=0, color='darkgrey') +
  #geom_vline(xintercept=lines, colour='green', linetype="dashed") +
  geom_pointrange(aes(ymin=ci_lo, ymax=ci_hi, group=DRUG, fatten=7, shape=as.character(tmp[tmp$COUNTS2>10 | tmp$HOMOPLASY2==TRUE,'HOMOPLASY2'])), position=position_jitterdodge()) +
  scale_shape_manual(values = c('True' = 16, 'False' = 1)) +
  scale_color_manual(values = c("#cccccc", "#b3ffff")) +
  #geom_text(size = sqrt(log(tmp[tmp$COUNTS2>10 | tmp$HOMOPLASY==TRUE,'COUNTS']))*2, aes(label=MUTS,hjust=0, vjust=-0.5,angle=90)) + #ifelse(output$ci_lo>output[output$muts=='INT','ci_hi'],as.character(output$muts),ifelse(output$ci_hi<output[output$muts=='INT','ci_lo'],as.character(output$muts),''))
  labs(x='Mutation', y='Main Effect with 95% CI') +
  theme_bw(base_size=12) +
  theme(legend.position = 'none', panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line = element_line(colour = "black"), panel.spacing.x = unit(.05, "lines"), panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#Figure 4
png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig4/inh_genes_w_RDRR.png', width = 1028, height =768)
ggplot(non_int[(non_int$DRUG %in% c('INH')) & non_int$GENE!='ahpC' & non_int$COUNTS>2 & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')) & non_int$BH_PVAL<0.05,],aes(x=factor(GENE),y=COEF, size=COUNTS, weight=COUNTS)) +
  theme_classic() + labs(shape='Homoplasy') +
  scale_x_discrete() +
  annotate("rect", xmin = c(1.5,3.5,5.5), xmax = c(2.5,4.5,6.5), ymin = -Inf, ymax = Inf, alpha = 0.3, fill = "grey") +
  scale_fill_brewer(palette="Dark2") +
  geom_point(aes(group=GENE, shape=as.character(non_int[(non_int$DRUG %in% c('INH')) & non_int$GENE!='ahpC' & non_int$COUNTS>2 & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')) & non_int$BH_PVAL<0.05,'HOMOPLASY2'])), position=position_jitterdodge(jitter.width = .25)) +
  scale_shape_manual(values = c('True' = 16, 'False' = 1)) +
  scale_size_continuous(range = c(4,16)) +
  geom_hline(yintercept = log2(as.numeric(ecoffs[ecoffs$DRUG=='INH','ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG=='INH','COEF'][1], color='red', size=1) +
  geom_hline(yintercept = 0) + xlab('Gene') + ylab('Effect on log2MIC') + theme(panel.grid.major = element_line(colour="lightgrey", size=0.5), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(text=element_text(size=18, family="Arial Rounded MT Bold", color='black'), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank())
dev.off()

png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig4/emb_genes_w_RDRR.png', width = 1028, height =768)
ggplot(non_int[(non_int$DRUG %in% c('EMB')) & non_int$COUNTS>2 & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')) & non_int$BH_PVAL<0.05,],aes(x=factor(GENE),y=COEF, size=COUNTS, weight=COUNTS)) +
  theme_classic() + labs(shape='Homoplasy') +
  scale_x_discrete() +
  annotate("rect", xmin = c(1.5,3.5), xmax = c(2.5,4.5), ymin = -Inf, ymax = Inf, alpha = 0.3, fill = "grey") +
  scale_fill_brewer(palette="Dark2") +
  geom_point(aes(group=GENE, shape=as.character(non_int[(non_int$DRUG %in% c('EMB')) & non_int$COUNTS>2 & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')) & non_int$BH_PVAL<0.05,'HOMOPLASY2'])), position=position_jitterdodge(jitter.width = .25)) +
  scale_shape_manual(values = c('True' = 16, 'False' = 1)) +
  scale_size_continuous(range = c(4,16)) +
  geom_hline(yintercept = log2(as.numeric(ecoffs[ecoffs$DRUG=='EMB','ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG=='EMB','COEF'][1], color='red', size=1) +
  geom_hline(yintercept = 0) + xlab('Gene') + ylab('Effect on log2MIC') + theme(panel.grid.major = element_line(colour="lightgrey", size=0.5), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(text=element_text(size=18, family="Arial Rounded MT Bold", color='black'), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank())
dev.off()

#print mutations by gene for INH
non_int[non_int$DRUG=='INH' & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>%
  group_by(GENE) %>%
  summarise(no_rows = length(MUTS))
#print mutations below ECOFF by gene for INH
non_int[non_int$DRUG=='INH' & non_int$ci_lo<non_int$ECOFF & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>%
  group_by(GENE) %>%
  summarise(no_rows = length(MUTS))

#print mutations by gene for EMB
non_int[non_int$DRUG=='EMB' & non_int$COEF>0 & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>%
  group_by(GENE) %>%
  summarise(no_rows = length(MUTS))
#print mutations below ECOFF by gene for EMB
non_int[non_int$DRUG=='EMB' & non_int$ci_lo<non_int$ECOFF & non_int$COEF>0 & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>%
  group_by(GENE) %>%
  summarise(no_rows = length(MUTS))


#Figure 5
png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig5/FQ_genes.png', width = 1028, height =768)
ggplot(non_int[(non_int$DRUG %in% c('LEV','MXF')) & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')) & non_int$BH_PVAL<0.05,],aes(x=factor(GENE),y=COEF,fill=DRUG, size=COUNTS, weight=COUNTS)) +
  geom_violin() + theme_classic() + labs(shape='Homoplasy') +
  scale_fill_manual(values=c("pink", "#66b9be")) +
  geom_point(aes(group=DRUG, shape=as.character(non_int[(non_int$DRUG %in% c('LEV','MXF')) & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')) & non_int$BH_PVAL<0.05,'HOMOPLASY2'])), position=position_jitterdodge(jitter.width = .15)) +
  scale_shape_manual(values = c('True' = 16, 'False' = 1)) +
  scale_size_continuous(range = c(4,16)) +
  geom_hline(yintercept = log2(as.numeric(ecoffs[ecoffs$DRUG=='LEV','ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG=='LEV','COEF'][1], color='pink', size=1) +
  geom_hline(yintercept = log2(as.numeric(ecoffs[ecoffs$DRUG=='MXF','ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG=='MXF','COEF'][1], color='#66b9be', size=1) +
  geom_hline(yintercept = 0) + xlab('Gene') + ylab('Effect on log2MIC')
dev.off()

#print mutations by gene for FQ
non_int[non_int$DRUG %in% c('LEV','MXF') & non_int$COEF>0 & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>%
  group_by(GENE,DRUG) %>%
  summarise(no_rows = length(MUTS))
#print mutations below ECOFF by gene for FQ
non_int[non_int$DRUG %in% c('LEV','MXF') & non_int$ci_lo<non_int$ECOFF & non_int$COEF>0 & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>%
  group_by(GENE,DRUG) %>%
  summarise(no_rows = length(MUTS))

png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig5/AG_genes.png', width = 1028, height =768)
ggplot(non_int[(non_int$DRUG %in% c('AMI','KAN')) & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')) & non_int$BH_PVAL<0.05,],aes(x=factor(GENE),y=COEF,color=DRUG, size=COUNTS, weight=COUNTS)) +
  theme_classic() + #geom_violin() +
  scale_x_discrete() + labs(shape='Homoplasy') +
  annotate("rect", xmin = c(1.5,3.5,5.5), xmax = c(2.5,4.5,6.5), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "grey") +
  geom_point(aes(group=DRUG, shape=as.character(non_int[(non_int$DRUG %in% c('AMI','KAN')) & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')) & non_int$BH_PVAL<0.05,'HOMOPLASY2'])), position=position_jitterdodge(jitter.width = .15)) +
  scale_shape_manual(values = c('True' = 19, 'False' = 1)) +
  scale_size_continuous(range = c(4,16)) +
  geom_hline(yintercept = log2(as.numeric(ecoffs[ecoffs$DRUG=='AMI','ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG=='AMI','COEF'][1], color='#de7569', size=1) +
  geom_hline(yintercept = log2(as.numeric(ecoffs[ecoffs$DRUG=='KAN','ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG=='KAN','COEF'][1], color='#66b9be', size=1) +
  geom_hline(yintercept = 0) + xlab('Gene') + ylab('Effect on log2MIC')+ theme(panel.grid.major = element_line(colour="lightgrey", size=0.5), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
dev.off()

#print mutations by gene for AG
non_int[non_int$DRUG %in% c('KAN','AMI') & non_int$COEF>0 & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>%
  group_by(GENE,DRUG) %>%
  summarise(no_rows = length(MUTS))
#print mutations below ECOFF by gene for AG
non_int[non_int$DRUG %in% c('KAN','AMI') & non_int$ci_lo<non_int$ECOFF & non_int$COEF>0 & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>%
  group_by(GENE,DRUG) %>%
  summarise(no_rows = length(MUTS))


#Supplemental figs
png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig6/BDQ_genes.png', width = 1028, height =768)
ggplot(non_int[(non_int$DRUG %in% c('BDQ','CFZ')) & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')) & non_int$BH_PVAL<0.05,],aes(x=factor(GENE),y=COEF,color=DRUG, size=COUNTS, weight=COUNTS)) +
  theme_classic() + labs(shape='Homoplasy') +
  scale_x_discrete() +
  annotate("rect", xmin = c(1.5,3.5,5.5), xmax = c(2.5,4.5,6.5), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "grey") +
  geom_point(aes(group=DRUG, shape=as.character(non_int[(non_int$DRUG %in% c('BDQ','CFZ')) & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')) & non_int$BH_PVAL<0.05,'HOMOPLASY2'])), position=position_jitterdodge(jitter.width = .35)) +
  scale_shape_manual(values = c('True' = 16, 'False' = 1)) +
  geom_hline(yintercept = log2(as.numeric(ecoffs[ecoffs$DRUG=='BDQ','ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG=='BDQ','COEF'][1], color='#de7569') +
  geom_hline(yintercept = log2(as.numeric(ecoffs[ecoffs$DRUG=='CFZ','ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG=='CFZ','COEF'][1], color='#66b9be') +
  geom_hline(yintercept = 0) + xlab('Gene') + ylab('Effect on log2MIC') +
  geom_text_repel(size = 3, force = 5, family="Arial Rounded MT Bold", aes(label=ifelse(COUNTS>20,MUTS,''),hjust=0, vjust=0), color='black') +
  theme(text=element_text(size=18, family="Arial Rounded MT Bold", color='black'), axis.title.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

#print mutations by gene for BDQ/CFZ
non_int[non_int$DRUG %in% c('BDQ','CFZ') & non_int$COEF>0 & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>%
  group_by(GENE,DRUG) %>%
  summarise(no_rows = length(MUTS))

#print mutations by gene for BDQ/CFZ
non_int[non_int$DRUG %in% c('BDQ','CFZ') & non_int$COEF>0 & non_int$ci_lo<non_int$ECOFF & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>%
  group_by(GENE,DRUG) %>%
  summarise(no_rows = length(MUTS))

png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig6/DLM_LZD_genes.png', width = 1028, height =768)
ggplot(non_int[(non_int$DRUG %in% c('LZD','DLM')) & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')) & non_int$BH_PVAL<0.05,],aes(x=factor(GENE),y=COEF,color=DRUG, size=COUNTS, weight=COUNTS)) +
  theme_classic() + labs(shape='Homoplasy') +
  scale_x_discrete() +
  annotate("rect", xmin = c(1.5,3.5,5.5), xmax = c(2.5,4.5,6.5), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "grey") +
  geom_point(aes(group=DRUG, shape=as.character(non_int[(non_int$DRUG %in% c('DLM','LZD')) & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')) & non_int$BH_PVAL<0.05,'HOMOPLASY2'])), position=position_jitterdodge(jitter.width = .35)) +
  scale_shape_manual(values = c('True' = 16, 'False' = 1)) +
  geom_hline(yintercept = log2(as.numeric(ecoffs[ecoffs$DRUG=='DLM','ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG=='DLM','COEF'][1], color='#de7569') +
  geom_hline(yintercept = log2(as.numeric(ecoffs[ecoffs$DRUG=='LZD','ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG=='LZD','COEF'][1], color='#66b9be') +
  geom_hline(yintercept = 0) + xlab('Gene') + ylab('Effect on log2MIC') +
  geom_text_repel(size = 3, force = 5, family="Arial Rounded MT Bold", aes(label=ifelse(COUNTS>20,MUTS,''),hjust=0, vjust=0), color='black') +
  theme(text=element_text(size=18, family="Arial Rounded MT Bold", color='black'), axis.title.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

#print mutations by gene for DLM
non_int[non_int$DRUG=='DLM' & non_int$COEF>0 & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>%
  group_by(GENE) %>%
  summarise(no_rows = length(MUTS))

#print mutations by gene for DLM
non_int[non_int$DRUG=='DLM' & non_int$COEF>0 & non_int$ci_lo<non_int$ECOFF & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>%
  group_by(GENE) %>%
  summarise(no_rows = length(MUTS))

png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig6/ETH_genes.png', width = 1028, height =768)
ggplot(non_int[(non_int$DRUG %in% c('ETH')) & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')) & non_int$BH_PVAL<0.05,],aes(x=factor(GENE),y=COEF,color=GENE, size=COUNTS, weight=COUNTS)) +
  theme_classic() + labs(shape='Homoplasy') +
  scale_x_discrete() +
  annotate("rect", xmin = c(1.5,3.5,5.5), xmax = c(2.5,4.5,6.5), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "grey") +
  geom_point(aes(group=DRUG, shape=as.character(non_int[(non_int$DRUG %in% c('ETH')) & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')) & non_int$BH_PVAL<0.05,'HOMOPLASY2'])), position=position_jitterdodge(jitter.width = .35, seed = 1)) +
  scale_shape_manual(values = c('True' = 16, 'False' = 1)) +
  geom_hline(yintercept = log2(as.numeric(ecoffs[ecoffs$DRUG=='ETH','ECOFF_CONC']))-out[out$MUTS=='_cons' & out$DRUG=='ETH','COEF'][1], color='#de7569') +
  geom_hline(yintercept = 0) + xlab('Gene') + ylab('Effect on log2MIC') +
  #geom_text_repel(size = 3, force=3, position=position_jitter(seed=1,height = 0), family="Arial Rounded MT Bold", aes(label=ifelse(COUNTS>50,MUTS,''),hjust=0, vjust=0), color='black') +
  theme(text=element_text(size=18, family="Arial Rounded MT Bold", color='black'), axis.title.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

#print mutations by gene for ETH
non_int[non_int$DRUG=='ETH' & non_int$COEF>0 & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>%
  group_by(GENE) %>%
  summarise(no_rows = length(MUTS))

#print mutations by gene for ETH
non_int[non_int$DRUG=='ETH' & non_int$COEF>0 & non_int$ci_lo<non_int$ECOFF & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')),] %>%
  group_by(GENE) %>%
  summarise(no_rows = length(MUTS))

#Sensitive mutations plot
sens <- sig[sig$COEF<0,]
sens <- sens[sens$GENE!='DRUG' & sens$GENE!='SITEID' & sens$GENE!='LINEAGE',]
sing_sens <- sens
#sing_sens <- sens[-grep(':',sens$MUTS),]
sing_sens$PROM <- 'SNP'
sing_sens[sing_sens$POSITION<0,'PROM'] = 'PROM'
sing_sens$NON <- 0
sing_sens[grep('!',sing_sens$MUTS),'PROM'] = 'NONSENSE'
sing_sens$INDEL <- 0
sing_sens[grep('indel',sing_sens$MUTS),'PROM'] = 'INDEL'
sing_sens$PROM <- factor(sing_sens$PROM)

sing_sens[sing_sens$BH_PVAL<0.05 & sing_sens$COUNTS>3,] %>%
  group_by(DRUG) %>%
  summarise(no_rows = length(DRUG))


png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig6/sens_wide.png', width = 1028, height =768)
ggplot(sing_sens[sing_sens$COUNTS>=3,],aes(label=MUTS, x=factor(DRUG, levels=c('INH','ETH','RIF','RFB','EMB','LEV','MXF','AMI','KAN','LZD','DLM','BDQ','CFZ')),y=COEF,color=factor(PROM,c('PROM','SNP','NONSENSE','INDEL')), size=COUNTS)) +
  theme_classic() + xlab('DRUG') + ylab('Effect on log2MIC') + labs(color = "Mutation Type", size='Mutation count', shape='Homoplasy') +
  scale_x_discrete() +
  annotate("rect", xmin = c(1.5,3.5,5.5,7.5,9.5), xmax = c(2.5,4.5,6.5,8.5,10.5), ymin = -Inf, ymax = Inf, alpha = 0.3, fill = "grey") +
  geom_point(aes(group=DRUG, shape=as.character(sing_sens[sing_sens$COUNTS>=3,'HOMOPLASY2'])), position=position_jitterdodge()) +
  scale_shape_manual(values = c('True' = 16, 'False' = 1))
  #geom_text_repel(size = 4.5, aes(label=ifelse(HOMOPLASY=='True',MUTS,''),hjust=0, vjust=0), color='black') +
  #geom_text_repel(size = 5, aes(label=ifelse(COUNTS>50,MUTS,''),hjust=0, vjust=0), color='black')
dev.off()

#read in interaction outputs to get final results
flag=0
for (drug in drug_list) { 
  
  #read in input data file prepared in Jupyter to get mutation names
  dat <- read.csv(file = paste('/Users/joshuacarter/Desktop/cryptic-analysis/Interval_regression/Output_files/int/R/',drug,'/',drug,'_main_effects_homoplasy.csv',sep=''), header = TRUE)
  dat$DRUG = drug
  if (flag==0) {
    out <- dat
  } else {out <- rbind(out,dat)}
  flag <- 1
  
}
sig <- out[out$BH_PVAL<0.05,]

#create interaction plots
pheno_out <- c()
for (drug in drug_list) { 
  #read input files to get phenotypes
  pheno <- read.csv(paste('/Users/joshuacarter/Desktop/cryptic-analysis/Interval_regression/Input_files/',drug,'_interaction.csv',sep=''), header = TRUE)
  pheno_out <- append(pheno_out,max(pheno$UPLOG2MIC))
}
pheno_out <- data.frame(drug_list,pheno_out)
pheno_out$DRUG <- pheno_out$drug_list

int <- sig[sig$BH_PVAL<0.05,]
int <- int[grep(':',int$MUTS),]
int <- int[int$GENE!='DRUG' & int$GENE!='SITEID' & int$GENE!='LINEAGE',]

temp <- regmatches(int$MUTS, regexpr(":", int$MUTS), invert = TRUE)
mat  <- matrix(unlist(temp), ncol=2, byrow=TRUE)
int$MUT1 <- mat[,1]
int$MUT2 <- mat[,2]
tmp <- sig
tmp$MUT1 <- tmp$MUTS
tmp$MUT2 <- tmp$MUTS
tmp$COEF1 <- tmp$COEF
tmp$COEF2 <- tmp$COEF
int <- merge(int, tmp %>% select('COEF1','DRUG','MUT1'),by=c('DRUG','MUT1'))
int <- merge(int, tmp %>% select('COEF2','DRUG','MUT2'),by=c('DRUG','MUT2'))
int$baseline <- 0
for (drug in drug_list) {
  int[int$DRUG==drug,'baseline'] <- out[out$MUTS=='_cons' & out$DRUG==drug,'COEF'][1]
}
int$ADD_COEF <- int$COEF1 + int$COEF2
int$COMB_COEF <- int$COEF1 + int$COEF2 + int$COEF
int <- merge(int,pheno_out %>% select('DRUG','pheno_out'),by='DRUG')
int$shape <- ifelse(int$COEF1>0,ifelse(int$COEF2>0,'16','1'),'1')


#plot negative interactions to check if just adjusting for censoring
ggplot(int[int$COMB_COEF>int$pheno_out,],aes(label=MUTS,x=pheno_out-COMB_COEF,y=COEF,color=DRUG, size=COUNTS)) +
  theme_classic() + xlab('Difference between pure additivity and highest MIC measured') + ylab('Inverse of effect on log2MIC') + labs(color = "Mutation Type", size='Mutation count') +
  geom_point() + geom_pointrange(aes(ymin=abs(ci_lo), ymax=abs(ci_hi), size=1)) +
  geom_abline() + geom_abline(intercept = 1, color='grey') + geom_abline(intercept = 2.7, linetype='dashed') +
  geom_text_repel(size = 3, aes(label=ifelse(COUNTS>0,MUTS,''),hjust=0, vjust=0), color='black')

#Interaction plot
png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig6/int_noncensored.png', width = 1028, height =768)
ggplot(int[abs(int$COMB_COEF)<=int$pheno_out,],aes(label=MUTS, x=factor(DRUG, levels=c('INH','ETH','RIF','RFB','EMB','LEV','MXF','AMI','KAN','LZD','DLM','BDQ','CFZ')),y=COEF,color=DRUG, shape=shape)) +
  theme_classic() + xlab('DRUG') + ylab('Effect on log2MIC') + labs(color = "Mutation Type", size='Mutation count', shape='Res Status') +
  scale_shape_manual(values = c(1,16)) +
  geom_pointrange(position=position_jitterdodge(seed = 1), aes(ymin=ci_lo, ymax=ci_hi)) +
  geom_hline(yintercept=0,color='black') +
  geom_text_repel(size = 3, aes(label=ifelse(COUNTS>0,MUTS,''),hjust=0, vjust=0), color='black')
dev.off()

png('/Users/joshuacarter/Desktop/paper-cryptic-mic/Figures/Fig6/int_noncensored.png', width = 1028, height =768)
tmp <- bind_rows(int[int$GENE %in% c('rpoB','embA','embB') & ((int$MUT1 %in% c('embA_g-43c','embA_-44_indel','embB_M306I')) | (int$MUT2=='rpoB_D435G')),],sig[sig$DRUG %in% c('EMB','RIF') & sig$MUTS %in% c('embA_g-43c','embA_-44_indel','embB_M306V','rpoB_D435G','rpoB_L430P','embB_Q497R','embB_M306I','rpoB_S450L'),])
tmp[!complete.cases(tmp$COMB_COEF),'COMB_COEF'] <- tmp[!complete.cases(tmp$COMB_COEF),'COEF']
ggplot(tmp, aes(x=factor(DRUG), y=COMB_COEF, color=DRUG, size=COUNTS)) +
  theme_classic() + xlab('DRUG') + ylab('Effect on log2MIC') + labs(color = "Mutation Type", size='Mutation count') +
  geom_point(position=position_jitter(seed=1,height = 0)) +
  ylim(0,NA) +
  geom_text_repel(position=position_jitter(seed=1,height = 0), size=4, family="Arial Rounded MT Bold", color='black', aes(label=ifelse(COUNTS>0,MUTS,''),hjust=0, vjust=0)) +
  geom_point(aes(x=factor(DRUG),y=ADD_COEF),color='black') +
  theme(text=element_text(size=18, family="Arial Rounded MT Bold", color='black'), legend.position = 'none', axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

#Lineage plot
lineage <- sig[grep('lineage',sig$MUTS),]
ggplot(lineage[lineage$BH_PVAL<0.05,],aes(x=as.character(factor(DRUG, levels=c('INH','ETH','RIF','RFB','EMB','LEV','MXF','AMI','KAN','LZD','DLM','BDQ','CFZ'))), color=MUTS, y=COEF)) +
  theme_classic() + xlab('DRUG') + ylab('Effect on log2MIC') + labs(size='Mutation count', color='Lineage') +
  scale_x_discrete() +
  annotate("rect", xmin = c(0,2.5,4.5,6.5,8.5,10.5), xmax = c(1.5,3.5,5.5,7.5,9.5,11.5), ymin = -Inf, ymax = Inf, alpha = 0.3, fill = "grey") +
  geom_pointrange(position=position_jitterdodge(seed = 1), aes(ymin=ci_lo, ymax=ci_hi)) + geom_point(position=position_jitterdodge(seed = 1)) +
  geom_hline(yintercept=0,color='black') +
  scale_y_continuous(limits=c(-1.8,1), breaks = round(seq(-1.8,1,by=0.25),1))
  #scale_y_continuous(limits=c(-1.5,1), breaks = round(seq(min(lineage[lineage$BH_PVAL<0.05,'COEF']), max(lineage[lineage$BH_PVAL<0.05,'COEF']), by = 0.5),1))

lineage['COMP_COEF'] <- 0
for (drug in drug_list) {
  lineage[lineage$DRUG==drug,'COMP_COEF'] = lineage[lineage$DRUG==drug,'COEF'] / median(non_int[non_int$DRUG==drug & (!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')) & non_int$COEF>0,'COEF'])
}

non_int[(!non_int$GENE %in% c('DRUG','LINEAGE','SITEID')) & non_int$COEF>0,] %>%
  group_by(DRUG) %>%
  as.data.frame(summarise(no_rows = median(COEF)))


#SITE plot
siteid <- sig[c(grep('siteid',sig$MUTS), grep('lineage',sig$MUTS)),]
ggplot(siteid[siteid$BH_PVAL<0.05,],aes(x=as.character(factor(DRUG, levels=c('INH','ETH','RIF','RFB','EMB','LEV','MXF','AMI','KAN','LZD','DLM','BDQ','CFZ'))), color=MUTS, y=COEF)) +
  theme_classic() + xlab('DRUG') + ylab('Effect on log2MIC') + labs(color='Variable name', shape='Site (Technical) effect') +
  scale_x_discrete() +
  geom_hline(yintercept = 1, color='grey', linetype='dashed') +
  geom_hline(yintercept = -1, color='grey', linetype='dashed') +
  annotate("rect", xmin = c(0,2.5,4.5,6.5,8.5,10.5), xmax = c(1.5,3.5,5.5,7.5,9.5,11.5), ymin = -Inf, ymax = Inf, alpha = 0.3, fill = "grey") +
  geom_pointrange(position=position_jitterdodge(seed = 1), alpha = ifelse(siteid[siteid$BH_PVAL<0.05,'GENE']=='SITEID',ifelse(abs(siteid$COEF)>=1,1,0.2),1), aes(ymin=ci_lo, ymax=ci_hi, shape =  ifelse(siteid[siteid$BH_PVAL<0.05,'GENE']=='SITEID','True','False'))) +
  #geom_point(aes(shape=as.character(ifelse(siteid[siteid$BH_PVAL<0.05,'GENE']=='SITEID','True','False')), alpha = ifelse(siteid[siteid$BH_PVAL<0.05,'GENE']=='SITEID',0,0.5)), position=position_jitterdodge(seed = 1)) +
  geom_hline(yintercept=0,color='black') +
  scale_shape_manual(values = c('True' = 1, 'False' = 17)) +
  scale_alpha(range = c(0.2,1)) +
  scale_y_continuous(limits=c(-2,2), breaks = round(seq(-2,2,by=0.33),1))
#scale_y_continuous(limits=c(-1.5,1), breaks = round(seq(min(lineage[lineage$BH_PVAL<0.05,'COEF']), max(lineage[lineage$BH_PVAL<0.05,'COEF']), by = 0.5),1))



#### Disputed mutations effect on other MICs

INH <- read.csv('/Users/joshuacarter/Downloads/Interval_regression/Input_files/INH_interaction.csv')
MXF <- read.csv('/Users/joshuacarter/Downloads/Interval_regression/Input_files/MXF_interaction.csv')
RIF <- read.csv('/Users/joshuacarter/Downloads/Interval_regression/Input_files/RIF_no_interaction_equal_end_short_final.csv')
