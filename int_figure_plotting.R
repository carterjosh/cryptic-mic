library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)

ecoffs <- read.csv('/Users/joshuacarter/Dropbox/paper-cryptic-mic/ecoffs.csv')

clus <- 'cluster100'
pop_control <- TRUE
no_int <- TRUE
int <- FALSE
drug_epi <- FALSE
drug_list <- c('ETH','BDQ','CFZ') #'RIF','RFB','INH','ETH','EMB','KAN','AMI','BDQ','CFZ','LEV','LZD','MXF','DLM'
#replace drug list for drugs with interacting mutants
if (int==TRUE) {int_list <- c('INH', 'EMB', 'CFZ', 'RFB', 'ETH', 'MXF', 'LEV', 'KAN', 'RIF')}

if (int==TRUE) {path <- '/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/int/'}
if (int==TRUE) {path2 <- '_interaction.csv'}
if (int==TRUE) {path3 <- paste('STATA/',clus,'/',sep='')}
if (int==TRUE) {path4 <- 'int'}


if (no_int==TRUE) {path <- '/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/no_int/'}
if (no_int==TRUE) {path2 <- '_no_interaction_equal_end.csv'}
if (no_int==TRUE) {path3 <- paste('STATA/equal_end_range/',clus,'/',sep='')}
if (no_int==TRUE) {path4 <- 'no_int'}
if (no_int==TRUE) {path5 <- '_no_interaction_equal_end_short_final.csv'}

if (drug_epi==TRUE) {path <- '/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/int/'}
if (drug_epi==TRUE) {path2 <- '_drug_epi.csv'}
if (drug_epi==TRUE) {path3 <- 'STATA/drug_epi/'}
if (drug_epi==TRUE) {path4 <- 'drug_epi'}

if (pop_control==FALSE) {path2 <- '_no_interaction_equal_end_no_pop.csv'}

for (drug in drug_list) { 
  #create directory if it does not exist
  dir.create(paste(path,'R/',drug,sep=''))
  
  #check to see if drug had interactions and if not just use base data
  if (int==TRUE && drug %in% int_list) {
    #read in input data file prepared in Jupyter to get mutation names
    dat <- read.csv(file = paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Input_files/',drug,path2,sep=''), header = TRUE)
    #read in STATA output from excel and clean labels etc. for dataframe
    out <- data.frame(read.csv(file = paste(path,path3,drug,path2,sep='')))
  }else{
    #read in input data file prepared in Jupyter to get mutation names
    dat <- read.csv(file = paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Input_files/',drug,path5,sep=''), header = TRUE)
    #read in STATA output from excel and clean labels etc. for dataframe
    if (pop_control==TRUE) {out <- data.frame(read.csv(file = paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/no_int/STATA/equal_end_range/',clus,'/',drug,'_no_interaction_equal_end.csv',sep='')))
    }else {out <- data.frame(read.csv(file = paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/no_int/STATA/equal_end_range/',clus,'/',drug,'_no_interaction_equal_end_no_pop.csv',sep='')))}
  }
  muts <- colnames(dat)
  if (drug_epi==TRUE) {
    tmp_num <- match("UNIQUEID",muts)
  }
  if (int==TRUE && drug %in% int_list) {
    tmp_num <- match("LINEAGE",muts)
  }else{tmp_num <- match("LINEAGE",muts)}
  
  if (no_int==TRUE) {
    tmp_num <- match("LINEAGE",muts)
  }
  muts <- muts[1:tmp_num-1]
  mut_counts <- colSums(dat[,1:tmp_num-1])
  mut_counts <- append(mut_counts,rep(0,length(out$var)-tmp_num+1))
  muts[grep("__",muts)] = gsub("__", ":", muts[grep("__",muts)])
  temp <- regmatches(muts, regexpr("_", muts), invert = TRUE)
  mat  <- matrix(unlist(temp), ncol=2, byrow=TRUE)
  df   <- as.data.frame(mat)
  colnames(df) <- c('GENE','POSITION')
  #if (int==TRUE && drug %in% int_list) {
  temp <- regmatches(df$POSITION[grep(":",df$POSITION)], regexpr(":", df$POSITION[grep(":",df$POSITION)]), invert = TRUE)
  mat  <- matrix(unlist(temp), ncol=2, byrow=TRUE)
  df$POSITION[grep(":",df$POSITION)] <- mat[,1]
  #}
  df$POSITION <- gsub("[^0-9.]", "", df$POSITION)
  df$POSITION <- as.numeric(gsub("^\\.", "-", df$POSITION))
  out$var <- substring(as.character(out$var), 9)
  out[2:7] <- sapply(out[2:7],as.numeric)
  if (pop_control==TRUE) {
    site_lineage <- c(rep("LINEAGE", 5),rep("SITEID", 10))
    addgene <- c(site_lineage,rep("DRUG",(length(out$var)-(tmp_num+length(site_lineage)))+1))
  }else{
    addgene <- rep("DRUG",(length(out$var)-tmp_num+1))
  }
  addpos <- c(1:(length(out$var)-tmp_num+1))
  tmp <- data.frame(addgene,addpos)
  colnames(tmp) <- c('GENE','POSITION')
  df <- rbind(df,tmp)
  muts[grep("[actg]\\.",muts)] <- gsub("\\.", "-", muts[grep("[actg]\\.",muts)])
  muts[grep("[0-9]\\.",muts)] <- gsub("\\.", "!", muts[grep("[0-9]\\.",muts)])
  muts[grep("_\\.",muts)] <- gsub("\\.", "-", muts[grep("_\\.",muts)])
  out[1:tmp_num-1,1]<- muts
  
  #create dataframe and confidence intervals
  output <- data.frame(out$var,df$GENE,as.numeric(df$POSITION),out$coef,out$stderr,out$tstat,out$pval,mut_counts)
  colnames(output) <- c('MUTS','GENE','POSITION','COEF','SE','Z','P','COUNTS')
  output$corrected_se <- output$SE*1.96  #*qnorm(1-.05/(length(muts)-1))
  output$other_se <- output$SE^2 / abs(output$COEF)
  output$ci_hi <- output$COEF+output$corrected_se
  output$ci_lo <- output$COEF-output$corrected_se
  output$BH_PVAL <- p.adjust(output$P, method = "BH")
  
  #catch stop codon mutations
  output[(abs(output$POSITION)>100) & (str_detect(output$MUTS, '-')),'MUTS'] <- gsub("-", "!", output[(abs(output$POSITION)>100) & (str_detect(output$MUTS, '-')),'MUTS'])
  output[abs(output$POSITION)>100,'POSITION'] = abs(output[abs(output$POSITION)>100,'POSITION'])
  
  
  #create main effects output dataframe, sort by gene, position, and effect and write to file
  output <- output[order(output$GENE,output$POSITION,-output$COEF),]
  output <- output[complete.cases(output),]
  if (pop_control==TRUE) {write.csv(output, file = paste(path,'R/',drug,'/',drug,'_',path4,'_',clus,'_main_effects.csv',sep=''),row.names=FALSE)
  }else{write.csv(output, file = paste(path,'R/',drug,'/',drug,'_',path4,'_',clus,'_main_effects_no_pop.csv',sep=''),row.names=FALSE)}
  baseline <- output[output$MUTS=='_cons','COEF']

  #read in homoplasy output
  if (int==TRUE) {
    output <- read.csv(paste(path,'R/',drug,'/',drug,'_main_effects_homoplasy.csv',sep=''),header=TRUE)
  }
  output <- read.csv(paste(path,'R/',drug,'/',drug,'no_int_main_effects_homoplasy.csv',sep=''),header=TRUE)
  tmp <- filter(output, grepl('indel', MUTS) & POSITION>0)
  tmp$POSITION <- tmp$POSITION / 3
  output <- output[!grepl("indel", output$MUTS),]
  output <- rbind(output,tmp)
  output <- output[order(output$GENE,output$POSITION,-output$COEF),]
  output <- output[output$MUTS!='ns[cluster100])',]
  output <- output[output$MUTS!='_cons',]
  output <- output[output$MUTS!='og2mic)',]
  output$MUTS <- gsub("[.]", "-", output$MUTS)
  output$MUTS <- gsub('-$', "!", output$MUTS)
  #removing anything with absurd standard error
  output <- output[output$corrected_se<10,]
  #remove mutations that were omitted due to collinearity
  output <- output[output$COEF!=0,]
  #keep only significant effects for plotting
  output <- output[!(output$ci_hi>0 & output$ci_lo<0),]
  output <- output[output$BH_PVAL<0.05,]
  
  #keep only epistatic effects > 1
  #output <- output[!abs(output$COEF)<0.5,]
  
  #plot mutations and save
  genes <- data.frame(output[!rev(duplicated(rev(output$GENE))),'GENE'], stringsAsFactors = FALSE)
  colnames(genes) <- c('GENE')
  genes <- filter(genes, !(GENE %in% c('Lineage','Site','DRUG')))[['GENE']]
  genes <- append(genes,c('DRUG','Lineage','Site'),after=length(genes))
  output$GENE <- factor(output$GENE, levels = genes)

  #convert to real MICs
  output$micCOEF <- 2^(output$COEF + baseline)
  output$micci_hi <- 2^(output$COEF + output$SE*1.96 + baseline)
  output$micci_lo <- 2^(output$COEF - output$SE*1.96 + baseline)
  baseline <- 2^baseline
  
  #plot ordered by abundance
  png(paste(path,'R/',drug,'/',drug,'_',path4,'_main_effects_by_abund.png',sep=''), width = 4096, height =1536)
  print(ggplot(output, aes(x=COUNTS, y=COEF)) + 
          scale_x_continuous(trans='log10') +
          geom_hline(yintercept=round(min(output$ci_lo)):round(max(output$ci_hi)), color='lightgrey') +
          #geom_hline(yintercept=log2(as.numeric(ecoffs[ecoffs$DRUG==drug,'ECOFF_CONC'])), color='darkgrey', linetype='dashed') +
          geom_hline(yintercept=0, color='darkgrey') +
          #geom_vline(xintercept=lines, colour='green', linetype="dashed") +
          geom_pointrange(aes(ymin=ci_lo, ymax=ci_hi), color=ifelse(output$ci_lo>0,'red',ifelse(output$ci_hi<0,'blue','grey')), alpha=ifelse(output$COEF==0,0,0.3)) +
          geom_text(size = 4, aes(label=MUTS,hjust=0, vjust=0,angle=90)) + #ifelse(output$ci_lo>output[output$muts=='INT','ci_hi'],as.character(output$muts),ifelse(output$ci_hi<output[output$muts=='INT','ci_lo'],as.character(output$muts),''))
          labs(x='# of occurences', y='Main Effect with 95% CI') +
          theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line = element_line(colour = "black"), panel.spacing.x = unit(.05, "lines"), panel.background = element_blank()))
  dev.off()
  
  #plot ordered by coefficient
  tmp <- output[order(-output$COEF),]
  counts_lab <- output[order(-output$COEF),"COUNTS"]
  counts_5_lab <- tmp[tmp$COUNTS>10 | tmp$HOMOPLASY!="False","COUNTS"]
  output$COUNTS2 <- output$COUNTS
  output$COUNTS <- jitter(output$COUNTS)
  
  png(paste(path,'R/',drug,'/',drug,'_',path4,'_main_effects_MIC.png',sep=''), width = 4096, height =1536)
  print(ggplot(output, aes(x=reorder(COUNTS, -micCOEF), y=micCOEF)) + 
    scale_x_discrete(labels = counts_lab) +
    scale_y_continuous(labels=comma, breaks=breaks_log(n=9, base=2), trans='log2') +
    geom_hline(yintercept=as.numeric(ecoffs[ecoffs$DRUG==drug & ecoffs$CURRENT_ECOFF==TRUE,'ECOFF_CONC']), color='tan', alpha=1, linetype='dashed') +
    geom_text(aes(7,as.numeric(ecoffs[ecoffs$DRUG==drug,'ECOFF_CONC']),label = paste('CRyPTIC ECOFF',round(as.numeric(ecoffs[ecoffs$DRUG==drug,'ECOFF_CONC']),2)), vjust = -.5)) +
    geom_hline(yintercept=baseline, color='darkgrey') +
    geom_text(aes(7,baseline,label = paste('Baseline MIC',round(baseline,2)), vjust = -.5)) +
    geom_pointrange(aes(ymin=micci_lo, ymax=micci_hi, fatten=sqrt(abs(COUNTS))), color=ifelse(output$micci_lo>baseline,'red',ifelse(output$micci_hi<baseline,'blue','grey')), alpha=ifelse(output$micCOEF==baseline,0,0.3)) +
    geom_text(size = 4, aes(label=MUTS,hjust=0, vjust=0,angle=90)) +
    labs(x='# of occurences', y='Solo MIC with 95% CI') +
    theme(panel.grid.major.y = element_line(colour = "lightgrey"),panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line = element_line(colour = "black"), panel.spacing.x = unit(.05, "lines"), panel.background = element_blank()) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
  dev.off()
  
  png(paste(path,'R/',drug,'/',drug,'_',path4,'_main_effects.png',sep=''), width = 4096, height =1536)
  print(ggplot(output, aes(x=reorder(COUNTS, -COEF), y=COEF)) + 
          scale_x_discrete(labels = counts_lab) +
          geom_hline(yintercept=round(min(output$ci_lo)):round(max(output$ci_hi)), color='lightgrey') +
          geom_hline(yintercept=log2(as.numeric(ecoffs[ecoffs$DRUG==drug,'ECOFF_CONC']))-log2(baseline), color='tan', alpha=1, linetype='dashed') +
          geom_text(aes(7,log2(as.numeric(ecoffs[ecoffs$DRUG==drug,'ECOFF_CONC']))-log2(baseline),label = 'CRyPTIC ECOFF', vjust = -.5)) +
          geom_hline(yintercept=0, color='darkgrey') +
          #geom_vline(xintercept=lines, colour='green', linetype="dashed") +
          geom_pointrange(aes(ymin=ci_lo, ymax=ci_hi, fatten=sqrt(abs(COUNTS))), color=ifelse(output$ci_lo>0,'red',ifelse(output$ci_hi<0,'blue','grey')), alpha=ifelse(output$COEF==0,0,0.3)) +
          geom_text(size = 4, aes(label=MUTS,hjust=0, vjust=0,angle=90)) + #ifelse(output$ci_lo>output[output$muts=='INT','ci_hi'],as.character(output$muts),ifelse(output$ci_hi<output[output$muts=='INT','ci_lo'],as.character(output$muts),''))
          labs(x='# of occurences', y='Main Effect with 95% CI') +
          theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line = element_line(colour = "black"), panel.spacing.x = unit(.05, "lines"), panel.background = element_blank()) +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
  dev.off()
  
  png(paste(path,'R/',drug,'/',drug,'_',path4,'_main_effects_5count.png',sep=''), width = 1024, height =768)
  print(ggplot(output[output$COUNTS2>10 | output$HOMOPLASY!="False",], aes(x=reorder(COUNTS, -COEF), y=COEF)) + 
          scale_x_discrete(labels = counts_5_lab, limits=levels(output[output$COUNTS2>10 | output$HOMOPLASY!="False",'COUNTS'])) +
          geom_hline(yintercept=round(min(output[output$COUNTS2>10 | output$HOMOPLASY!="False",'ci_lo'])):round(max(output[output$COUNTS2>10 | output$HOMOPLASY!="False",'ci_hi'])), color='lightgrey') +
          geom_hline(yintercept=log2(as.numeric(ecoffs[ecoffs$DRUG==drug,'ECOFF_CONC']))-log2(baseline), color='tan', alpha=1, linetype='dashed') +
          geom_text(aes(7,log2(as.numeric(ecoffs[ecoffs$DRUG==drug,'ECOFF_CONC']))-log2(baseline),label = 'CRyPTIC ECOFF', vjust = -.5)) +
          geom_hline(yintercept=0, color='darkgrey') +
          #geom_vline(xintercept=lines, colour='green', linetype="dashed") +
          geom_pointrange(aes(ymin=ci_lo, ymax=ci_hi, fatten=7, shape=as.character(output[output$COUNTS2>10 | output$HOMOPLASY!="False",'HOMOPLASY2'])), color=ifelse(output[output$COUNTS2>10 | output$HOMOPLASY!="False",'ci_lo']>0,'red',ifelse(output[output$COUNTS2>10 | output$HOMOPLASY!="False",'ci_hi']<0,'blue','grey')), alpha=ifelse(output[output$COUNTS2>10 | output$HOMOPLASY!="False",'COEF']==0,0,1)) +
          scale_shape_manual(values = c('True' = 16, 'False' = 1)) +
          geom_text(size = sqrt(log(output[output$COUNTS2>10 | output$HOMOPLASY!="False",'COUNTS']))*2, aes(label=MUTS,hjust=0, vjust=-0.5,angle=90)) + #ifelse(output$ci_lo>output[output$muts=='INT','ci_hi'],as.character(output$muts),ifelse(output$ci_hi<output[output$muts=='INT','ci_lo'],as.character(output$muts),''))
          labs(x='# of occurences', y='Main Effect with 95% CI') +
          theme_bw(base_size=16) +
          theme(legend.position = 'none', panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line = element_line(colour = "black"), panel.spacing.x = unit(.05, "lines"), panel.background = element_blank()) +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
  dev.off()
  
  #pad genes with few points so they don't get cut off in facets
  output$GENE <- as.character(output$GENE)
  for (gene in output[!rev(duplicated(rev(output$GENE))),'GENE']) {
    vec_length = length(filter(output, GENE %in% gene)[,'MUTS'])
    if (vec_length < 10) {
      short = 10 - vec_length
      for (i in c(1:round(short/2))) {
        output[nrow(output)+1,] <- 0
        output[output$MUTS==0,"MUTS"] = gene
        output$MUTS <- as.character(output$MUTS)
        output[output$MUTS==gene,"GENE"] = gene
        output$GENE <- as.character(output$GENE)
        output[output$MUTS==gene,"ci_lo"] = 0
        output[output$MUTS==gene,"ci_hi"] = 0
        output[output$MUTS==gene,"COEF"] = 0
        output[output$MUTS==gene,"POSITION"] = -500
      }
      for (i in c(1:(10- vec_length - round(short/2)))) {
        output[nrow(output)+1,] <- 0
        output[output$MUTS==0 & output$POSITION!=-500,"MUTS"] = gene
        output$MUTS <- as.character(output$MUTS)
        output[output$MUTS==gene & output$POSITION!=-500,"GENE"] = gene
        output$GENE <- as.character(output$GENE)
        output[output$MUTS==gene & output$POSITION!=-500,"ci_lo"] = 0
        output[output$MUTS==gene & output$POSITION!=-500,"ci_hi"] = 0
        output[output$MUTS==gene & output$POSITION!=-500,"COEF"] = 0
        output[output$MUTS==gene & output$POSITION!=-500,"POSITION"] = 10000
      }
    }
  }
  output[output$MUTS==output$GENE,'MUTS'] = ""
  output <- output[order(output$GENE,output$POSITION,-output$COEF),]
  
  #plot significant mutations in order of effect with CIs
  output$x <- NA
  for (gene in output[!rev(duplicated(rev(output$GENE))),'GENE']) {
    proportion = (length(output[output$GENE==gene,'x']) / length(output$x))
    output[output$GENE==gene,'x'] <- c(1:length(output[output$GENE==gene,'MUTS'])) * (5000/length(output[output$GENE==gene,'MUTS'])) * proportion
  }
  
  png(paste(path,'R/',drug,'/',drug,'_',path4,'_main_effects_by_gene.png',sep=''), width = 4096, height =1536)
  print(ggplot(output, aes(x=x, y=COEF)) + 
          facet_grid(.~GENE, scales="free_x", space="free_x") +
          geom_hline(yintercept=round(min(output$ci_lo)):round(max(output$ci_hi)), color='lightgrey') +
          #geom_hline(yintercept=log2(as.numeric(ecoffs[ecoffs$DRUG==drug,'ECOFF_CONC'])), color='darkgrey', linetype='dashed') +
          geom_hline(yintercept=0, color='darkgrey') +
          #geom_vline(xintercept=lines, colour='green', linetype="dashed") +
          geom_pointrange(aes(ymin=ci_lo, ymax=ci_hi, fatten=sqrt(abs(COUNTS))), color=ifelse(output$ci_lo>0,'red',ifelse(output$ci_hi<0,'blue','grey')), alpha=ifelse(output$COEF==0,0,0.3)) +
          geom_text(size = 4, aes(label=MUTS,hjust=0, vjust=0,angle=90)) + #ifelse(output$ci_lo>output[output$muts=='INT','ci_hi'],as.character(output$muts),ifelse(output$ci_hi<output[output$muts=='INT','ci_lo'],as.character(output$muts),''))
          labs(x='Mutations', y='Main Effect with 95% CI') +
          theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line = element_line(colour = "black"), panel.spacing.x = unit(.05, "lines"), panel.background = element_blank()))
  dev.off()
}