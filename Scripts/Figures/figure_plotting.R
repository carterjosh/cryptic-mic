library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(comprehenr)

se = 'normal'
path = ifelse(se=='normal','_normalSE','')
for (cluster in c('cluster100')) { #'cluster12','cluster25','cluster50',
for (drug in c('RIF','RFB','INH','EMB','ETH','KAN','AMI','BDQ','CFZ','LEV','LZD','MXF','DLM')) {
  #create directory if it does not exist
  dir.create(paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/no_int/R/',cluster,sep=''))
  dir.create(paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/no_int/R/',cluster,'/',drug,sep=''))
  
  #read in input data file prepared in Jupyter to get mutation names
  dat <- read.csv(file = paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Input_files/',drug,'_no_interaction_equal_end.csv',sep=''), header = TRUE)
  
  #read in STATA output from excel and clean labels etc. for dataframe
  out <- data.frame(read.csv(paste("/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/no_int/STATA/equal_end_range/",cluster,'/',drug,'_no_interaction_equal_end','.csv',sep=''),header = TRUE)) #,path
  addmuts <- c('Lineage 1','Lineage 2','Lineage 3','Lineage 4','Lineage 6','Site 2','Site 3','Site 4','Site 5','Site 6','Site 8','Site 10','Site 11','Site 14','Site 20','Baseline','var(_cons)','var(e.LOG2MIC)')
  addgene <- c('Lineage', 'Lineage', 'Lineage', 'Lineage','Lineage', 'Site', 'Site', 'Site', 'Site', 'Site', 'Site', 'Site', 'Site', 'Site', 'Site', 'Baseline', 'var', 'var')
  addpos <- c('1','2','3','4','6','2','3','4','5','6','8','10','11','14','20','21','22','23')
  muts <- colnames(dat)
  muts <- muts[1:(length(muts)-10)]
  temp <- regmatches(muts, regexpr("_", muts), invert = TRUE)
  mat  <- matrix(unlist(temp), ncol=2, byrow=TRUE)
  df   <- as.data.frame(mat)
  colnames(df) <- c('GENE','POSITION')
  df$POSITION <- gsub("[^0-9.]", "", df$POSITION)
  df$POSITION <- as.numeric(gsub("^\\.", "-", df$POSITION))
  tmp <- data.frame(addgene,addpos)
  colnames(tmp) <- c('GENE','POSITION')
  df <- rbind(df,tmp)
  out$muts <- append(muts,addmuts,after=length(muts))
  
  #create dataframe and confidence intervals
  output <- data.frame(out$muts,df$GENE,as.numeric(df$POSITION),out$coef,out$stderr,out$tstat)
  colnames(output) <- c('MUTS','GENE','POSITION','COEF','SE','Z')
  output$corrected_se <- output$SE*1.96  #*qnorm(1-.05/(length(muts)-1))
  output$ci_hi <- output$COEF+output$corrected_se
  output$ci_lo <- output$COEF-output$corrected_se
  output$PVAL <- 2*pnorm(-abs(output$Z))
  output$BH_PVAL <- p.adjust(output$PVAL, method = "BH")
  
  #create main effects output dataframe, sort by gene, position, and effect and write to file
  output <- output[order(output$GENE,output$POSITION,-output$COEF),]
  output <- output[complete.cases(output),]
  write.csv(output, file = paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/no_int/R/',cluster,'/',drug,'/',drug,'_main_effects',path,'.csv',sep=''),row.names=FALSE)
  tmp <- filter(output, grepl('indel', MUTS) & POSITION>0)
  tmp$POSITION <- tmp$POSITION / 3
  output <- output[!grepl("indel", output$MUTS),]
  output <- rbind(output,tmp)
  output <- output[order(output$GENE,output$POSITION,-output$COEF),]
  output <- output[output$MUTS!='Baseline',]
  output <- output[output$MUTS!='var(_cons)',]
  output <- output[output$MUTS!='var(e.LOG2MIC)',]
  output$MUTS <- gsub("[.]", "-", output$MUTS)
  output$MUTS <- gsub('-$', "!", output$MUTS)
  #removing anything with absurd standard error
  output <- output[output$corrected_se<10,]
  #remove mutations that were omitted due to collinearity
  output <- output[output$COEF!=0,]
  #keep only significant effects for plotting
  output <- output[!(output$ci_hi>0 & output$ci_lo<0),]
  output <- output[output$BH_PVAL<0.05,]
  
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
  
  genes <- data.frame(output[!rev(duplicated(rev(output$GENE))),'GENE'], stringsAsFactors = FALSE)
  colnames(genes) <- c('GENE')
  genes <- filter(genes, !(GENE %in% c('Lineage','Site')))[['GENE']]
  genes <- append(genes,c('Lineage','Site'),after=length(genes))
  output$GENE <- factor(output$GENE, levels = genes)

  #output$x <- c(1:length(output$MUTS)) * (2500/length(output$MUTS))
  #lines <- output[!rev(duplicated(rev(output$GENE))),'x'] + 5
  
  gene_labels <- output[!duplicated(output$GENE),'GENE']
  label_colors <- to_vec(for (i in 1:length(gene_labels)) ifelse (i %% 2 == 0, 'blue', 'grey'))
  
  png(paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/no_int/R/',cluster,'/',drug,'/',drug,'_main_effects',path,'.png',sep=''), width = 4096, height =1536)
  print(ggplot(output, aes(x=output$x, y=output$COEF)) + 
          facet_grid(.~GENE, scales="free_x", space="free_x") +
          geom_hline(yintercept=round(min(output$ci_lo)):round(max(output$ci_hi)), color='lightgrey') +
          #geom_hline(yintercept=log2(as.numeric(ecoffs[ecoffs$DRUG==drug,'ECOFF_CONC'])), color='darkgrey', linetype='dashed') +
          geom_hline(yintercept=0, color='darkgrey') +
          #geom_vline(xintercept=lines, colour='green', linetype="dashed") +
          geom_pointrange(aes(ymin=output$ci_lo, ymax=output$ci_hi), color=ifelse(output$ci_lo>0,'red',ifelse(output$ci_hi<0,'blue','grey')), alpha=ifelse(output$COEF==0,0,0.3)) +
          geom_text(size = 4, aes(label=output$MUTS,hjust=0, vjust=0,angle=90)) + #ifelse(output$ci_lo>output[output$muts=='INT','ci_hi'],as.character(output$muts),ifelse(output$ci_hi<output[output$muts=='INT','ci_lo'],as.character(output$muts),''))
          labs(x='Mutations', y='log2(MIC) Change') +
          scale_y_continuous(breaks = round(min(output$ci_lo)):round(max(output$ci_hi))) +
          theme(panel.border = element_blank(), panel.background = element_blank(),
              panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line = element_line(colour = "black"),
                axis.ticks.x=element_blank(), axis.text.x=element_blank(),panel.spacing.x = unit(0, "lines"),
                strip.background = element_rect(color="black", fill=label_colors, size=1.5, linetype="solid")))
  dev.off()
}
}