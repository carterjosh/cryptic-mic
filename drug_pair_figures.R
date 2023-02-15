library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(plotly)
library(ggrepel)
library(comprehenr)

drugA = 'RIF'
drugB = 'RFB'
drug_pair = paste(drugA,'-',drugB,sep='')
drug_epi <- FALSE
int <- FALSE
no_int <- TRUE
pop_control <- FALSE

se = 'normal'
#path = ifelse(se=='normal','_normalSE','')
path <- '/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/no_int/'
path3 <- 'STATA/equal_end_range/cluster100/'
for (cluster in c('cluster100')) {
  #for (drug in c('RIF','RFB','INH','ETH','EMB','KAN','AMI','BDQ','CFZ','LEV','LZD','MXF','DLM')) { 
    #create directory if it does not exist
    dir.create(paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/no_int/R/',cluster,sep=''))
    dir.create(paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/no_int/R/',cluster,'/drug_pairs',sep=''))
    
    count=0
    for (drug in c(drugA,drugB)) {
      #read in input data file prepared in Jupyter to get mutation names
      dat <- read.csv(file = paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Input_files/',drug,'_no_interaction_equal_end_short_final.csv',sep=''), header = TRUE)
      
      #read in STATA output from excel and clean labels etc. for dataframe
      out <- data.frame(read.csv(file = paste(path,path3,drug,'_no_interaction_equal_end.csv',sep='')))
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
      output <- data.frame(out$var,df$GENE,as.numeric(df$POSITION),out$coef,out$stderr,out$tstat,out$pval)
      colnames(output) <- c('MUTS','GENE','POSITION','COEF','SE','Z','P')
      output$corrected_se <- output$SE*1.96  #*qnorm(1-.05/(length(muts)-1))
      output$ci_hi <- output$COEF+output$corrected_se
      output$ci_lo <- output$COEF-output$corrected_se
      output$BH_PVAL <- p.adjust(output$P, method = "BH")
      
      #catch stop codon mutations
      output[(abs(output$POSITION)>100) & (str_detect(output$MUTS, '-')),'MUTS'] <- gsub("-", "!", output[(abs(output$POSITION)>100) & (str_detect(output$MUTS, '-')),'MUTS'])
      output[abs(output$POSITION)>100,'POSITION'] = abs(output[abs(output$POSITION)>100,'POSITION'])
      
      
      #create main effects output dataframe, sort by gene, position, and effect and write to file
      output <- output[order(output$GENE,output$POSITION,-output$COEF),]
      output <- output[complete.cases(output),]
      #write.csv(output, file = paste(path,'R/',drug,'/',drug,'_',path4,'_main_effects.csv',sep=''),row.names=FALSE)
      tmp <- filter(output, grepl('indel', MUTS) & POSITION>0)
      tmp$POSITION <- tmp$POSITION / 3
      output <- output[!grepl("indel", output$MUTS),]
      output <- rbind(output,tmp)
      output <- output[order(output$GENE,output$POSITION,-output$COEF),]
      output <- output[output$MUTS!='ns[cluster])',]
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
      
      if (count==0) {outputA <- output} else {outputB <- output}
      count = count+1
    }
    
    outputB$COEF2 <- outputB$COEF
    outputB$ci_hi2 <- outputB$ci_hi
    outputB$ci_lo2 <- outputB$ci_lo
    outputB <- outputB[c('MUTS','COEF2','ci_hi2','ci_lo2')]
    output <- merge(outputA, outputB, by="MUTS")
    
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
    
    png(paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/no_int/R/',cluster,'/drug_pairs/',drug_pair,'_main_effects.png',sep=''), width = 4096, height =1536)
    print(ggplot(output, aes(x=output$COEF2, y=output$COEF, lab)) + 
            facet_grid(.~GENE, scales="free_x", space="free_x") +
            geom_hline(yintercept=round(min(output$ci_lo)):round(max(output$ci_hi)), color='lightgrey') +
            #geom_hline(yintercept=log2(as.numeric(ecoffs[ecoffs$DRUG==drug,'ECOFF_CONC'])), color='darkgrey', linetype='dashed') +
            geom_hline(yintercept=0, color='darkgrey') +
            geom_abline(slope=1,intercept=0, color="darkgrey") +
            geom_abline(slope=1,intercept=-1, linetype="dashed", color="lightgrey") +
            geom_abline(slope=1,intercept=1, linetype="dashed", color="lightgrey") +
            #geom_vline(xintercept=lines, colour='green', linetype="dashed") +
            geom_pointrange(aes(ymin=output$ci_lo, ymax=output$ci_hi), color=ifelse(output$ci_lo>0,'red',ifelse(output$ci_hi<0,'blue','grey')), alpha=ifelse(output$COEF==0,0,ifelse(abs(output$COEF-output$COEF2)>=1,0.4,0.1))) +
            geom_errorbarh(aes(xmin=output$ci_lo2, xmax=output$ci_hi2), color=ifelse(output$ci_lo>0,'red',ifelse(output$ci_hi<0,'blue','grey')), alpha=ifelse(output$COEF==0,0,ifelse(abs(output$COEF-output$COEF2)>=1,0.4,0.1))) +
            #geom_text(size = 4, aes(label=output$MUTS,hjust=0, vjust=0,angle=90)) + #ifelse(output$ci_lo>output[output$muts=='INT','ci_hi'],as.character(output$muts),ifelse(output$ci_hi<output[output$muts=='INT','ci_lo'],as.character(output$muts),''))
            labs(x=paste0(drugB,' LOG2MIC'), y=paste0(drugA,' LOG2MIC')) +
            scale_y_continuous(breaks = round(min(output$ci_lo)):round(max(output$ci_hi))) +
            geom_text_repel(aes(label = ifelse((abs((output$COEF-output$COEF2))>=1),output$MUTS,''))) +
            theme(panel.border = element_blank(), panel.background = element_blank(),
                  panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.line = element_line(colour = "black"),
                  panel.spacing.x = unit(0.5, "lines"),
                  strip.background = element_rect(color="black", fill=label_colors, size=1.5, linetype="solid")))
    dev.off()
#  }
}