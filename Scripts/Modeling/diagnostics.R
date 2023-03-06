library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(survival)
library(lmtest)
library(lawstat)

no_int <- TRUE
int <- FALSE
mixed <- FALSE

if (mixed==FALSE) {path3 = 'int_reg/'}
if (mixed==TRUE) {path3 = ''}

if (no_int==TRUE) {path2 = '_no_interaction_equal_end.csv'}
if (int==TRUE) {path2 = '_interaction_equal_end.csv'}

if (no_int==TRUE) {drug_list <- c('RIF','RFB','INH','ETH','EMB','KAN','AMI','BDQ','CFZ','LEV','LZD','MXF','DLM')}
if (int==TRUE) {drug_list <- c('INH','ETH','EMB','KAN','LEV','MXF')}

for (drug in drug_list) {
if (no_int==TRUE) {path = paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/no_int/R/',path3,drug,'/diagnostics/',sep='')}
if (int==TRUE) {path = paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/int/R/',path3,drug,'/diagnostics/',sep='')}
  
#read in input file with mutations, UNIQUEIDs, LINEAGE, SITE, and cluster info
dat <- read.csv(file = paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Input_files/',drug,path2,sep=''), header = TRUE)

#read in fitted predictions
if (no_int==TRUE) {out = read.csv(file = paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/no_int/STATA/equal_end_range/cluster100/',path3,drug,'_fitted_predictions.csv',sep=''), header = TRUE)}
if (int==TRUE) {out = read.csv(file = paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/int/STATA/cluster100/',path3,drug,'_fitted_predictions.csv',sep=''), header = TRUE)}


#generate midpoints of interval ranges
midpoints <- out$uplog2mic-((out$uplog2mic - out$log2mic)/2)

#create psuedo-r^2 to see how much variation we are explaining
pseudo_r_squared <- with(out, cor(cbind(eta_fix,eta_cond, log2mic, uplog2mic, midpoints)))^2
dir.create(paste(path,sep=''))
write.csv(pseudo_r_squared,paste(path,drug,'_correlations.csv',sep=''))

#plot predicted vs actual (midpoint)
png(paste(path,drug,'_pred_vs_actual.png',sep=''))
print(ggplot(out,aes(eta_cond,midpoints)) +
  geom_point(shape = 21, alpha = 0.3, size = 3) + #position=position_jitter(h=0.1, w=0.1), 
  geom_abline() +
  xlab('Predicted') + ylab('Actual (midpoint)') + xlim(-15,10) + ylim(-15,10) + theme_bw())
dev.off()


#plot studentized residuals from interval midpoint
png(paste(path,drug,'_studentized_resid.png',sep=''))
print(ggplot(out,aes(eta_cond,(midpoints-eta_cond)/sd(midpoints-eta_cond))) +
  geom_point(shape = 21, alpha = 0.3, size = 3) + #position=position_jitter(h=0.1, w=0.1), 
  geom_hline(yintercept = 0) +
  xlab('Predicted') + ylab('Studentized residuals (midpoint)') + theme_bw())
dev.off()

#qqplot of midpoint residuals
png(paste(path,drug,'_qqplot.png',sep=''))
qqnorm(midpoints-out$eta_cond)
qqline(midpoints-out$eta_cond)
dev.off()

#histogram of studentized residuals
png(paste(path,drug,'_resid_histo.png',sep=''))
hist((midpoints-out$eta_cond)/sd(midpoints-out$eta_cond), xlab = 'Studentized Residuals', main = 'Studentized Residual Distribution')
dev.off()

#regress pred vs actual and run likelihood ratio test to ensure random effects are beneficial
#linear_full <- lm(log2mic ~ eta_cond - 1, data=out)
#summary(linear_full)
#linear_fix <- lm(log2mic ~ eta_fix - 1, data=out)
#summary(linear_fix)

#full_fit<-survreg(Surv(out$log2mic, out$uplog2mic, event = rep(3, nrow(out)),type='interval')~out$eta_cond - 1, dist='gaussian')
#summary(full_fit)
#fix_fit<-survreg(Surv(out$log2mic, out$uplog2mic, event = rep(3, nrow(out)),type='interval')~out$eta_fix - 1, dist='gaussian')
#summary(fix_fit)

#lrtest(linear_fix, linear_full)
#lrtest(fix_fit, full_fit)
}

for (drug in drug_list) {

dat <- read.csv(file = paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Input_files/',drug,path2,sep=''), header = TRUE)
  
#read in fitted predictions
if (no_int==TRUE) {out = read.csv(file = paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/no_int/STATA/equal_end_range/cluster100/',path3,drug,'_fitted_predictions.csv',sep=''), header = TRUE)}
if (int==TRUE) {out = read.csv(file = paste('/Users/joshuacarter/Dropbox/cryptic-analysis/Interval_regression/Output_files/int/STATA/cluster100/',path3,drug,'_fitted_predictions.csv',sep=''), header = TRUE)}
  
  
#generate midpoints of interval ranges
midpoints <- out$uplog2mic-((out$uplog2mic - out$log2mic)/2)
  
  
tmp_num <- match("SITEID",colnames(dat))
dat$LINEAGE <- as.factor(dat$LINEAGE)
dat$SITEID <- as.factor(dat$SITEID)

test <- lm((midpoints-out$eta_fix)^2 ~ ., data = dat[1:tmp_num])
test_hets <- data.frame(summary(test)$coef[summary(test)$coef[,4] <= .05, 4])
paste(tolower(row.names(test_hets)), collapse=" ")

hets <- c()
for (mut in colnames(dat[1:tmp_num])) {
  tmp <- levene.test(midpoints-out$eta_fix,dat[,mut])
  if (tmp$p.value<0.05) {
    hets <- append(hets,tolower(mut))
  }
}
paste(hets, collapse = " ")

}
