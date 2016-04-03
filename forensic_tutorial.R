#MARMAL-AID
#SvH; 02-04-2016

#TODO:  1) Check whether I can put the local depository to the GeekCloud 
#       2) Make a function for feature selection
#LC:    1) Server version maximally takes only 5000 probes instead of 450k!
#       2) Files komen in paren, sample.Rdata & sample_nm_imp.Rdata
#       3) Filtering on regression coefficients seems logical, because it's not only high variance
#          but also a relationship with age.
#       4) Alle TCGA kan er uitgetieft worden, is allemaal kanker samples.
#       5) "Beta" is a confusiong name for "fraction methylated".

#### FUNCTIONS ####
load_beta <- function(use_local,n){
  #get data from server or local drive;n is the number of samples you want
  data(meth450)
  meth450 = subset(meth450,meth450$CHR!= 'X' & meth450$CHR!= 'Y')
  probes = row.names(meth450)

  if (use_local == T){
    #only select samples that are in directory and have an Age annotation
    samples = list.files('./Local_Methylation_DB/Rdata')
    samples = unique(gsub('_nm_imp','',gsub('.Rdata','',samples)))
    samples = annotation[samples,][!is.na(annotation[samples,"Age"]),1][1:n]
    
    #somehow loading problems with certain files, due to the nm_inp suffix
    beta=getbeta(samples,probes,marmalaiddir="~/Desktop/mBioinformatics/2.2 Forensic Epigenetics/Local_Methylation_DB");beep(4)
  }else{
    samples = annotation[which(annotation$DISEASE=="Healthy" & annotation$TISSUE=="Blood" & !is.na(annotation$Age)),1][1:n]
    u=url("http://marmal-aid.org/probes.Rdata")
    load(u)
    close(u)
    beta=getbeta(samples[1:5],probes);beep(4)
  }
  return(beta)
}

#### MAIN ####
setwd('~/Desktop/mBioinformatics/2.2 Forensic Epigenetics')
library(marmalaid)
library(IlluminaHumanMethylation450k.db)
library(beepr)

#Load beta and age for n samples
data(annotation_v1.1)
annotation$Age <- as.numeric(annotation$Age);print('Coerced NAs are OK!') #overcomes problems later when selecting samples;drop about a 100 healthy-bloods
beta <- load_beta(use_local = TRUE,n=200);rm(meth450)
annotation <- annotation[annotation$Id %in% colnames(beta),]
age = annotation[colnames(beta),"Age"] #Age vector
hist(age)
CpG_annotation <- as.list(IlluminaHumanMethylation450kSYMBOL[mappedkeys(IlluminaHumanMethylation450kSYMBOL)])

############################
##### Feature Selection ####
############################

fs_methods = c('variance','correlation')
fs_choice = fs_methods[2]

if (fs_choice == 'variance'){
  #prefilter on variance to reduce amount of CpG sites
  #assume that age correlated sites also have relative high variance
  print('Feature selection by "variance"')
  CpG_site_variance <- apply(beta,1,var)
  hist(CpG_site_variance)
  CpG_subset <- subset(CpG_site_variance,CpG_site_variance > 0.01)
  beta_sub <- beta[row.names(beta) %in% names(CpG_subset),]

}else if(fs_choice == 'correlation'){
  CpG_site_correlation <- apply(beta,1, function(x) cor(age,x))
  hist(abs(CpG_site_correlation))
  CpG_subset <- subset(CpG_site_correlation,abs(CpG_site_correlation) > 0.8)
  beta_sub <- beta[row.names(beta) %in% names(CpG_subset),]
  print('Feature selection by "correlation"')

}else{print('Feature selection by "other"')}

##########################
##### Model Building #####
##########################


#plot n CpG-sites over all samples
for (n in 1:10){
  Sys.sleep(1)
  #n=5
plot(age,y=beta_sub[n,],ylim=c(0,1),col='blue',main=paste(row.names(beta_sub)[n],CpG_annotation[names(CpG_subset)][n]),xlab='Age',ylab='fraction methylated')
}










