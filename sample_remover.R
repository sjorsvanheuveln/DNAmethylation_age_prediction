#Sample Remover from Local MarmalAid
#SvH; 02-04-2016

#Info: Removes sample files that are not blood tissue data

setwd('~/Desktop/mBioinformatics/2.2 Forensic Epigenetics/Local_Methylation_DB')
samples <- read.csv('Samples.csv')
non_blood_IDs <- samples[samples$TISSUE!='Blood' | samples$DISEASE!='Healthy','Id']
files <- paste(non_blood_IDs,'_nm_imp.Rdata',sep='')
files2 <- paste(non_blood_IDs,'.Rdata',sep='')
files <- c(files,files2)
remove_list <- files[files %in% list.files('./Rdata/')]

for (file in remove_list){
  print(paste(file,'removed'))
  unlink(paste('./Rdata/',file,sep=''))
}
