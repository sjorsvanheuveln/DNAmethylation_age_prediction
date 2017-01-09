#Functions for MFM.R
#12-10-2016

#### Functions ####
## Main Functions ##
load.preset <- function(beta,type='new'){
  new_set <- c("age","cg02872426","cg22156456","cg16867657","cg16054275","cg12934382","cg18473521","cg08097417","cg08262002","cg06874016","cg25410668","cg03224418","cg07553761","cg11807280","cg13959344","cg06784991")
  old_set <- c("age","cg02872426","cg22156456","cg16867657","cg16054275","cg12934382","cg18473521","cg07955995","cg08262002","cg06874016","cg25410668","cg03224418","cg07553761","cg11807280","cg13959344","cg16015712")#,"cg02694715","cg10007452"
  if (type=='lasso'){
    print('Selecting Feature PreSet with LASSO')
    beta <- NA.replace(beta,'median')
    set.seed(1) #required as cv.glmnet creates random folds
    registerDoMC(cores=3)
    #maybe iterate here 10 times on cv.glmnet and catch only the overlapping features?
    CV <- cv.glmnet(x=as.matrix(beta[,-1]),y=beta$age,alpha=1,standardize=T,parallel=T)
    plot(CV)
    sites <- order(coef(CV,CV$lambda.1se))
    sites <- sites@Dimnames[[1]][sites@i+1][-1]
    beep(28)
    return(sites)}
  else if (type=='old'){return(old_set)}
  else if (type=='new'){return(new_set)}
  else if (type=='all'){return(colnames(beta))}
  else{return(new_set)}
}
load_data <- function(data_name = 'regular_equal',cap=T,exclude_outliers=F,age_low,age_high,preset){
  files <- list.files(path = './Data/Methylation_data/',recursive=T)
  print(paste('Loading',data_name))
  
  if (data_name == 'Hannum' | data_name == 'hannum'){
    load(paste('./Data/Methylation_data/',files[9],sep=''))
    load(paste('./Data/Methylation_data/',files[10],sep=''))
    beta <- rbind(test_Ha,train_Ha)
    beta$sex <- NULL; beta$ethnicity <- NULL
  }else if (data_name =='regular_equal'){
    beta <- load(paste("Data/Methylation_data/other dataset/beta_wb_eq_RF.Rdata"))
    beta <- get(beta)
    row.names(beta) <- unlist(lapply(strsplit(row.names(beta),'_'),function(x){x[1]}))
  }
  if ((data_name == 'Hannum' | data_name == 'hannum') & preset=='all'){beep(28);return(beta)}
  #subset data
  PreSet <- load.preset(beta,preset)
  beta_preset <- beta[beta$age > age_low & beta$age <= age_high,PreSet] #betadata from biomarker preselection
  beta_preset <- NA.replace(beta_preset,'median')
  if (cap){beta_preset <- beta_preset[cap_samples(beta_preset,5),]}
  
  #remove outliers
  if (exclude_outliers){
    print('detecting outliers ...')
    outliers <- get_outlier_wrapper(beta_preset,n,k,missing,poly,BFE)
    print(paste('removing',length(outliers[[1]]),'outliers'))
    beta_preset <- beta_preset[which(!row.names(beta_preset) %in% outliers[[1]]),]
  }
  beep(28)
  return(beta_preset)
}
train_model <- function(beta_train,model_type,poly,BFE){
  if (model_type=='MLR'){model <- lm(age ~., data=beta_train)}
  else if (model_type=='MLR_polyELOVL'){if ('cg16867657' %in% colnames(beta_train)){model <- lm(age ~ poly(cg16867657,poly) + . - cg16867657, data=beta_train)}else if ('fake' %in% colnames(beta_train)){model <- lm(age ~ poly(fake,poly) + . - fake, data=beta_train)}else{model <- lm(age ~., data=beta_train)}}
  else if (model_type=='MLR_polyAll'){model <- lm(as.formula(paste('age ~',paste('poly(',colnames(beta_train[-1]),',',poly,')',collapse = ' + '))),data=beta_train)}
  else if (model_type=='RF'){model <- randomForest(x=beta_train[,-1],y=beta_train[,1],ntree=100,replace=T,importance=F,localImp=F,nPerm=3,corr.bias=F,keep.inbag=F,norm.votes=F)}
  else if (model_type=='SVR'){#tuneResult <- tune(svm, age ~ .,  data = beta_train, ranges = list(epsilon = seq(0.1,0.3,0.05), cost = seq(0.2,1,0.2), kernel='linear'));
                              model <- svm(age ~., data=beta_train, epsilon = 0.16, cost = 1, kernel='linear')}
  else if (model_type=='PCR'){model <- pcr(age ~., data=beta_train, validation = "CV", scale=F)}
  else if (model_type=='PLS'){model <- plsr(age ~., data=beta_train, validation = "CV", scale=F)}
  else if (model_type=='MARS'){model <- earth(age ~ ., data=beta_train, nk=12, nprune=8)} #10,5 for normal, 13/12,8 hannum
  else if (model_type=='GAM'){model <- gamboost(age ~ .,data=beta_train)}
  else{print('WARNING! No valid model chosen!');stop()}
  if (BFE){model <- stepAIC(model,direction='backward',trace=F)}
  return(model)
}
get_models <- function(printmodels=T,exclude_model=NULL){
  models <- c('MLR','MLR_polyELOVL','MLR_polyAll','RF','SVR','PCR','PLS','MARS','GAM') 
  models <- models[!models %in% exclude_model]
  if (printmodels){
    for (i in 1:length(models)){
      print(paste(i,': ',models[i]),sep='')
    }
  }
  return(models)
}
simulate_plot <- function(deltaMAD,model_id,n,k,AFM.performance,beta_preset,AFM.prediction,AFM.errors,poly){
  ## Plot MFM results ##
  #calclulate statistics
  shapiro.pvalue <- unlist(lapply(deltaMAD,function(x){if (diff(range(x))==0){NA}else{shapiro.test(x)$p.value}})) #test for normal distribution
  shapiro.padj <- round(p.adjust(shapiro.pvalue,method='BH'),4) #shapiro.test is not normal when significant
  wilcox.pvalue <- unlist(lapply(deltaMAD,function(x){wilcox.test(x,mu=0)$p.value})) #non-parametric t-test
  wilcox.padj <- round(p.adjust(wilcox.pvalue,method='BH'),4)
  ttest.pvalue <- unlist(lapply(deltaMAD,function(x){t.test(x,mu=0)$p.value})) #parametric t-test
  ttest.padj <- round(p.adjust(ttest.pvalue,method='BH'),4)
  
  #plot
  print('Plotting MFM')
  layout(matrix(c(1,1,1,2,2,3,3,3), nrow = 8, ncol = 1, byrow = T)) #indicate per plot how much rows it can take
  par(mar=c(1,5,2,1))
  boxplot(deltaMAD,las=2,col='red',ylim=c(-2,4),main=paste('deltaMAD ',model_id,' (n=',n,', k=',k,')',sep=''),ylab='deltaMAD',xaxt="n")
  abline(h=0,lty=2)
  barplot(shapiro.padj,col='blue',main=paste('Shapiro-Wilk p.adjusted BH'),ylim=c(0,1),xaxt="n",ylab='p.value',las=2);abline(h=0.05,lty=2)
  par(mar=c(7.5,5,1,1))
  #barplot(wilcox.padj,col='purple',main=paste(paste('Wilcoxon p.adjusted BH'),model_id),ylim=c(0,1),ylab='p.value',las=2);abline(h=0.05,lty=2)
  barplot(rbind(wilcox.padj,ttest.padj),beside=T,col=c('orange','purple'),las=2,main=paste(paste('Wilcoxon/ttest p.adjusted BH')),ylim=c(0,1),ylab='p.value');abline(h=0.05,lty=2)
  legend("topright",fill=c("orange","purple"), c("Wilcox", "ttest"),bty='n',horiz = TRUE,cex=0.8)
  
  #copy to file and reset plot parameters
  path <- paste('/Images/',model_id,'MFM.pdf',sep='');dev.copy2pdf(file = paste('.',path,sep=''),width = 5, height = 7);dev.off();par(mfrow=c(1,1),mar=c(4,4,2,1))
  
  ## Plot AFM performance ##
  
  RMSE <- lapply(AFM.errors,rmse); RMSE.mean <- round(mean(unlist(RMSE)),2)
  MAD <- lapply(AFM.errors,mae); MAD.mean <- round(mean(unlist(MAD)),2)
  
  predictions <- unlist(AFM.prediction)
  predictions <- split(predictions,names(predictions))
  pred_df <- data.frame(matrix(unlist(predictions),nrow=nrow(beta_preset),byrow=T));row.names(pred_df) <- names(predictions)
  pred_df$age <- beta_preset[row.names(pred_df),'age']
  setDT(pred_df, keep.rownames = TRUE)[]
  pred_df=melt(pred_df,id.vars=c("rn","age"))
  
  print('Plotting AFM Performance')
  p <- ggplot(pred_df,aes(x=age,y=value,group=rn)) + geom_boxplot(color = 'black',position="identity",alpha=0.6) + theme_linedraw() + ylim(0, 90) + scale_y_continuous(breaks = seq(0, 100, 10)) #alpha = transparancy
  p + labs(list(title = paste('Performance',model_id),x='Real Age',y='Predicted Age')) + geom_abline(intercept=0,slope=1,col='red') + theme(plot.title = element_text(size=22)) +
    annotate("text", x= 20, y=78, label= paste("RMSE =",RMSE.mean,'±',ci(RMSE),'(95% CI)'),hjust = 0,size=6) +
    annotate("text", x= 20, y=72, label= paste("MAD =",MAD.mean,'±',ci(MAD),'(95% CI)'),hjust = 0,size=6)
  ggsave(filename=paste('./Images/',model_id,'_RvsP.pdf',sep=''), plot = last_plot(),width=7,height=7)
}
simulate <- function(beta_preset,n,k,model_type,missing,poly,BFE,bin_frac,bin_start){
  deltaMAD_list <- list()
  AFM.performance <- list()
  AFM.prediction <- list()
  AFM.errors <- list()
  set.seed(1)
  print(paste('Running Simulation for',model_type))
  
  for (j in 1:n){
    index <- createFolds(beta_preset[,1], k = k)
    deltaMAD_list[[j]] <- list()
    AFM.performance[[j]] <- vector()
    AFM.prediction[[j]] <- list()
    AFM.errors[[j]] <- list()
    
    for (i in 1:k){
      print(paste(j,i))
      beta_train <- beta_preset[-unlist(index[i]),]
      beta_test <-  beta_preset[unlist(index[i]),]
      
      #Train,Predict and Evaluate AFM
      AFM <- train_model(beta_train,model_type,poly,BFE)
      AFM.prediction[[j]][[i]] <- AFM_predict(AFM,model_type,beta_test)
      #do binned prediction if bin_frac > 0
      if (bin_frac > 0){AFM.prediction[[j]][[i]] <- binPredict(AFM.prediction[[j]][[i]],beta_train,beta_test,model_type,poly,BFE,bin_frac,bin_start)}
      AFM.errors[[j]][[i]] <- AFM.prediction[[j]][[i]]-beta_test$age
      AFM.mad <- mae(AFM.errors[[j]][[i]])
      #AFM.rmse <- rmse(AFM.errors)
      AFM.performance[[j]][i] <- AFM.mad
      
      #Train,Predict and Evaluate MFM
      if (missing){
        MFM.prediction <- MFM_predict(beta_train,beta_test,model_type,poly)
        MFM.errors <- lapply(MFM.prediction,function(x){x-beta_test$age})
        MFM.mad <- mae(MFM.errors)
        #MFM.rmse <- rmse(MFM.errors)
        
        #calculate deltaMAD
        deltaMAD <- MFM.mad-AFM.mad
        deltaMAD_list[[j]][[i]] <- deltaMAD}
    }
  }
  if (missing){
    deltaMAD_list <- unlist(deltaMAD_list)
    deltaMAD_list <- split(deltaMAD_list,unique(names(deltaMAD_list)))}
  
  #beep(13);Sys.sleep(1)
  return(list(deltaMAD_list,AFM.performance,AFM.prediction,AFM.errors))
}
MFM_predict <- function(beta_train,beta_test,model_type,poly){
  #train and predict with all MFM models (dependent on n/o features) and returns vector of predictions
  MFM_prediction <- list()
  for (feature in colnames(beta_train)[-1]){
    beta_missing <- beta_train[,-which(colnames(beta_train)==feature)]
    MFM <- train_model(beta_missing,model_type,poly,BFE)
    MFM_prediction[[feature]] <- predict(MFM,beta_test)
  }
  return(MFM_prediction)
}
AFM_predict <- function(AFM, model_type,beta_test){
  #a wrapper for predict() which allows generating the same predict output for different models
  if (model_type %in% get_models(F)[1:5]){prediction <- predict(AFM,beta_test)}
  else if (model_type %in% get_models(F)[6:7]){prediction  <- predict(AFM,newdata=beta_test)[,,which(AFM$validation$PRESS==min(AFM$validation$PRESS))]}
  else if (model_type=='MARS'){prediction <- as.vector(predict(AFM,newdata=beta_test)); names(prediction) <- row.names(beta_test)}
  else if (model_type=='GAM'){prediction <- predict(AFM,beta_test)}
  else{print('WARNING! No valid model chosen!');stop()}
  return(prediction)
}
binPredict <- function(pre_predictions,beta_train,beta_test,model_type,poly,BFE,bin_frac,bin_start){
  #uses a first prediction to re-predict based on closest (bin_range) age samples
  #binsize should be a fraction between 0 ≤ binsize ≤ 1
  bin_predictions <- c()
  bin_range <- floor(nrow(beta_train)*bin_frac) #set number of samples in bin-range
  for (i in 1:length(pre_predictions)){
    sample <- pre_predictions[i]
    if (sample < bin_start){bin_predictions[(length(bin_predictions)+1)] <- sample;next} # only bin for high age
    bin_samples <- order((beta_train$age - sample)^2)[1:bin_range]
    bin_train <- beta_train[bin_samples,]
    bin_model <- train_model(bin_train,model_type,poly,BFE)
    bin_predictions[(length(bin_predictions)+1)] <- AFM_predict(bin_model,model_type,beta_test[names(sample),])
  }
  names(bin_predictions) <- names(pre_predictions)
  return(bin_predictions)
}

## Wrappers ##
bin_performance_wrapper <- function(beta_preset,n,k,model_type,missing,poly,BFE,bin_start){
  #wrapper for simulate function to optimize bin_frac and age-cut off
  #iterate over models?
  bin_frac_vector <- seq(0.1,1,0.1)
  performances <- list()
  for (bin_frac in bin_frac_vector){
    performances[[as.character(bin_frac)]] <- unlist(simulate(beta_preset,n,k,model_type,missing,poly,BFE,bin_frac,bin_start)[[2]])
  }
  pred_df <- data.frame(matrix(unlist(performances),nrow=100,byrow=F));colnames(pred_df) <- names(performances)
  pred_df$id <- seq(1,100,by=1)
  pred_df <- melt(pred_df, id=("id"))
  
  p <- ggplot(pred_df,aes(x=as.numeric(as.character(variable)),y=value,group=variable)) + geom_boxplot(color = 'black',position="identity",alpha=0.6) + theme_light() + ylim(2,6)#alpha = transparancy
  p <- p + labs(list(title = paste(model_type,': BinPredict at ',bin_start,' years',sep=''),x='bin_fraction',y='MAD')) + theme(plot.title = element_text(size=22)) +
    stat_summary(fun.y=mean, geom="line", aes(group=1),col='red') + geom_hline(aes(yintercept=mean(performances$`1`)),col='blue',linetype="dashed")
  
  t.test(performances[[which.min(lapply(performances,mean))]],performances$`1`)$p.value
  print(paste('deltaMAD =',round(min(unlist(lapply(performances,mean))- mean(performances$`1`)),4),'at bin_frac =',names(which.min(lapply(performances,mean)))))
  return(p)
}
performance_wrapper <- function(beta_preset,n,k,strat,poly,BFE,exclude_model=NULL,doprint=F){
  #gives mean MAD/RMSE and confidence interval for all models
  
  #stratification methods
  if (strat=='random'){ages <- 1:nrow(beta_preset);if (doprint){print('Random Folds')}}
  else if (strat=='nearest5'){ages <- mround(beta_preset$age,5);if (doprint){print('Stratifying by closest 5 Folds')}}
  else if (strat=='nearest10'){ages <- mround(beta_preset$age,10);if (doprint){print('Stratifying by closest 10 Folds')}}
  else{ages <- beta_preset[,1];if (doprint){print('Stratifying by Age Folds')}} #normal mode;
  
  modelMAD <- list()
  for (model in get_models(F,exclude_model)){
    
    set.seed(1)
    MAD.performance <- list()
    RMSE.performance <- list()
    predictions <- list()
    
    for (i in 1:n){
      index <- createFolds(ages, k = k) #create folds
      MAD.performance[[i]] <- vector()
      RMSE.performance[[i]] <- vector()
      
      for (j in 1:k){
        beta_train <- beta_preset[-unlist(index[j]),] 
        beta_test <- beta_preset[unlist(index[j]),]
        
        #print(paste(opt,i,j))
        LM <- train_model(beta_train,model,poly,BFE)
        predictions[[j]]  <- AFM_predict(LM,model,beta_test)
        MAD.performance[[i]][j] <- mae(predictions[[j]]-beta_test$age)
        RMSE.performance[[i]][j] <- rmse(predictions[[j]]-beta_test$age)
      }
    }
    MAD <- unlist(MAD.performance)
    RMSE <- unlist(RMSE.performance)
    modelMAD[[model]] <- MAD
    if (doprint){print(paste(model,round(mean(MAD),2),ci(MAD),round(mean(RMSE),2),ci(RMSE)))}
  }
  return(modelMAD)
}
strat_wrapper <- function(beta_preset,n,k){
  #loops over performance wrapper with different stratification methods
  strat = c('random','nearest5','nearest10','age')
  stratMAD <- list()
  for (method in strat){
    stratMAD[[method]] <- performance_wrapper(beta_preset,n,k,method)}
  strat_plotter(stratMAD)
  return(stratMAD) #or perform analysis and plotting! (in this function)
}
strat_plotter <- function(strat_exp){
  #plots the result from stratwrapper
  par(mfrow=c(2,4),mar=c(7,2,2,2),oma=c(0,0,2,0))
  for (i in 1:8){
    strat_data <- unlist(strat_exp,recursive = F)[seq(i,32,8)]
    data <- data.frame(matrix(unlist(strat_data), nrow=100, byrow=F));colnames(data) <- unlist(lapply(strsplit(names(strat_data),'\\.'), function(x){x[1]}))
    setDT(data, keep.rownames = TRUE)[]
    data <- melt(data,id.vars=c("rn"))
    model_type <- unlist(strsplit(names(strat_data)[1],'\\.'))[2]
    boxplot(data$value ~ data$variable,las=2,col='red',main=model_type,cex.main=0.9)
    strat.aov <- aov(data$value ~ data$variable)
    summary(strat.aov) #aov is for balanced groups, which we have here
    print(TukeyHSD(strat.aov))
  }
  title(main=paste('Stratification Methods (in MAD)'),outer=T,cex=4);
  par(mar=c(5.1, 4.1, 4.1, 2.1),mfrow=c(1,1))
  dev.copy2pdf(file = './Images/Stratification/strat_method.pdf',width = 7, height = 7);dev.off()
}

## Outlier Functions ##
get_outliers <- function(AFM.errors,model_type){
  #calculates outliers from errors
  total_errors <- unlist(AFM.errors)
  total_errors <- split(total_errors,names(total_errors))
  mean_sample_error <- unlist(lapply(total_errors,mean))
  hist(mean_sample_error,col='red',breaks=100,main=paste(model_type,'Error Distribution'))
  outliers <- names(which(abs(mean_sample_error) > mean(mean_sample_error)+3*sd(mean_sample_error)))
  return(outliers)
}
get_outlier_wrapper <- function(beta_preset,n,k,missing,poly,BFE){
  #wrapper for get_outliers enabling iteration over model_types
  outlier_list <- list()
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  for (model_type in get_models(F)){
    print(model_type)
    AFM.errors <- simulate(beta_preset,n,k,model_type,missing,poly,BFE,bin_frac,bin_start)[[4]]
    outlier_list[[model_type]] <- get_outliers(AFM.errors,model_type)
  }
  par(mar=(c(6.2,2,2,1)))
  barplot(table(unlist(outlier_list)),las=2,col='red',main='Outlier Frequency')
  return(outlier_list)
}
plot_outlier <- function(beta_preset,outlier){
  #plot all variables with
  sample= outlier[2]
  par(oma=c(0,0,2,0))
  par(mfrow=c(3,5))
  for (i in 2:ncol(beta_preset)){
    plot(beta_preset$age,beta_preset[,i],col=ifelse(row.names(beta_preset)==sample,"red", "black"),cex=ifelse(row.names(beta_preset)==sample, 2, 1),pch=ifelse(row.names(beta_preset)==sample, 19, 1),main=colnames(beta_preset)[i],xlab='age')}
  title(main=paste('Outlier',sample),outer=T)
  path <- paste('/Images/Outliers/outlier_',sample,'.pdf',sep='');dev.copy2pdf(file = paste('.',path,sep=''),width = 7, height = 7);dev.off();par(mfrow=c(1,1),mar=c(4,4,2,1))
  print('Image saved to file!')
}
## Optimizer ##
optimize_model <- function(optimize_vector,beta_preset,n,k,...){
  #for user to optimize parameters
  total_perf <- list()
  for (opt in optimize_vector){
    if (length(list(...)) > 1){tune = expand.grid(list(...)[-1]);print(tune)}
    set.seed(1)
    MAD.performance <- list()
    RMSE.performance <- list()
    predictions <- list()
    all_predictions <- list()
    
    for (i in 1:n){
      index <- createFolds(beta_preset[,1], k = k) #create folds, NOTE that we stratify on age!!!
      MAD.performance[[i]] <- vector()
      RMSE.performance[[i]] <- vector()
      all_predictions[[i]] <- list()
      
      for (j in 1:length(index)){
        beta_train <- beta_preset[-unlist(index[j]),] 
        beta_test <- beta_preset[unlist(index[j]),]
        
        print(paste(opt,i,j))
        if (length(list(...)) > 1){LM <- train(age~.,data=beta_train,tuneGrid=tune,method = list(...)[[1]])}
        else{LM <- train(age~.,data=beta_train,method = list(...)[[1]])}
        
        
        predictions[[j]]  <- as.vector(predict(LM,newdata=beta_test,n.trees=450)); names(predictions[[j]]) <- row.names(beta_test)
        
        MAD.performance[[i]][j] <- mae(predictions[[j]]-beta_test$age)
        RMSE.performance[[i]][j] <- rmse(predictions[[j]]-beta_test$age)
      }
    }
    total_perf[[as.character(opt)]]<- unlist(MAD.performance)
    MAD <- unlist(MAD.performance); print(paste(mean(MAD),ci(MAD)))
    RMSE <- unlist(RMSE.performance); print(paste(mean(RMSE),ci(RMSE)))
  }
  
  pred_df <- data.frame(matrix(unlist(total_perf),nrow=100,byrow=F));colnames(pred_df) <- names(total_perf)
  pred_df$id <- seq(1,100,by=1)
  pred_df <- melt(pred_df, id=("id"))
  
  p <- ggplot(pred_df,aes(x=as.numeric(as.character(variable)),y=value,group=variable)) + geom_boxplot(color = 'black',position="identity",alpha=0.6) + theme_linedraw() #+ xlim(0,0.11) #alpha = transparancy
  p + labs(list(title = paste('Performance'),x='opt',y='MAD')) + theme(plot.title = element_text(size=22)) +
    stat_summary(fun.y=mean, geom="line", aes(group=1),col='red')  #+ scale_x_log10()
  beep(30)
}

## Other ##
cap_samples <- function(data,cap){
  #uniformize the data by capping the number of samples per age
  index <- c()
  set.seed(1) #make sure same samples get selected each time
  age <- data$age
  for (i in unique(age)){
    samples <- which(age == i)
    if (length(samples) <= cap){index <- c(index,samples)}
    else{index <- c(index,sample(samples,cap))}}
  index <- sort(index[!is.na(index)])
  return(index)
}

mround <- function(x,base){ 
  #bins numeric vector to closest 'base'
  return(base*round(x/base))
}

formula_poly_maker <- function(poly, beta_train){
  formula <- ''
  for (i in 1:poly){
    formula <- paste(formula,'+',paste('I(',colnames(beta_train[-1]),'^',i,')',collapse = ' + ',sep=''))
  }
  formula <- as.formula(paste('age ~',formula))
  return(formula)
}

auto <- function(vector){
  vector <- vector - mean(vector)
  vector <- vector/sd(vector)
  return(vector)
}
setup <- function(type,cores){
  #loads the necessary packages
  library(beepr2)
  library(caret)
  library(glmnet)
  library(randomForest) #randomForest
  library(e1071) #support vector regression
  library(pls) #pls and pcr
  library(earth) #Multivariate Adaptive Regression Splines
  library(entropy)
  library(MASS)

  if (type=='parallel'){
    library(foreach)
    library(doParallel)
    cl <- makeCluster(cores)
    registerDoParallel(cores=cl)
    return(cl)
  }
}

#### Residual Analysis ####
analyze_errors <- function(beta_preset,AFM.errors){
  error_by_sample <- split(unlist(AFM.errors),row.names(beta_preset))
  mean_error_by_sample <- unlist(lapply(split(unlist(AFM.errors),row.names(beta_preset)),mean))
  age <- beta_preset$age; names(age) <- row.names(beta_preset)
  plot(mean_error_by_sample~age[names(mean_error_by_sample)]);abline(h=0,lty=2,col='blue')
  plot(unlist(error_by_sample)~age[names(unlist(error_by_sample),use.names=F)]);abline(h=0,lty=2,col='blue')
}


## Performance Measures ##
rmse <- function(errors){
  #use a vector if differences, e.g. (observations - predicted)
  if (is.list(errors)){
    RMSE <- unlist(lapply(errors, function(x){sqrt(sum(x^2)/length(x))}))}
  else{
    RMSE <- sqrt(sum(errors^2)/length(errors))}
  return(RMSE)
}
mae <- function(errors){
  if (is.list(errors)){
    MAE <- unlist(lapply(errors, function(x){mean(abs(x))}))
  }
  else{
    MAE <- mean(abs(errors))
  }
  return(MAE)
}
ci <- function(values){ 
  #calculates the distance from the mean to the ends of the confidence interval
  #so mean +/- CI is the confidence interval
  #values should be a vector
  if (is.list(values)){values=unlist(values)}
  n = length(values)
  sd = sd(values)
  CI = round(qnorm(0.975)*sd/sqrt(n),2)
  return(CI)}

## NA replacers ##
NA.mean <- function(x) {
  #replaces NAs in vector x by the mean of the vector
  x[is.na(x)] <- mean(x[!is.na(x)])
  return(x)
}
NA.median <- function(x) {
  #replaces NAs in vector x by the median of the vector
  #this function is equal to the na.roughfix in the randomForest package
  x[is.na(x)] <- median(x[!is.na(x)])
  return(x)
}
NA.replace <- function(beta,mode){
  #replaces NAs in a dataframe by a mode of choice
  if (mode == 'mean'){
    beta <- apply(beta,2, NA.mean)
  }
  else if (mode == 'median'){
    beta <- apply(beta,2,NA.median)
  }
  return(as.data.frame(beta))
}
doubleCV <- function(beta){
  #feature_selection <- sort(abs(coef(CV)[,1][coef(CV)[,1]!=0]),decreasing=T)
}

## Noise Simul ##
add_noise <- function(vector,mean=0,sd,include_outliers,iter,rescale=T){
  #adds normal noise + outlier (3*sd) to the vector
  set.seed(iter)
  noise <- rnorm(length(vector),mean,sd)
  noised_vector <- vector + noise
  if (include_outliers){outliers <- as.logical(rbinom(length(vector),1,0.1)) #select percentage of data
  outlier_noise <- rnorm(sum(outliers),0,3*sd)
  noised_vector[outliers] <- noised_vector[outliers] + outlier_noise}
  if (rescale){noised_vector <- (noised_vector - mean(noised_vector)) * 0.4/(max(noised_vector) - min(noised_vector)) + 0.64}
  return(noised_vector)
}
parsim <- function(data,max,inc,cols,include_outliers,poly,exclude_model,j,doprint=F){
  results <- list()
  
  for (i in seq(0,max,inc)){
    #print('');print(i);print('')
    fake <- add_noise(data$age,0,i,include_outliers,j)
    #plot(data$age,fake,col='blue')
    fake_data <- as.data.frame(cbind(data[,c(1,cols)],fake));colnames(fake_data)[1] <- 'age'
    out <- performance_wrapper(fake_data,n = 10,k=10,strat = 'normal',poly = 5,BFE=F,exclude_model,doprint=F)
    results[[(length(results)+1)]] <- unlist(lapply(out,mean))
  }
  df <- as.data.frame(do.call(rbind, lapply(results, c)));df$SD <- seq(0,max,inc)
  return(df)
}
test <- function(data,j){
  data <- data$garbage1 * j 
  data <- as.data.frame(data)
  return(data)
}
add_df <- function(a,b){
  return(a+b)
}
valsim <- function(data,max=10,inc=1,cols=NULL,include_outliers = T,poly=5,exclude_model='RF',iter=1,doprint=F,title){
  #simulate with fake data, cols argument is to add additional features from original data
  
  print("Parallel Simulation, please wait...")
  if (is.null(cols) & is.null(exclude_model)){stop('for univariate predictor please exclude RF model')}
  df <- foreach(j=1:iter,.combine=add_df) %dopar% {parsim(data,max,inc,cols,include_outliers,poly,exclude_model,j,doprint=F)}
  df <- df/iter
  colnames(df)[2:3] <- paste(colnames(df)[2:3],poly)
  
  #create title and file name
  name <- title
  #if (!is.null(cols)){name <- paste(c(name,colnames(data)[cols]),collapse = '.');title <- paste(c(title,colnames(data)[cols]),collapse=' +')}
  #if (include_outliers){name <- paste(name,'_out',sep='');title <- paste(title,' + Outliers',sep='')}
  #name <- paste(name,'_n',iter,sep='');title <-  paste(title,' n = ',iter,sep='')
  
  #calculate MSE of fit on fake ELOVL2
  age <- data$age
  MSE <- list()
  for (i in seq(0,max,inc)){
    MSE[[i+1]] <- vector()
    for (j in 1:iter){
      set.seed(iter)
      fake <- add_noise(data$age,0,i,include_outliers,j)
      MSE[[i+1]][j] <- mean(lm(age~fake,data=as.data.frame(cbind(age,fake)))$residuals^2)
    }
    
  }
  MSE <- unlist(lapply(MSE,mean))
  df$MSE <- MSE
  
  #Plot results
  mycolors <- rainbow(8); mycolors[4] <- 'black';names(mycolors) <- get_models(F)
  mycolors <-  mycolors[!names(mycolors) %in% exclude_model]
  for (i in 1:(length(colnames(df))-1)){
    if (i==1){plot(df$MSE,df[,i],type = 'l',col=mycolors[i],xlab='MSE of lm(age ~ fake+noise)',ylab='MAD',main=title,ylim=c(0,8),xlim=c(0,60))}#'Fake ELOVL2 + Rest'
    else{lines(df$MSE,df[,i],col=mycolors[i])}}
  abline(v=mean(lm(age~cg16867657,data=data)$residuals^2),col='red',lty=2)
  legend('bottomright',legend = head(colnames(df),-2),cex=0.6,fill=mycolors)
  dev.copy2pdf(file = paste('./Images/ValSim/',name,'.pdf',sep=''),width = 12, height = 9);dev.off()
  beep(28)
}


########################
## Feature Subsetting ##
########################

CV <- function(beta,k=10){
  
  performance <- list()
  
  for (i in 1:10){
    set.seed(i)
    index <- createFolds(beta[,1], k = k)
    performance[[i]] <- vector()
    features[[i]] <- list()
    
    for (j in 1:k){
      train <- beta[-unlist(index[j]),] 
      test <- beta[unlist(index[j]),]
      print(paste(i,j))
      
      model <- lm(age~.,train)
      performance[[i]][j] <- mae(test$age - predict(model,test))
    }
  }
  return(performance)
}
mutual_entropy <- function(x){
  #age vector is required and CpGsite as input
  y2d <- discretize2d(x,age,numBins1=20,numBins2=20)
  MI = entropy(rowSums(y2d)) + entropy(colSums(y2d)) - entropy(y2d)
  return(MI)
}

FPS <- function(beta,fs_choice,threshold){
  #fs_methods = c('sd','cor','r^2','regression','mutual_entropy')
  age <- beta[,1]
  
  if (fs_choice == 'sd'){
    print('Feature selection by standard deviation')
    CpG_sd <- apply(beta[,-1],2,sd)
    hist(CpG_sd,breaks=100,col='red')
    CpG_subset <- subset(CpG_sd,CpG_sd > threshold) #0.02
    
  }else if(fs_choice == 'cor'){
    print('Feature selection by correlation')
    CpG_correlation <- apply(beta[,-1],2, function(x) cor(beta$age,x))
    hist(abs(CpG_correlation),breaks=100,col='red')
    CpG_subset <- subset(CpG_correlation,abs(CpG_correlation) > threshold) #0.1
    
  }else if(fs_choice == 'r^2'){
    clusterExport(cl, "age") #needed for cluster
    print('Feature selection by correlation')
    CpG_R_squared <- unlist(parLapply(cl,beta[,-1], function(x) {
      summary(lm(x ~ age))$adj.r.squared}))
    hist(abs(CpG_R_squared),col='red',breaks=100)
    CpG_subset <- subset(CpG_R_squared,abs(CpG_R_squared) > threshold) #0.01
    
  }else if(fs_choice == 'regression'){
    print('Feature selection by "regression coefficients"')
    clusterExport(cl, "age") #needed for cluster
    CpG_regression_coefficients <- unlist(parLapply(cl,beta[,-1], function(x){lm(x ~ age)$coefficients[2]}))
    hist(abs(CpG_regression_coefficients),col='red',breaks=100)
    CpG_subset <- subset(CpG_regression_coefficients,abs(CpG_regression_coefficients) > threshold)#0.0002
  
  }else if(fs_choice == 'poly'){
    print('Standardizing')
    beta[,-1] <- as.data.frame(parLapply(beta[,-1],auto))
    print('Feature selection by MSE of polynomial fit')
    attach(beta)
    CpG_MSE <- foreach(i = 1:ncol(beta[,-1]), .combine=c) %dopar% { mean(lm(as.formula(paste('age ~ poly(',colnames(beta[i+1]),',5)')))$residuals^2)}
    hist(CpG_MSE,col='red',breaks=100)
    CpG_subset <- subset(CpG_MSE,abs(CpG_MSE) > threshold)#0.0002
    
  }else if(fs_choice=='mutual_entropy'){
    clusterExport(cl,list('age','entropy','discretize2d'))
    CpG_ME <- unlist(parLapply(cl,beta,mutual_entropy))
    names(CpG_ME) <- row.names(beta)
    hist(CpG_ME,breaks=100,col='red')
    CpG_subset <- subset(CpG_ME,CpG_ME > threshold)
  
  }else if(fs_choice=='fcbf'){
    #fact correlation based filtering
    library(Biocomb)  
    z <- select.fast.filter(matrix=as.matrix(beta[,c(2:ncol(beta),1)]),disc.method='MDL',threshold=0.5,attrs.nominal=numeric())
  
  }else if(fs_choice=='markov_blanket'){
    library(bnlearn)
    z <- fast.iamb(x=beta[,2:100],cluster=cl,optimized=T)
  
  }else if(fs_choice=='mRMR'){
    #minimum redundancy maximum relevance
    print('Feature selection by Minumum Redundance Maximum Relevance')
    library(mRMRe)
    beta$age <- as.numeric(beta$age)
    fs <- new("mRMRe.Filter", data = mRMR.data(beta), target_indices = 1,levels = c(8,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
    CpG_subset <- colnames(beta)[solutions(fs)[[1]][,1]]
    
    }
  beep(28)
  return(CpG_subset)
}

get_lasso_coefs <- function(cv.glmnet){
  #returns features from cv.glmnet object
  coefs <- coef(cv.glmnet,s=cv.glmnet$lambda.1se) 
  features <- sort(coefs@Dimnames[[1]][coefs@i+1][-1]) #remove intercept
  return(features)
}

add_fake_and_antifake <- function(sub){
  #adds 2 new variables to data frame which both add up to the first column (which should be age)
  set.seed(1)
  #sub <- beta[,1:1000]
  a = 0.9
  fake_var <- rnorm(nrow(sub),1,1);fake_var <- fake_var + abs(min(fake_var));fake_var <- fake_var/max(fake_var)
  new_var <- (sub$age/100 - a*fake_var)/(1-a);new_var <- new_var+abs(min(new_var));new_var <- new_var/max(new_var)
  plot(sub$age,new_var,col='red',main='f2',xlab='Age',ylab='methylation fraction')
  plot(sub$age,fake_var,col='red',main='f1',xlab='Age',ylab='methylation fraction')
  
  sub <- cbind(sub,cbind(fake_var,new_var))
  model <- lm(age~.,data=sub[,c(1,ncol(sub)-1,ncol(sub))])
  #plot(sub$age, model$coefficients[3] *new_var + model$coefficients[2]*fake_var + model$coefficients[1], ylab='predicted age',xlab='real age',col='blue')
  return(sub)
}

single_wrapper <- function(beta_preset,n,k,strat,poly,BFE,model,doprint=F){
  #gives mean MAD/RMSE and confidence interval for all models
  
  #stratification methods
  if (strat=='random'){ages <- 1:nrow(beta_preset);if (doprint){print('Random Folds')}}
  else if (strat=='nearest5'){ages <- mround(beta_preset$age,5);if (doprint){print('Stratifying by closest 5 Folds')}}
  else if (strat=='nearest10'){ages <- mround(beta_preset$age,10);if (doprint){print('Stratifying by closest 10 Folds')}}
  else{ages <- beta_preset[,1];if (doprint){print('Stratifying by Age Folds')}} #normal mode;
  
  modelMAD <- list()
  
  set.seed(1)
  MAD.performance <- list()
  RMSE.performance <- list()
  predictions <- list()
  
  for (i in 1:n){
    index <- createFolds(ages, k = k) #create folds
    MAD.performance[[i]] <- vector()
    RMSE.performance[[i]] <- vector()
    
    for (j in 1:k){
      beta_train <- beta_preset[-unlist(index[j]),] 
      beta_test <- beta_preset[unlist(index[j]),]
      
      #print(paste(opt,i,j))
      LM <- train_model(beta_train,model,poly,BFE)
      predictions[[j]]  <- AFM_predict(LM,model,beta_test)
      MAD.performance[[i]][j] <- mae(predictions[[j]]-beta_test$age)
      RMSE.performance[[i]][j] <- rmse(predictions[[j]]-beta_test$age)
    }
  }
  MAD <- unlist(MAD.performance)
  RMSE <- unlist(RMSE.performance)
  modelMAD[[model]] <- MAD
  if (doprint){print(paste(model,round(mean(MAD),2),ci(MAD),round(mean(RMSE),2),ci(RMSE)))}
  return(modelMAD)
}
