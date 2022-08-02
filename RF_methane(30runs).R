#RF_methane(30runs).R

rm(list=ls())
library(caret)
path <- "~//dir//"
setwd(path)

###### Grab data
myData<-read.csv(paste0(path,"methanenew.csv")  ,header=T,sep=",",quote=NULL, comment="")
str(myData)
colnames(myData)

numericData<-myData[,-c(1:7)]
colnames(numericData)

# highCH4group <-numericData [numericData$CH4>=2000,]
discreteDat <- data.frame(myData$LocationAbbrev,
                          myData$Location, 
                          myData$LithotypeAbbrev, 
                          myData$GenericLithotype, 
                          myData$Lithotype, 
                          myData$Sample, 
                          myData$Ophiolite)          

# fit1 <- lm(CH4 ~., data = numericData)
# # plot(numericData$CH4..ppmv ~ numericData$Rh..0.2.ppb_NI.FINA., data = numericData)
# summary(fit1)
# pd <- predict(fit1, numericData, interval = "confidence")
# str(pd)
# plot(numericData$CH4,pd[,1])
# abline( lm(numericData$CH4~as.vector(pd[,1])))
# 
# ############################
# # caret: lm
# #################################
# # In R, the annual_pm~. means use annual_pm as the model response
# # and use all other variables as predictors
# lm1 <- train(CH4 ~., data = numericData, method = "lm")
# lm1$finalModel
# # Unrealistic R-squared based on final model
# summary(lm1$finalModel)$r.squared
# # Realistic R-squared based on resampling
# lm1
# # This table has 1 row
# lm1$results
# 
# ######################
# # lm + CV
# 
# # Set up a 10-fold cross validation
# tc <- trainControl(method = "boot_all")
# tc <- trainControl(method = "cv", number = 10)
# # Include the setup in your model
# lm1_cv <- train(CH4 ~., data = numericData, method = "lm",
#                 trControl = tc) # here
# lm1_cv
# 
# 
# ############################
# # caret: rf
# #################################
# rf1 <- train(CH4 ~., data = numericData, method = "rf")
# rf1$finalModel
# # This table has 3 rows, we explain later
# rf1$results
# 
# tg <- data.frame(mtry = seq(2, 10, by =2))
# rf1 <- train(CH4 ~., data = numericData,  method = "rf", tuneGrid = tg)
# rf1$results
# 
# ##############################
# model_list <- list(lm = lm1, rf = rf1)
# res <- resamples(model_list)
# summary(res)
# 
# ###############################################
# ggplot(numericData, aes(CH4, Rh..0.2.ppb_NI.FINA.)) + geom_point() + geom_smooth()
# 
# #########################################
# samp <- createDataPartition(numericData$CH4, p = 0.8, list = FALSE)
# training <- numericData[samp,]
# testing <- numericData[-samp,]
# 
# ############################################
# # Set seed for reproducibility
# set.seed(100)
# 
# # Set up the resampling, here repeated CV
# tr <- trainControl(method = "repeatedcv", number = 3, repeats = 5)
# 


#cor(basinDat[,-(1:2)])

# check that we get mean of 0 and sd of 1
#colMeans(scaled.wl)  # faster version of apply(scaled.dat, 2, mean)
#apply(scaled.wl, 2, sd)
set.seed(998)

modellist <- list()
getTrainPerf <- list()
varImp <- list()
rank <-NULL
finalModelPredicted <- list()
trainVSmodel <-list()
testVSmodel <-list()
trainPostResample <-list()
testPostResample <-list()
sumrank <- rep(0, (dim(numericData)[2]-1))
i<-1
while(i<=30){
  ### Using Well-log data###
  #set.seed(998)
  print(paste0("i = ", i))
  inTraining <- createDataPartition(numericData$CH4, p = .80, list = FALSE)
  training <- numericData[ inTraining,]
  testing  <- numericData[-inTraining,]
  
  #trainControl can be used to specifiy the type of resampling
  fitControl1 <- trainControl(## 10-fold CV
                             method = "repeatedcv",
                             number = 10,
                             ## repeated ten times
                             repeats = 300)
  
  fitControl2  <- trainControl(method = "cv", 
  			     		number = 10, 
                         	returnResamp = "all")
  
  fitControl3 <- trainControl(method = "LOOCV")
  fitControl4 <- trainControl(method = "boot_all")
  
  fitControlR  <- trainControl(method = "cv", number = 10, returnResamp = "all", search = "random")
  
  
  ### first fit############################
  #set.seed(998)
  # nnetGrid <- expand.grid(mtry = seq(1,10,1), ntree = c(200,500,700, 1000,2000) )
  nnetGrid <- expand.grid(mtry = seq(1,2,1) )
  
  
  # for (ntree in c(1000,1500,2000,2500)){
  ntree <- 600
  nnetFit <- train(CH4 ~., 
                    data = training,
                    method = "rf",
  		              #preProcess = c("center","scale"),
  		              # preProcess = "scale",
                    trControl = fitControl1,
                    tuneGrid = nnetGrid,
                    ntree = ntree,
                    verbose = FALSE, trace = F, linout = 1)
  # key <- toString(ntree)
  key <- toString(i)
  
  modellist[[key]] <- nnetFit
  getTrainPerf[[key]] <- getTrainPerf(nnetFit)
  trainRSq <- getTrainPerf(nnetFit)[2]
  varImp[[key]] <- cbind(varImp(nnetFit)$importance,rank(-varImp(nnetFit)$importance$Overall))
  rank <- cbind(rank, rank(-varImp(nnetFit)$importance$Overall))
  sumrank <- sumrank + rank(-varImp(nnetFit)$importance$Overall)
  #==========================================================================================
  png(filename=paste0("~//dir//","//RF_reg//0",i,"_errorModels.png"))
  par(xpd = FALSE) # Avoid clipping the text in some device
  plot(nnetFit) #plot the rmse of all models
  dev.off()
  #==========================================================================================
  # final value evaluation
  cData<-cbind(nnetFit$finalModel$predicted,training$CH4)
  finalModelPredicted[[key]] <- cData
  fit <- lm(nnetFit$finalModel$predicted ~ training$CH4, data = data.frame(cData))
  write(toString(summary(fit)), file=paste0("~//dir//","//RF_reg//0",i,"_finalModelFitSummary.txt"))
  png(filename=paste0("~//dir//","//RF_reg//0",i,"_finalmodel.png"))
  par(xpd = FALSE) # Avoid clipping the text in some device
  plot(nnetFit$finalModel$predicted ~ training$CH4, data = data.frame(cData), 
       xlim = c(0,max(c(nnetFit$finalModel$predicted,training$CH4))), 
       ylim= c(0,max(c(nnetFit$finalModel$predicted,training$CH4))), 
       xlab = "actual CH4", ylab ="predicted training CH4")
  abline(fit, col=2)# R2 <- 1 - (sum((training$CH4-nnetFit$finalModel$predicted)^2)/
  abline(fit, col=2)
  abline(coef = c(0,1), lty=2)
  dev.off()
  #==========================================================================================
  # predicting training data
  trainps <- predict(nnetFit , training)
  trfit <- lm(trainps ~ training$CH4)
  write(file=paste0("~//dir//","//RF_reg//0",i,"_finalModelFitSummary.txt"), toString(summary(trfit)))
  trainpsData<-cbind(trainps, training$CH4)
  trainVSmodel[[key]]<- trainpsData
  trainPostResample[[key]] <- round(postResample(pred=trainps, obs=training$CH4),3)
  png(filename=paste0("~//dir//","//RF_reg//0",i,"_trainedModel.png"))
  par(xpd = FALSE) # Avoid clipping the text in some device
  plot(trainps ~ training$CH4, xlim = c(0,max(c(trainps,training$CH4))), 
       ylim= c(0,max(c(trainps,training$CH4))), 
       xlab = "actual training CH4", ylab ="trained CH4")
  abline(trfit, col=2)
  abline(coef = c(0,1), lty=2)
  dev.off()
  #==========================================================================================
  # predicting test data
  ps <- predict(nnetFit , testing)
  ###Examine results ps ###########
  testfit <- lm(ps ~ testing$CH4, data = data.frame(cbind(ps,testing$CH4)))
  summary(testfit)
  write(file=paste0("~//dir//","//RF_reg//0",i,"_testFitSummary.txt"), toString(summary(testfit)))
  testpredobs<- cbind(ps,testing$CH4)
  testVSmodel[[key]]<- testpredobs
  testPostResample[[key]] <- round(postResample(pred=ps, obs=testing$CH4),3)
  testRSq <- round(postResample(pred=ps, obs=testing$CH4),3)[2]
  
  png(filename=paste0("~//dir//","//RF_reg//0",i,"_testingModel.png"))
  par(xpd = FALSE) # Avoid clipping the text in some device
  plot(ps ~ testing$CH4, data = cbind(ps,testing$CH4),
       xlim=c(0,max(ps,testing$CH4)),ylim=c(0,max(ps,testing$CH4)),
       xlab = "actual testing CH4", ylab ="predicted CH4")
  abline(testfit, col=2)
  abline(coef = c(0,1), lty=2)
  dev.off()
  
  if(trainRSq>testRSq) i <- i+1
}

#Compare results
results <- resamples(modellist)
summary(results)
write(file=paste0("~//dir//","//RF_reg//resampleModellistSummary.txt"), toString(modellist))

#==========================================================================================
write(file="~//dir//RF_reg//getTrainPerf.txt", toString(getTrainPerf))
write(file="~//dir//RF_reg//varImp.txt", toString(varImp))
write(file="~//dir//RF_reg//finalModelPredicted.txt", toString(finalModelPredicted))
write(file="~//dir//RF_reg//trainVSmodel.txt", toString(trainVSmodel))
write(file="~//dir//RF_reg//testVSmodel.txt", toString(testVSmodel))
write(file="~//dir//RF_reg//trainPostResample.txt", toString(trainPostResample))
write(file="~//dir//RF_reg//testPostResample.txt", toString(testPostResample))
write(file="~//dir//RF_reg//sumRank.txt", 
      cbind(row.names(varImp(nnetFit)$importance),sumrank))
write.table(file="~//dir//RF_reg//rank.csv", rank, sep = "\t")

