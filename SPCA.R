library(spls)   #dataset  


#load the functions file
source('~/Documents/ggargucsc_github/Supervised-PCA/functions.r')

#data is already standarized
data(prostate)

y <- prostate$y
x <- prostate$x
prostateData <- as.data.frame(cbind(y,x))
colnames(prostateData)[1] <- "y"



################################ main body  #####################################

############# apply pca on dataset and then do discriminant analysis #####
seqVar <- seq(from=0.70, to=0.99, by=0.04)
rep <- 100
errlda=dim(rep)
errqda=dim(rep)
merrlda=dim(length(seqVar))
merrqda=dim(length(seqVar))
n <- nrow(prostateData)
ntrain <- 71          #0.70*n
ntest <- 31          #0.30*n


for(i in 1:length(seqVar))
{
  #apply pca on training dataset(only X's) to get compact data
  compactData <- pcaData(prostateData[,-1], seqVar[i]) 
  finalCompactData <- as.data.frame(cbind(prostateData[,"y"],compactData))
  colnames(finalCompactData)[1] <- "y"
  
  for (k in 1:rep) 
  {
    #sample ntrain samples from n 
    train <- sample(1:n, ntrain)
    
    #linear discriminant analysis on compact data
    library(MASS)
    mLda=ldaData(finalCompactData[train,], finalCompactData[-train,])
    errlda[k]=(ntest-sum(diag(mLda)))/ntest
    
    #quadratic discriminant analysis on compact data
    mQda=qdaData(finalCompactData[train,], finalCompactData[-train,])
    errqda[k]=(ntest-sum(diag(mQda)))/ntest
    
  }
  (merrlda[i]=mean(errlda))
  (merrqda[i]=mean(errqda))
}
seqPCs <- c(7,12,19,29,43, 63)
errorRate <- as.data.frame(cbind(seqVar,merrlda,merrqda))
colnames(errorRate) <- c("x", "y1", "y2")
plot(y1 ~ x, data=errorRate, type="b", col="red", ylim=range(0,1),xlim=range(0.70,0.98),xlab="proportion of variance (proportional to number of PCs)", ylab="Test Error")
points(y2 ~ x, errorRate,pch = 2, type="b", col="green")
legend(0.85, 0.8,c("PCA-LDA", "PCA-QDA"), lty=c(1,1),lwd=c(1.5,1.5), col=c("red", "green"), bty="n",bg ="white")

###############################################################################
########## apply spca on dataset, choosing different threshold values #########

seqThreshold <- seq(from=0, to=4.2, by=0.3)
seqVar <- seq(from=0.70, to=0.99, by=0.04)
rep <- 100
errlda=dim(rep)
errqda=dim(rep)
merrlda=dim(length(seqVar))
merrqda=dim(length(seqVar))
n <- nrow(prostateData)
ntrain <- 71          #0.70*n
ntest <- 31          #0.30*n

for(i in 1:length(seqVar))
{
  reducedData <- spcaData(prostateData,3.9) 
  scompactData <- as.data.frame(cbind(prostateData[,"y"], pcaData(reducedData, seqVar[i])))
  #finalCompactData <- as.data.frame(cbind(prostateData[,"y"],compactData))
  colnames(scompactData)[1] <- "y"
  
  for (k in 1:rep) 
  {
    #sample ntrain samples from n 
    train <- sample(1:n, ntrain)
    
    #linear discriminant analysis on compact data
    library(MASS)
    mLda=ldaData(scompactData[train,], scompactData[-train,])
    errlda[k]=(ntest-sum(diag(mLda)))/ntest
    
    #quadratic discriminant analysis on compact data
    mQda=qdaData(scompactData[train,], scompactData[-train,])
    errqda[k]=(ntest-sum(diag(mQda)))/ntest
    #table(scompactData$y[-train],predict(m1,scompactData[-train,])$class)
    
  }
  (merrlda[i]=mean(errlda))
  (merrqda[i]=mean(errqda))
}
errorRate <- as.data.frame(cbind(seqVar,merrlda,merrqda))
colnames(errorRate) <- c("x", "y1", "y2")
plot(y1 ~ x, data=errorRate, type="b", col="red", ylim=range(0,1),xlab="proportion of variance (proportional to number of PCs)", ylab="Test Error")
points(y2 ~ x, errorRate,pch = 2, type="b", col="green")
legend(0.87, 0.8,c("SPCA-LDA", "SPCA-QDA"), lty=c(1,1),lwd=c(2.5,2.5), col=c("red", "green"), bty="n",bg ="white")

######################################################################################
