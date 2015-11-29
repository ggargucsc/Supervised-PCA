library(spls)   #dataset  


#data is already standarized
data(prostate)

y <- prostate$y
x <- prostate$x
prostateData <- as.data.frame(cbind(y,x))
colnames(prostateData)[1] <- "y"



#function to reduce dimensions of original data using PCA 
#returns compact data (PCs Scores)
pcaData <- function(data, k)
{
  
  #number of principal components
  n <- min(nrow(data), ncol(data))
  pca <- prcomp(data)
  lambda <- pca$sdev^2
  lambdaProp <- lambda/sum(lambda)
  numPCs <- length(which(cumsum(lambdaProp)<k)) #k
  
  compactData <- pca$x[, 1:numPCs] 
  
  
  return(compactData)
}

#function to implement supervised pca, returns reduced data matrix, and then we apply pca
spcaData <- function(data, m)
{
  #data <- prostateData
  mylogit <- NULL
  n <- ncol(data)-1
  
  for(i in 1:n)
  {
    temp <- coef(summary(glm(data$y ~ data[,(i+1)], family = "binomial")))[2]
    mylogit <- c(mylogit,temp)
  }
  
 
  index <- which(abs(mylogit)>m)
  reducedX <- data[, c(index+1)]
  
  return(reducedX)
}


#classification using Linear Discriminant Analysis
ldaData <- function(dataTrain, dataTest)
{
  n <- nrow(dataTrain)
 
  #prior probabilities for each of the two classes
  prior <- table(dataTrain[,1])/n
  
 
  #dividing dataset based on classes 
  X1 <- dataTrain[dataTrain[,1]==0,][,-1]
  X2 <- dataTrain[dataTrain[,1]==1,][,-1]
  
  meanX1 <- as.matrix(apply(X1, 2, mean))
  meanX2 <- as.matrix(apply(X2, 2, mean))
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  
  #estimating sigma 
  sigmaHat <- ((n1-1)*var(X1)+(n2-1)*var(X2))/n
  invSigHat <- ginv(sigmaHat)
  
  n.test <- nrow(dataTest)
  LDAscores <- array(NA, dim=c(n.test,2))
  
  #LDA scores
  for(i in 1:n.test)
  {
    x <- t(as.matrix(dataTest[i,-1]))
    LDAscores[i, 1] <- 1/2*t(x-meanX1)%*%invSigHat%*%(x-meanX1)-log(prior[1])
    LDAscores[i, 2] <- 1/2*t(x-meanX2)%*%invSigHat%*%(x-meanX2)-log(prior[2])
  }
  
  LDApred <- apply(LDAscores, 1, which.min)
  LDApred <- as.factor(LDApred)
  
  tableFinal <- table(dataTest[,1],LDApred)
  
  
  return(tableFinal)  
}

#classification using Quadratic Discriminant Analysis
qdaData <- function(dataTrain, dataTest)
{
  n <- nrow(dataTrain)
  
  #prior probabilities for each of the two classes
  prior <- table(dataTrain[,1])/n
 
  X1 <- dataTrain[dataTrain[,1]==0,][,-1]
  X2 <- dataTrain[dataTrain[,1]==1,][,-1]
  
  meanX1 <- as.matrix(apply(X1, 2, mean))
  meanX2 <- as.matrix(apply(X2, 2, mean))
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  
  #estimating sigma 
  sig1Hat <- var(X1)*(n1-1)/n1
  #print(sig1Hat)
  invSig1Hat <- ginv(sig1Hat)
  
  sig2Hat <- var(X2)*(n2-1)/n2
  invSig2Hat <- ginv(sig2Hat)
  
  
  n.test <- nrow(dataTest)
  LDAscores <- array(NA, dim=c(n.test,2))
  #LDA scores
  for(i in 1:n.test)
  {
    x <- t(as.matrix(dataTest[i,-1]))
    LDAscores[i, 1] <- 1/2*t(x-meanX1)%*%invSig1Hat%*%(x-meanX1)-log(prior[1])
    LDAscores[i, 2] <- 1/2*t(x-meanX2)%*%invSig2Hat%*%(x-meanX2)-log(prior[2])
  }
  
  LDApred <- apply(LDAscores, 1, which.min)
  LDApred <- as.factor(LDApred)
  #finalCheck <- data.frame(cbind(dataTest[,1], LDApred))
  tableFinal <- table(dataTest[,1],LDApred)
  
  
  return(tableFinal)  
}


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
