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
