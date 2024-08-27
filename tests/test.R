setwd("P:/VanGinkelKarch3/")
library(mice)
library(MASS)
source("RsquareMI.r")
#Data properties
mu = c(2,5,10,20)
sigma = matrix(c(5,0,0,3,0,5,1,0,0,1,5,0,3,0,0,10),nrow = 4, ncol = 4)
sigmaerror = 0.6
muerror = 0
NumberOfRsquareds = 2
Populationbs = c(0.2,0.1,0.1,-0.1,0)
#Populationbs = c(0,0,0,0,0)
n=500
percent=0.25
NumberOfImp = 25

#Simulating data
sigmaYhat = 0
ii = 1
for (i in 1:nrow(sigma)){
  ii = ii + 1
  jj = 1 
  for (j in 1:ncol(sigma)){
    jj = jj + 1
    sigmaYhat = sigmaYhat + Populationbs[ii]*Populationbs[jj]*sigma[i,j]
  }
}
  
sigmaY <- sigmaYhat + sigmaerror
#Seed = 234 
#set.seed(Seed)
X = mvrnorm(n = n, mu=mu, Sigma=sigma)        
error = mvrnorm(n=n,mu=0, Sigma=sigmaerror)
colnames(X) = c("X1","X2","X3","X4")
intercept = rep(1,times=n)
interceptX = cbind(intercept, X)
Y = rep(0,times=n)
for (j in 1:length(Populationbs)) {Y = Y + Populationbs[j]*interceptX[,j]}
Y = Y + error
colnames(Y) = "Y"
XY = cbind(X,Y)     
   

#Simulating missing data
Numberofmissing = percent*(nrow(XY)*ncol(XY))
#Seed = 235 
#set.seed(Seed) 
randommatrix = runif((nrow(XY)*(ncol(XY)-1))) 
randommatrix = matrix(randommatrix,nrow(XY),(ncol(XY)-1))
missingmatrix = matrix(1,nrow(XY),(ncol(XY)-1)) 
weightmatrix <- matrix(1, nrow(XY), (ncol(XY)-1)) 
randommatrix <- randommatrix*weightmatrix
stoploop = FALSE
for (k in 1:Numberofmissing){
  for (ii in 1:nrow(randommatrix)){
    for (jj in 1:ncol(randommatrix)){
      if (randommatrix[ii,jj] == max(randommatrix) & any(randommatrix[ii,-jj] != 0)) {
        missingmatrix[ii,jj] = 0
        randommatrix[ii,jj] = 0 
        stoploop = TRUE
        break             
      }
    }
    if (stoploop == TRUE) {
      stoploop = FALSE
      break
    }
  }
}
XYmis <- XY
XYmis[,2:ncol(XYmis)][missingmatrix==0] = NA
seed <- 554
 
XYimp  <- mice(XYmis, maxit = 10, m = NumberOfImp, printFlag = FALSE, seed = seed)

model = as.formula("Y~X1+X2+X3+X4")

results = RsquareSP(XYimp, model, cor=TRUE, beta=TRUE) 

results

        