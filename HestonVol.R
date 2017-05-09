#This code generates a time series following the Heston Model of pre-defined length, initial stock and vol value.
#It also generates a White Noise of the same length, and a pre-defined standard deviation.
#It adds these together, generates a periodogram, then runs a Newton Method on an MLE to generate the ratios we shrink this by.
#After, it multiplies the periodogram by the given ratio. It runs for a pre-defined number of times, and saves all relevant data into
#a matrix. All histograms and calculations that are in the project write-up were done in the console after running.

#To give the log of a time series that follows the Heston Model. It has length "len", and the exponential as initial stock and vol StockStart, VolStart.
Heston.Series <- function(len, StockStart, VolStart)
{
  k <- 5
  a <- 0.04
  gamma <- 0.5
  mu <- 0.05
  correlation <- -0.5
  stock <- StockStart
  logstock <- log(stock)
  vol <- VolStart
  loglist <- numeric(len)
  loglist[1] <- logstock
  dt <- 1/(len*252)
  for(t in (2:len))
  {
    UnCorrVar <- numeric(2)
    UnCorrVar[1] <- rnorm(1)
    UnCorrVar[2] <- rnorm(1)
    CorrVar <- numeric(2)
    CorrVar <- Correlate.Variables(UnCorrVar, correlation)
    loglist[t] <- loglist[t-1] + mu*dt - 0.5*vol*dt + sqrt(vol*dt)*CorrVar[1]
    vol <- vol + k*(a-vol)*dt + gamma*sqrt(vol*dt)*CorrVar[2]
    if(vol<0){
      cat(vol)
    }
  }
  return(loglist)
  
}

#This function takes in a vector of two uncorrelated standard normal variables, and correlates them with correlation coefficient "correlation"
#using a Cholesky Decomposition
Correlate.Variables <- function(UnCorrVar, correlation)
{
  NewVar <- numeric(2)
  NewVar[1] <- UnCorrVar[1]
  NewVar[2] <- correlation*UnCorrVar[1] + sqrt(1-correlation^2)*UnCorrVar[2]
  return(NewVar)
}

#Generates a white noise of length "len" with standard deviation "stdev"
White.Noise <- function(len, stdev)
{
  BList <- numeric(len)
  BList <- rnorm(len, sd = stdev)
  return(BList)
}

#This is only here to make vectorising later easy, to avoid the use of for loops.
fk <-vector(length = len)
for(p in 1:len){
  fk[p] = p
}

#A series of coefficients that appear a lot throughout the project. We define them here.
sinsummands <-vector(length = len)
for(p in 1:len){
  sinsummands[p] = 2*sin(pi*fk[p]*delt)
  sinsummands[p] = (Mod(sinsummands[p]))^2
}

#This function is the main Whittle Likelihood, what is being maximised in the MLE. In these and the methods below, we detail how to find the
#gradient and Hessian of this function, both needed to run the Newton algorithm. Each derivative is indexed before the start of the method.
#All functions are reparameterised to the whole real line, as described in the write-up.
whittlelik <- function(sigsq)
{
  idx <- 1:((len+1)/2-1)
  likelihood <- vector(length = (len/2)-1)
  p <- 1:(((len+1)/2)-1)
  likelihood[p] = -log(sigsq[1] + sigsq[2]*sinsummands[p]) - modulus[(p+floor(len/2))]/(sigsq[1] + sigsq[2]*sinsummands[p])
  whittle <- sum(likelihood) + sigsq[1] + sigsq[2]
  return(whittle)
}
#dwhittle/dsigx^2
XGrad <- function(sigsq)
{
  likelihood <- vector(length = (len/2)-1)
  p <- 1:((len+1)/2-1)
  #for(p in 1:(((len+1)/2)-1)){
  likelihood[p] = -1/(sigsq[1] + sigsq[2]*sinsummands[p]) + modulus[(p+floor(len/2))]/((sigsq[1] + sigsq[2]*sinsummands[p]))^2
  #}
  xgrad <- sum(likelihood)
  xgrad <- xgrad * sigsq[1] + 1
  return(xgrad)
}
#dwhittle/dsige^2
EGrad <- function(sigsq)
{
  likelihood <- vector(length = (len/2)-1)
  p <- 1:((len+1)/2-1)
  #for(p in (1:((len+1)/2)-1)){
  likelihood[p] = -sinsummands[p]/(sigsq[1] + sigsq[2]*sinsummands[p]) + modulus[(p+floor(len/2))]*sinsummands[p]/(sigsq[1] + sigsq[2]*sinsummands[p])^2
  #}
  xgrad <- sum(likelihood)
  xgrad <- xgrad * sigsq[2] + 1
  return(xgrad)
}
#Forming the Grad vector
grad <- function(sigsq)
{
  gradient <- matrix(nrow = 2, ncol = 1)
  gradient[1] <- XGrad(sigsq)
  gradient[2] <- EGrad(sigsq)
  #cat("Gradient: ")
  #print(gradient)
  return(gradient)
}
#d2whittle/d(sigx^2)^2
HessianDoubleX <- function(sigsq)
{
  likelihood <- vector(length = (len/2)-1)
  p <- 1:((len+1)/2-1)
  likelihood[p] = 1/(sigsq[1] + sigsq[2]*sinsummands[p])^2 - 2*modulus[(p+floor(len/2))]/(sigsq[1] + sigsq[2]*sinsummands[p])^3
  DoubleX <- sum(likelihood)
  DoubleX <- DoubleX*sigsq[1]^2 + XGrad(sigsq)
  return(DoubleX)
}
#d2whittle/d(sigx^2)(sige^2)
HessianMixed <- function(sigsq)
{
  likelihood <- vector(length = (len/2)-1)
  p <- 1:((len+1)/2-1)
  likelihood[p] = sinsummands[p]/(sigsq[1] + sigsq[2]*sinsummands[p])^2 - 2*sinsummands[p]*modulus[(p+floor(len/2))]/(sigsq[1] + sigsq[2]*sinsummands[p])^3
  Mixed <- sum(likelihood)
  Mixed <- Mixed * prod( sigsq )
  return(Mixed)
}
#d2whittle/d(sige^2)^2
HessianDoubleY <- function(sigsq)
{
  likelihood <- vector(length = (len/2)-1)
  p <- 1:((len+1)/2-1)
  likelihood[p] = (sinsummands[p])^2/(sigsq[1] + sigsq[2]*sinsummands[p])^2 - 2*modulus[(p+floor(len/2))]*(sinsummands[p])^2/(sigsq[1] + sigsq[2]*sinsummands[p])^3
  DoubleY <- sum(likelihood)
  DoubleY <- DoubleY*sigsq[2]^2 + EGrad(sigsq)
  return(DoubleY)
}
#This function forms the Hessian matrix and actually returns the inverse, what is used in the Newton Method.
Hessian <- function(sigsq)
{
  DoubleX <- HessianDoubleX(sigsq)
  MixedXY <- HessianMixed(sigsq)
  DoubleY <- HessianDoubleY(sigsq)
  HMatrix <- matrix(c(DoubleX, MixedXY, MixedXY, DoubleY), nrow = 2, ncol = 2)
  HMatrix <- solve(HMatrix)
  return(HMatrix)
}

#The main MLE, running a Newton Method to maximising the whittlelik function.
#The answer returned is the original parameterisation, on the positive real axis.
MLE <- function(sigsq)
{
  count <- 0
  xold <- sigsq
  xold <- log( xold )  
  xnew <- xold - 0.5*Hessian( exp(xold) )%*%grad( exp(xold) )
  print(xnew)
  while(abs( whittlelik(exp(xold)) - whittlelik(exp(xnew))) > .00000001  ){
    #print(whittlelik(exp(xold)))
    #print(exp(xold))
    xold <- xnew
    xnew <- xold - 0.5*Hessian( exp(xold) )%*%grad( exp(xold) )
  }
  print(whittlelik(exp(xnew)),digits=20)
  return( exp(xnew) )
}


newsinsummands <-vector(length = len-1)
for(p in 1:(len-1)){
  newsinsummands[p] = 2*sin(pi*(p-1-len/2)*delt)
  newsinsummands[p] = (Mod(newsinsummands[p]))^2
}

#These are the L_k coefficients described in the project. We multiply the periodogram by these at each individual frequency k
ShrinkageRatio <- function(maxSigSq, sinsummands)
{
  ratios <- vector(length=(len-1))
  for(k in 1:(len-1)){
    ratios[k] <- maxSigSq[1]/(maxSigSq[1] + maxSigSq[2]*newsinsummands[k])
  }
  plot(ratios[1:len])
  return(ratios)
}

#Returns our final new estimate of the volatility
NewEstimator <- function(modulus, Shrinks){
  volatility <- 0
  volatility <- sum(Shrinks*modulus)
  return(volatility)
}

len <- 10800
trials <- 10
NoiseSD <- 0.0005
trialCount <- 0
StockStart <- 1
VolStart <- 0.04
delt <- 1/len
averageB <- numeric(length = (len-1))
averageS <- numeric(length = (len-1))
EstimatorData <- matrix(data = NA, nrow = 5, ncol = trials)
#The matrix above will be used to store: sum of squares estimate, frequency estimate, realised vol and parameter estimates. 

for(i in 1: trials){
  trialCount <- trialCount+1
  loglist <- Heston.Series(len, StockStart, VolStart)
  StockList <- exp(loglist)
  S.Difference <- numeric(len-1)
  W.Difference <- numeric(len-1)
  N.Difference <- numeric(len-1)
  j <- 2:(len)
  S.Difference[j-1] <- (StockList[j] - StockList[j-1])
  Realised.Vol <- sum(S.Difference^2)
  WhiteNoise <- White.Noise(len, NoiseSD)
  N.Difference[j-1] <- WhiteNoise[j] - WhiteNoise[j-1]
  WorkingList <- StockList+WhiteNoise
  W.Difference[j-1] <- WorkingList[j] - WorkingList[j-1]
  W.Vol <- sum(W.Difference^2) # The sum of squares estimate is here
  plot( StockList, xlab="time", ylab="spot priceS", type="l")
  plot( WorkingList, xlab="time", ylab="spot priceW", type="l")


  modulus <- numeric(length = (len-1))
  modulus <- (fft(W.Difference))
  modulus <- fftshift(modulus) #Note that fft moves the zero frequency to the start of the vector. So we need fftshift to move it to the center, as in DFT representation
  modulus <- Mod(modulus)^2/len #fft has no 1/sqrt(len) factor, so we deal with that here
  
  plot(modulus)
  #plot(modulus.E)
  #plot(modulus.S)
  
  j <- 1:(len-1)

  sigsq <- matrix(nrow = 2, ncol = 1)
  sigsq[1] <- 1.5e-08
  sigsq[2] <- 0.00000025
  maxSigSq <- MLE(sigsq)

  Shrinks <- ShrinkageRatio(maxSigSq, newsinsummands)

  New.Vol <- NewEstimator(modulus, Shrinks)
  EstimatorData[1, trialCount] <- W.Vol
  EstimatorData[2, trialCount] <- Realised.Vol
  EstimatorData[3, trialCount] <- New.Vol
  EstimatorData[4, trialCount] <- maxSigSq[1]
  EstimatorData[5, trialCount] <- maxSigSq[2]
  print(W.Vol)
  print(Realised.Vol)
  cat("Adjusted Vol")
  print(New.Vol)
  print(i)
  plot(Shrinks*modulus)
}