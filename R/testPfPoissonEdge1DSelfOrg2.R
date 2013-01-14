testPfPoissonEdge1DSelfOrg2 <- function(){
  x <- c(rep(0,100),rep(0,49),1,rep(0,49),1,rep(0,50),1,0,0,0,1,rep(0,100))
#   y <- 100*cumsum(x) + rpois(n=length(x),lambda=40) # + rnorm(n=length(x),mean=seq(from=-10,to=10,along.with=x),sd=5)
  y <- rpois(n=length(x), lambda=100 + 100*cumsum(x))
#   y <- rpois(n=length(x),lambda=seq(from=10,to=100,along.with=x))
  trueData <- list(x=x, y=y)
  
  n <- 1000
  output <- pfPoissonEdge1DSelfOrg2(observations=trueData$y,numOfParticles=n,lambdaMax=10000,sigma_gamma=0.01,gamma0=1,
                                   sigma_alpha=0.001, alpha0=0.1,beta0=1,lambda0=1)
  
  plot(trueData$x)
  points(output$mean, pch=20, col="blue")
  plot(y)
  meanLambda <- output$meanLambda
  points(meanLambda, col="blue", pch=20)
  points(meanLambda + sqrt(output$varLambda), col="green", pch=20)
  points(meanLambda - sqrt(output$varLambda), col="green", pch=20)
  
}

