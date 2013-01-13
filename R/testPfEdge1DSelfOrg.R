testPfEdge1DSelfOrg <- function(){
  x <- c(rep(0,100),rep(0,49),1,rep(0,49),1,rep(0,50),1,0,0,0,1,rep(0,100))
  y <- -100*cumsum(x) + rpois(n=length(x),lambda=40)
  trueData <- list(x=x, y=y)
  
  n <- 10000
  output <- pfEdge1DSelfOrg(observations=trueData$y,numOfParticles=n,alpha0=0.01,beta0=100,gamma0=20)
  
  plot(trueData$x)
  points(output$mean, pch=20, col="blue")
}


