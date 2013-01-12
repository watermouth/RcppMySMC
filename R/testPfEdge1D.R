testPfEdge1D <- function(){
  x <- c(rep(0,49),1,rep(0,49),1,rep(0,50))
  y <- 100*cumsum(x) + rpois(n=length(x),lambda=40)
  trueData <- list(x=x, y=y)

  n <- 10000
  output <- pfEdge1D(observations=trueData$y,numOfParticles=n,alpha=0.01,beta=100,gamma=30)
  
  plot(trueData$x)
  points(output$mean, pch=20, col="blue")
}



