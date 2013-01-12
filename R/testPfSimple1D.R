testPfSimple1D <- function(){
  stateDiff <- rnorm(n=100,mean=0,sd=1)
  x <- cumsum(stateDiff)
  y <- x + rnorm(n=length(x),mean=0,sd=2)
  trueData <- list(x=x, y=y)
  
  n <- 1000
  output <- pfSimple1D(observations=trueData$y,numOfParticles=n,std_x=2,std_y=2)
  p1sigma <- output$mean + output$sd
  m1sigma <- output$mean - output$sd
  
  plot(trueData$x)
  points(output$mean, pch=20)
  points(p1sigma, col="green", pch=20)
  points(m1sigma, col="green", pch=20)
  
  print(length(which(m1sigma < trueData$x & trueData$x < p1sigma)) / length(trueData$x))
}

