#' @title Simulation of gbm
#' @description Simulation of gbm with euler method.
#' @param a the paramter of gbm
#' @param b the paramter of gbm
#' @param x0 the start point of gbm
#' @param t the time of gbm
#' @param delta_t the step of simulation
#' @return a path of gbm
#' @importFrom stats rnorm
#' @export
sde_euler<-function(a,b,x0,t,delta_t){
  n<-floor(t/delta_t)
  x<-vector()
  x[1]=x0
  for (i in 2:n) {
    r<-rnorm(1,0,sqrt(delta_t))
    x[i]<-x[i-1]+a*x[i-1]*delta_t+b*x[i-1]*r
  }
  return(x)
}


#' @title Simulation of gbm
#' @description Simulation of gbm with milstein method.
#' @param a the paramter of gbm
#' @param b the paramter of gbm
#' @param x0 the start point of gbm
#' @param t the time of gbm
#' @param delta_t the step of simulation
#' @return a path of gbm
#' @importFrom stats rnorm
#' @export
sde_milstein<-function(a,b,x0,t,delta_t){
  n<-floor(t/delta_t)
  x<-vector()
  x[1]=x0
  for (i in 2:n) {
    r<-rnorm(1,0,sqrt(delta_t))
    R<-0.5*b^2*x[i-1]*(r^2-delta_t)
    x[i]<-x[i-1]+a*x[i-1]*delta_t+b*x[i-1]*r+R
  }
  return(x)
}

#' @title Simulation of gbm
#' @description Simulation of gbm with more R in taylor series.
#' @param a the paramter of gbm
#' @param b the paramter of gbm
#' @param x0 the start point of gbm
#' @param t the time of gbm
#' @param delta_t the step of simulation
#' @return a path of gbm
#' @importFrom stats rnorm
#' @export
sde_taylor<-function(a,b,x0,t,delta_t){
  n<-floor(t/delta_t)
  x<-vector()
  x[1]=x0
  for (i in 2:n) {
    r<-rnorm(1,0,sqrt(delta_t))
    R1<-0.5*b^2*x[i-1]*(r^2-delta_t)
    R2<-0.5*a^2*x[i-1]*delta_t^2+a*b*x[i-1]*r*delta_t
    x[i]<-x[i-1]+a*x[i-1]*delta_t +b*x[i-1]*r+R1+R2
  }
  return(x)
}

#' @title Convergence order of the simulation
#' @description To get the convergence order of simulation with different method.
#' @param a the paramter of gbm
#' @param b the paramter of gbm
#' @param x0 the start point of gbm
#' @param t the time of gbm
#' @param type the method of simulation
#' @param N the number of paths
#' @importFrom stats rnorm lm
#' @importFrom graphics lines
#' @importFrom utils tail
#' @return cnvergence order
#' @export
converge_order<- function(N,a,b,x0,t,type) {
  EXT<-x0*exp(a*t)
  delta<-2^c(-1,-2,-3,-4,-5,-6)
  xt<-vector()
  
  if(type=="euler"){
    for (i in 1:6) {
      delta_t<-delta[i]
      xt[i]<-mean(replicate(N,tail(sde_euler(a,b,x0,t,delta_t),1)))
    }
  }else if(type=="milstein"){
    for (i in 1:6) {
      delta_t<-delta[i]
      xt[i]<-mean(replicate(N,tail(sde_milstein(a,b,x0,t,delta_t),1)))
    }
  }else if(type=="taylor"){
    for (i in 1:6) {
      delta_t<-delta[i]
      xt[i]<-mean(replicate(N,tail(sde_taylor(a,b,x0,t,delta_t),1)))
  }
  }
  err<-abs(xt-EXT)
  delta<-log(delta)
  err<-log(err)
  plot(delta,err,xlab = "log(delta_t)",ylab = "log(errors)")
  lines(delta,err)
  order<-lm(formula = err~delta)
  order<-as.numeric(order$coefficients[2])
  return(order)
  }
  
  