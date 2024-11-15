#' Title
#'
#' @param n number of biased sample
#' @param m number of unbiased sample
#' @param N population sizze
#' @param rnull function to generate samples from null distribution, e.g., rexp, rnorm
#' @param bias_type format of bias function
#' @param b a parameter in bias function, usually larger b represents larger bias
#' @param ... additional arguments used in 'rnull', e.g., mean, sd
#'
#' @return a sample of biased (+unbiased) sample
#' @export
#'
#' @examples
#' RnG_UL(n=10,m=2,rnull=rexp,bias_type='971',b=3)

RnG_UL=function(n,m=0,N=-Inf,rnull=rexp,bias_type='971',b=0,...){
  if(N>0) return(NA)

  r=rnull(3*n,...)
  u=runif(3*n)

  y=switch(bias_type,
           '971'=(10*r+1)/(10*r+1+b)*0.8,
           '972'=NA,
           'linear'=NA,
           'mary'=(10*r*(r>0)+1)/(10*r*(r>0)+1+b)*0.8+0.2,
           NA
  )

  r=r[u<y]
  y=rnull(m,...)
  return(list(x=r[1:n],y=y))
}
