#' @title Basic draft of Test1
#' @description Performs the P-test in the paper "(1999) Testing uniformity versus a monotone density".
#'   Calculates a test statistic (Lambda) and a p-value.
#'
#' @param x numeric vector. The input data of observed sample, assumed to be sorted.
#' @param setting character. Specifies the null distribution. Currently only
#'   `'u'` (Uniform) is supported (default 'u').
#' @param c numeric. Tuning parameter for penalty terms (default 0.3).
#' @param B integer. Number of Monte Carlo replicates for p-value calculation
#'   (default 1000). Higher values increase precision but take longer.
#'
#' @return Test results, an object of class `"htest"` containing the following components:
#'   \item{statistic}{The value of the test statistic (Lambda).}
#'   \item{p.value}{The p-value for the test, computed via simulation.}
#'   \item{method}{A character string indicating the test performed.}
#'   \item{data.name}{A character string giving the name of the data.}
#'   \item{parameters}{A named numeric vector containing `c` (tuning parameter).}
#'   \item{null.value}{Description of the null hypothesis distribution.}
#'
#' @export
#' @import pbapply
#' @import stats
#'
#' @examples
#' set.seed(123)
#' # Data close to uniform
#' uniform_data <- sort(runif(50))
#' test_uniform <- test1(uniform_data, c = 0.3, B = 1000)
#' print(test_uniform)
#'
#' # Data with a monotone trend (e.g., Beta(2,1))
#' monotone_data <- sort(rbeta(50, 2, 1))
#' test_monotone <- test1(monotone_data, c = 0.3, B = 1000)
#' print(test_monotone)
#'
#' @references Jiayang Sun. Michael Woodroofe. "Testing uniformity versus a monotone density."
#' Ann. Statist. 27 (1) 338 - 360, February 1999. https://doi.org/10.1214/aos/1018031114
#'
test1=function(x,setting='u',c=0.3,B=1000){
  n=length(x)
  op=pboptions(type="timer")

  if(setting=='u'){ # set null distribution to be U[0,1]
    Fx=function(x){return(x)}
    info='U[0,1]'
    null_sample=function(n){return(runif(n))}
  }else {
    # Placeholder for future settings
    stop("Setting not implemented.")
  }


  gamma3=function(x,err0,an,beta,n){ # use (4) to obtain omega_n_hat, and consequently estimate gamma
    g1=1;g=Inf
    while(abs(g1-g)>err0){
      g=g1
      omega_n_s=rep(NA,n)
      omega_n_s[1]=(n+an)/(g*(1-Fx(x[1]))+beta) # H0: uniform dist F(x[1])
      omega_n_s[2:n]=sapply(2:n,function(i){return((n-i+1)/(g*(1-Fx(x[i]))+beta))}) # i changes
      omega_n=max(omega_n_s)/n
      g1=1+an/n-beta*omega_n
    }
    return(g1)
  }

  est3=function(x,g,an,beta,i,n){ # estimate each omega
    outs=rep(NA,n-i+1)
    if(i>1) an=0 # no need to add an
    if(i==n) return(list(i0=n,out=(n-i+1+an)/(g*(1-Fx(x[i]))+beta)/n))
    outs[1:(n-i)]=sapply(i:(n-1),function(j){return((j-i+1+an)/(g*(Fx(x[j+1])-Fx(x[i]))))}) # j changes
    outs[n-i+1]=(n-i+1+an)/(g*(1-Fx(x[i]))+beta)
    outn=min(outs)
    i0=which(outs==outn)+i-1 # i0 is the index of n1+n2+n3+...+n_i
    ## Note: the reason of returning i0 is that: if i<k<j reaches the max(min(.)) for k,
    ## then this couple of (i,j) also reaches the max(min(.)) for any k<=j. (The property of the greatest convex minorant.)
    ## That's why we use i0 to record the position of the j reaching max(min(.)) for k.
    return(list(i0=i0,out=outn/n)) # out is the estimated omega_k
  }

  pmle3=function(x,n,err0,an,beta){ # estimate each omega_hats
    omegas=index=vector()
    g=gamma3(x,err0,an,beta,n) # estimate gamma
    iout=0;i=1;sum=0
    while(i<=n){
      iout=iout+1
      result=est3(x,g,an,beta,i,n) # estimate omegas[k]=omega_k and index i0
      i0=result[[1]];omegas[iout]=result[[2]]
      sum=sum+omegas[iout]*(Fx(x[i0+1])-Fx(x[i])) # measure the appropriation of our estimation
      index[iout]=i0 # index[i] = n1+n2+n3+...+n_i
      i=i0+1 # continue with the next group
      #cat('omegas:',omegas,'index:',index,'i0:',i0,'i:',i,'iout:',iout,'\n')
    }
    if(abs(1-sum)>1e-4) warning("1 =/= sum =",sum)
    return(list(omegas=omegas,index=index,iout=iout))
  }

  f=function(nouse=NA,n=150,err0=0.00001,c=0.3,dist){ # function to generate test statistic for null/alternative
    cat('null=',info,', n=',n,', c=',c,sep='')
    an=c*sqrt(n) # an=alpha*n, i.e., alpha=c/sqrt(n)=beta
    beta=c/sqrt(n)
    x=c(sort(dist(n)),1) # x[0]=0, x[n+1]=1 ###
    result=pmle3(x,n,err0,an,beta) # estimate index[i] = n1+n2+n3+...+n_i, omegas[k]=omega_hat_k, iout=m
    omegas=result[[1]];index=result[[2]];iout=result[[3]]
    p=an*log(omegas[1])-beta*n*(omegas[iout]-1) # penalty term
    rt=index[1]*log(omegas[1]) # non-penalty term
    if(iout>1) for(i in 2:iout) rt=rt+(index[i]-index[i-1])*log(omegas[i])
    return(p+rt)
  }

  f2=function(x,err0=0.00001,c=0.3){ # function to compute test statistic for input sample
    n=length(x)
    x=c(x,1)
    an=c*sqrt(n) # an=alpha*n, i.e., alpha=c/sqrt(n)=beta
    beta=c/sqrt(n)
    result=pmle3(x,n,err0,an,beta) # estimate index[i] = n1+n2+n3+...+n_i, omegas[k]=omega_hat_k, iout=m
    omegas=result[[1]];index=result[[2]];iout=result[[3]]
    p=an*log(omegas[1])-beta*n*(omegas[iout]-1) # penalty term
    rt=index[1]*log(omegas[1]) # non-penalty term
    if(iout>1) for(i in 2:iout) rt=rt+(index[i]-index[i-1])*log(omegas[i])
    return(p+rt)
  }

  g=function(n,C,B){ # function to compute null values (B simulations for each setting)
    pboptions(nout=B) # number of independent trials
    M=matrix(nrow=B,ncol=length(C))
    #cat('Computing critical value using the null distribution:\n')
    for(j in 1:length(C)){
      M[,j]=pbsapply(1:B,f,n=n,c=C[j],dist=null_sample)
    }
    colnames(M)=C
    #cat('critical values:\n')
    #print(M)
    return(M)
  }

  # critic=g(n=length(x),C=0.3)[1,1]
  null_stats=g(n=length(x),C=c,B=B)[,1]
  lambda=f2(x)
  # if(lambda>critic) print(paste('Since lambda = ',round(lambda,3),' > ',round(critic,3)
  #                               ,' = critical value ,we fail to reject the null (no bias) hypothesis.',sep=''))
  # else print(paste('Since lambda = ',round(lambda,3),' < ',round(critic,3)
  #                               ,' = critical value ,we reject the null (no bias) hypothesis.',sep=''))

  result <- list(
    statistic = structure(lambda, names = "Lambda"),
    parameter = structure(c(c), names = c("c")),
    p.value = sum(null_stats>lambda)/length(null_stats),
    # conf.int = structure(conf.int, conf.level = 0.95),
    # estimate = c('?'),
    null.value = structure(info, names = "distribution"),
    # alternative = 'monotone increasing bias',
    method = 'Penalized Likelihood Ratio Test',
    data.name = deparse(substitute(x))
  )

  class(result) <- "htest"
  message("Simulation complete.")
  return(result)
}
