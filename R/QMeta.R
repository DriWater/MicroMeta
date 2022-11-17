
.F.test <- function(x){

  x.stat = -2 * sum(log(x))
  return( 1 - pchisq(x.stat, df = 2 * length(x)) )
}

.simes.test <- function(x){

  return( min(length(x) * x/rank(x)) )

}


########################################
#                                      #
#           One Part Model             #
#                                      #
########################################


.Ei.beta <- function(m, p, beta, X.i, Y.i){

  Ei.out = rep(NA,m)

  for(j in 1:(m-1)){

    Ei.out[j] = exp(beta[((j-1)*p+1):(j*p)] %*% X.i)

  }


  Ei.out[m] = 1




  return (Ei.out)
}


## change on 03/28/2016
.fun.neg.loglik.beta <- function(beta, data){

  Y = data$Y; X = data$X;

  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)

  n.beta = (m - 1)*p
  loglik = 0

  if(length(beta)!=n.beta){

    warning("Dim of initial beta does not match the dim of covariates")

  }else{

    for(i in 1:n){

      E.i = .Ei.beta(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      #Y.pos.index = which(Y[i,]>0)
      #loglik = loglik + Y[i,Y.pos.index] %*% log(P.i[Y.pos.index])
      loglik = loglik + Y[i,] %*% log(P.i)
    }

  }

  return (-loglik)

}

.fun.neg.score.beta <- function(beta, data){

  Y = data$Y; X = data$X;

  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)

  n.beta = (m - 1)*p

  if(length(beta)!=n.beta){

    warning("Dim of initial beta does not match the dim of covariates")

  }else{

    Score.beta = rep(0, n.beta)
    nY = rowSums(Y)

    for(i in 1:n){

      E.i = .Ei.beta(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      Score.beta = Score.beta + kronecker( matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1))

    }

    return (-Score.beta)
  }

}
