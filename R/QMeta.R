
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

.fun.score.i.beta <- function(beta, data){

  Y = data$Y; X = data$X;

  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)

  n.beta = (m - 1)*p

  if(length(beta)!=n.beta){

    warning("Dim of initial beta does not match the dim of covariates")

  }else{

    Score.beta.i = matrix(0, n, n.beta)
    nY = rowSums(Y)

    for(i in 1:n){

      E.i = .Ei.beta(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i

      # add 03/28/2016
      #       if(sum.E.i==0){
      #         P.i = rep(0,m)
      #       }

      Score.beta.i[i,] =  kronecker( matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1) )

    }

    return (Score.beta.i)
  }



}

.fun.hessian.beta <- function(beta, data, save.list=FALSE){

  Y = data$Y; X = data$X

  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  n.beta = (m-1)*p

  if(length(beta)!=n.beta){
    print("Waring: dim of beta is not the same as beta\n")

  }else{

    Hessian.beta = matrix(0, nrow=n.beta, ncol=n.beta)
    nY = rowSums(Y)
    I.beta.list = list()

    for(i in 1:n){

      E.i = .Ei.beta(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i

      ## tmp.beta
      #tmp.beta =  (E.i[-m] %o% E.i[-m])*nY[i]/sum.E.i^2
      tmp.beta =  as.matrix(P.i[-m] %o% P.i[-m])
      diag(tmp.beta) = diag(tmp.beta) - P.i[-m]
      tmp.beta = nY[i] * tmp.beta
      #tmp.beta[is.na(tmp.beta)] = 0  ## add 03/28/2016

      Hessian.beta = Hessian.beta + kronecker( tmp.beta, ( X[i,] %o% X[i,] ) )

      if(save.list){
        I.beta.list[[i]] = tmp.beta

      }

    }


    if(save.list){

      return ( list(Hessian.beta=Hessian.beta, I.beta.list = I.beta.list) )

    }else{

      return (Hessian.beta)
    }

  }


}

.Score.test.stat <- function(Y, X, X.par.index){

  p = ncol(X)

  nY = rowSums(Y)

  n = nrow(Y)
  m = ncol(Y)
  n.beta = (m - 1)*p

  if(sum(X.par.index == 1)){

    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")

  }

  if(is.null(X.par.index) || n==0){
    score.stat.beta = NA
  }else{

    X.reduce = X[,-X.par.index, drop=FALSE]
    p.reduce = p - length(X.par.index)
    par.interest.index.beta =  kronecker( ((0:(m-2))*p), rep(1,length(X.par.index))) + X.par.index

    n.par.interest.beta = length(par.interest.index.beta)

    beta.ini.reduce = rep(0, (p.reduce*(m-1)))

    data.reduce.beta = list(Y=Y, X=X.reduce)

    est.reduce.beta = rep(NA, n.beta)
    est.reduce.beta[par.interest.index.beta] = 0
    # change 04/08/2016
    est.reduce.beta[-par.interest.index.beta] = optim(par=beta.ini.reduce, fn=.fun.neg.loglik.beta, gr=.fun.neg.score.beta, data = data.reduce.beta, method="BFGS")$par


    data.beta = list(Y=Y, X=X)
    Score.reduce.beta = .fun.score.i.beta(est.reduce.beta, data.beta)

    # for resampling: S.beta.list, I.beta.list
    S.beta.list = lapply(1:n, function(j) Score.reduce.beta[j, ((1:(m-1))*p-p+1)])
    tmp = .fun.hessian.beta(est.reduce.beta, data.beta, save.list=TRUE)
    I.beta.list = tmp$I.beta.list

    Hess.reduce.beta = tmp$Hessian.beta
    #Hess.reduce.beta =  .fun.hessian.beta.pos(est.reduce.beta, data.beta)

    # re-organized the score statistics and Hessian matrix
    Score.reduce.reorg = cbind( matrix(Score.reduce.beta[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
    Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.beta[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ),
                              cbind( matrix(Hess.reduce.beta[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))


    A = colSums(Score.reduce.reorg)[1:n.par.interest.beta]

    B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )


    B2 =  matrix(0, n.beta, n.beta)

    # change in 04/08/2016(warning! need to change this step in the resampling function too)
    for(i in 1:n){
      B2 = B2 + Score.reduce.reorg[i,] %o% Score.reduce.reorg[i,]
    }

    B = B1 %*% B2 %*% t(B1)
    score.stat.beta = A %*% ginv(B) %*% A


  }


  return(list(score.stat.beta=score.stat.beta, score.beta = A, est.cov=B, S.beta.list=S.beta.list, I.beta.list=I.beta.list))


}

