library(MASS)
library(data.table)
library(CompQuadForm)
library(geepack)
library(dplyr)
library(Rcpp)
library(RcppArmadillo)
library(brglm2)
library(psych)

sourceCpp("D:/QCAT_Meta_new/resampling.cpp")

.F.test <- function(x) {
  # Fisher's p-value combination
  x.stat <- -2 * sum(log(x))
  return(1 - pchisq(x.stat, df = 2 * length(x)))
}

.diag2 <- function(x){
  # transform the numeric into diag matrix
  if(length(x)>1){
    return(diag(x))
  }else{
    return(as.matrix(x))
    
  }
  
}

########################################
#                                      #
#           Positive Part Model             #
#                                      #
########################################


.Ei.beta <- function(m, p, beta, X.i, Y.i) {
  # calculate the exponential of beta times X
  Ei.out <- rep(NA, m)
  for (j in 1:(m - 1)) {
    Ei.out[j] = exp(crossprod(beta[((j-1)*p+1):(j*p)], X.i))
  }
  
  
  Ei.out[m] <- 1 # set the m-th taxa as reference
  
  
  
  
  return(Ei.out)
}

.fun.neg.loglik.beta <- function(beta, data) {
  # return  the negative log likelihood of estimated pi for each taxa
  Y <- data$Y
  X <- data$X
  
  n <- nrow(Y)
  m <- ncol(Y)
  p <- ncol(X)
  
  n.beta <- (m - 1) * p
  loglik <- 0
  
  # check the dimension of beta
  if (length(beta) != n.beta) {
    warning("Dim of initial beta does not match the dim of covariates")
  } else {
    for (i in 1:n) {
      E.i <- .Ei.beta(m, p, beta, X[i, ], Y[i, ])
      sum.E.i <- sum(E.i)
      P.i <- E.i / sum.E.i
      loglik <- loglik + crossprod(Y[i,], log(P.i))
    }
  }
  
  return(-loglik) # return the log likelihood of beta
}

.fun.neg.score.beta <- function(beta, data) {
  
  # Calculate the summation of score statistics based on the estimated the beta value
  
  Y <- data$Y
  X <- data$X
  
  n <- nrow(Y)
  m <- ncol(Y)
  p <- ncol(X)
  
  n.beta <- (m - 1) * p
  
  #check the dimension of beta
  if (length(beta) != n.beta) {
    warning("Dim of initial beta does not match the dim of covariates")
  } else {
    Score.beta <- rep(0, n.beta)
    nY <- rowSums(Y)
    
    #
    for (i in 1:n) {
      E.i <- .Ei.beta(m, p, beta, X[i, ], Y[i, ])
      sum.E.i <- sum(E.i)
      P.i <- E.i / sum.E.i
      Score.beta <- Score.beta + kronecker(matrix(Y[i, -m] - nY[i] * P.i[-m], ncol = 1), matrix(X[i, ], ncol = 1))
    }
    #  the negative score value for beta function (vector form) (for later optimize function use)
    return(-Score.beta)
  }
}

.fun.score.i.beta <- function(beta, data) {
  # return a vector of score statistics with each element corresponding to one subject.
  Y <- data$Y
  X <- data$X
  
  n <- nrow(Y)
  m <- ncol(Y)
  p <- ncol(X)
  
  n.beta <- (m - 1) * p
  # check the dimension of beta
  if (length(beta) != n.beta) {
    warning("Dim of initial beta does not match the dim of covariates")
  } else {
    Score.beta.i <- matrix(0, n, n.beta)
    nY <- rowSums(Y)
    
    for (i in 1:n) {
      E.i <- .Ei.beta(m, p, beta, X[i, ], Y[i, ])
      sum.E.i <- sum(E.i)
      P.i <- E.i / sum.E.i
      #  the score value for beta function (matrix form) (only for parameter of interest)
      Score.beta.i[i, ] <- kronecker(matrix(Y[i, -m] - nY[i] * P.i[-m], ncol = 1), matrix(X[i, ], ncol = 1))
    }
    
    return(Score.beta.i)
  }
}

.fun.hessian.beta <- function(beta, data, save.list = FALSE) {
  # calculate the hessian matrix for score statistics
  Y <- data$Y
  X <- data$X
  
  n <- nrow(Y)
  m <- ncol(Y)
  p <- ncol(X)
  n.beta <- (m - 1) * p
  
  # check the dimension of beta
  if (length(beta) != n.beta) {
    print("Waring: dim of beta is not the same as beta\n")
  } else {
    Hessian.beta <- matrix(0, nrow = n.beta, ncol = n.beta) # Initialize the Hessian matrix for beta
    nY <- rowSums(Y)
    I.beta.list <- list()
    
    for (i in 1:n) {
      E.i <- .Ei.beta(m, p, beta, X[i, ], Y[i, ])
      sum.E.i <- sum(E.i)
      P.i <- E.i / sum.E.i
      
      
      tmp.beta <- as.matrix(P.i[-m] %o% P.i[-m])
      diag(tmp.beta) <- diag(tmp.beta) - P.i[-m]
      tmp.beta <- nY[i] * tmp.beta
      
      Hessian.beta <- Hessian.beta + kronecker(tmp.beta, (X[i, ] %o% X[i, ]))
      
      if (save.list) {
        I.beta.list[[i]] <- tmp.beta
      }
    }
    
    
    if (save.list) {
      # if save.list = TRUE, this will be used for resampling test
      return(list(Hessian.beta = Hessian.beta, I.beta.list = I.beta.list))
    } else {
      return(Hessian.beta)
    }
  }
}

.Score.test.stat <- function(Y, X, X.par.index) {
  # generate the summary statistics for later meta analysis
  p <- ncol(X)
  nY <- rowSums(Y)
  n <- nrow(Y)
  m <- ncol(Y)
  n.beta <- (m - 1) * p
  
  if(is.null(colnames(Y))){
    colnames(Y) <- paste0(rep('V',m),c(1:m))
  }
  
  if (sum(X.par.index == 1)) {
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }
  
  if (is.null(X.par.index) || n == 0) {
    score.stat.beta <- NA
  } else {
    X.reduce <- X[, -X.par.index, drop = FALSE]
    p.reduce <- p - length(X.par.index)
    
    # the index of parameters of interest
    par.interest.index.beta <- kronecker(((0:(m - 2)) * p), rep(1, length(X.par.index))) + X.par.index
    
    n.par.interest.beta <- length(par.interest.index.beta)
    
    beta.ini.reduce <- rep(0, (p.reduce * (m - 1))) # initialize the those elements of beta which we are interested in
    
    data.reduce.beta <- list(Y = Y, X = X.reduce)
    
    est.reduce.beta <- rep(NA, n.beta)
    est.reduce.beta[par.interest.index.beta] <- 0
    
    # Use the optimize function to estimate beta value which we are interested in
    # est.reduce.beta[-par.interest.index.beta] <- optim(par = beta.ini.reduce, fn = .fun.neg.loglik.beta, gr = .fun.neg.score.beta, data = data.reduce.beta, method = "BFGS")$par
    suppressWarnings(est.reduce.beta[-par.interest.index.beta] <- c(t(coef(brmultinom(Y ~ X - 1, data = data.reduce.beta, type = "AS_mean", ref = m)))))
    
    data.beta <- list(Y = Y, X = X)
    
    # Calculate the Score value for beta parameter of interest
    Score.reduce.beta <- .fun.score.i.beta(est.reduce.beta, data.beta)
    
    # for resampling: S.beta.list, I.beta.list
    S.beta.list <- lapply(1:n, function(j) Score.reduce.beta[j, ((1:(m - 1)) * p - p + 1)])
    tmp <- .fun.hessian.beta(est.reduce.beta, data.beta, save.list = TRUE)
    I.beta.list <- tmp$I.beta.list
    
    Hess.reduce.beta <- tmp$Hessian.beta
    
    # re-organized the score statistics and Hessian matrix based on parameters of interest
    Score.reduce.reorg <- cbind(matrix(Score.reduce.beta[, par.interest.index.beta], ncol = n.par.interest.beta), matrix(Score.reduce.beta[, -par.interest.index.beta], ncol = n.beta - n.par.interest.beta))
    Hess.reduce.reorg <- rbind(
      cbind(matrix(Hess.reduce.beta[par.interest.index.beta, par.interest.index.beta], nrow = n.par.interest.beta), matrix(Hess.reduce.beta[par.interest.index.beta, -par.interest.index.beta], nrow = n.par.interest.beta)),
      cbind(matrix(Hess.reduce.beta[-par.interest.index.beta, par.interest.index.beta], nrow = n.beta - n.par.interest.beta), matrix(Hess.reduce.beta[-par.interest.index.beta, -par.interest.index.beta], nrow = n.beta - n.par.interest.beta))
    )
    
    # re-organize the test statistics
    A <- colSums(Score.reduce.reorg)[1:n.par.interest.beta]
    
    B1 <- ginv(Hess.reduce.reorg[(1:n.par.interest.beta), (1:n.par.interest.beta)] - crossprod(t(Hess.reduce.reorg[(1:n.par.interest.beta), ((n.par.interest.beta + 1):n.beta)]),
                                                                                               crossprod(t(ginv(Hess.reduce.reorg[((n.par.interest.beta + 1):n.beta), ((n.par.interest.beta + 1):n.beta)])), Hess.reduce.reorg[((n.par.interest.beta + 1):n.beta), (1:n.par.interest.beta)])))
    
    beta.hat <- crossprod(t(B1), A)
    
    U <- Score.reduce.reorg[ ,1:n.par.interest.beta] - crossprod(t(Score.reduce.reorg[ ,((n.par.interest.beta + 1):n.beta)]),
                                                                 tcrossprod(ginv(Hess.reduce.reorg[((n.par.interest.beta + 1):n.beta), ((n.par.interest.beta + 1):n.beta)]), Hess.reduce.reorg[(1:n.par.interest.beta), ((n.par.interest.beta + 1):n.beta)] ))
    
    B2 <- matrix(0, n.par.interest.beta, n.par.interest.beta)
    
    for (i in 1:n) {
      B2 <- B2 + U[i, ] %o% U[i, ]
    }
    
    cov.beta = crossprod(t(B1), crossprod(t(B2), B1))
  }
  
  # save those summary statistics for later use
  return(list(score.beta = beta.hat, est.cov = cov.beta, S.beta.list = S.beta.list, I.beta.list = I.beta.list))
}

.Score.test.stat.meta.4Gresampling <- function(X.perm.list, X.par.index, n.par.interest.beta, col.index.list, S.beta.list.meta, I.beta.list.meta, Method = "FE-MV") {
  stu.num <- length(X.perm.list) #  the total number of studies for meta analysis
  
  # initialize those statistics
  score.stat.beta <- NULL
  score.beta <- NULL
  est.cov <- NULL
  # initialize the score statistics and estimate covariance matrix for meta-analysis
  score.beta.meta <- rep(0, n.par.interest.beta)
  est.cov.meta <- matrix(0, nrow = n.par.interest.beta, ncol = n.par.interest.beta)
  
  for (j in c(1:stu.num)) {
    # resampling in each group
    S.beta.list <- S.beta.list.meta[[j]]
    I.beta.list <- I.beta.list.meta[[j]]
    X.perm <- X.perm.list[[j]]
    p <- ncol(X.perm)
    n <- nrow(X.perm)
    m.beta <- length(S.beta.list[[1]])
    n.beta <- m.beta * p
    par.interest.index.beta <- kronecker(((0:(m.beta - 1)) * p), rep(1, length(X.par.index))) + X.par.index
    n.par.beta.interest <- length(par.interest.index.beta)
    Score.reduce.beta.perm <- matrix(0, n, n.beta) ## initialize permuted Score statistics and Hessian matrix
    Hess.reduce.beta.perm <- matrix(0, n.beta, n.beta)
    for (i in 1:n) {
      ###################################################
      #                                                 #
      #         Beta part: resampling Score test        #
      #                                                 #
      ###################################################
      Score.reduce.beta.perm[i, ] <- Score.reduce.beta.perm[i, ] + kronecker(matrix(S.beta.list[[i]], ncol = 1), matrix(X.perm[i, ], ncol = 1))
      
      Hess.reduce.beta.perm <- Hess.reduce.beta.perm + kronecker(I.beta.list[[i]], (X.perm[i, ] %o% X.perm[i, ]))
      #     if(sum(is.na(Hess.reduce.beta.perm))>0){
      #       print(i); break;
      #
      #     }
    }
    ###################################################
    #                                                 #
    #         Beta part: resampling Score test        #
    #                                                 #
    ###################################################
    Score.reduce.beta.perm.reorg <- cbind(matrix(Score.reduce.beta.perm[, par.interest.index.beta], ncol = n.par.beta.interest), matrix(Score.reduce.beta.perm[, -par.interest.index.beta], ncol = n.beta - n.par.beta.interest))
    Hess.reduce.beta.perm.reorg <- rbind(
      cbind(matrix(Hess.reduce.beta.perm[par.interest.index.beta, par.interest.index.beta], nrow = n.par.beta.interest), matrix(Hess.reduce.beta.perm[par.interest.index.beta, -par.interest.index.beta], nrow = n.par.beta.interest)),
      cbind(matrix(Hess.reduce.beta.perm[-par.interest.index.beta, par.interest.index.beta], nrow = n.beta - n.par.beta.interest), matrix(Hess.reduce.beta.perm[-par.interest.index.beta, -par.interest.index.beta], nrow = n.beta - n.par.beta.interest))
    )
    
    
    # re-organized the score statistics and estimate covariance matrix based on the parameter of interest
    A <- colSums(Score.reduce.beta.perm.reorg)[1:n.par.beta.interest]
    
    B1 <- ginv(Hess.reduce.beta.perm.reorg[(1:n.par.beta.interest), (1:n.par.beta.interest)] - crossprod(t(Hess.reduce.beta.perm.reorg[(1:n.par.beta.interest), ((n.par.beta.interest + 1):n.beta)]),
                                                                                                         crossprod(t(ginv(Hess.reduce.beta.perm.reorg[((n.par.beta.interest + 1):n.beta), ((n.par.beta.interest + 1):n.beta)])), Hess.reduce.beta.perm.reorg[((n.par.beta.interest + 1):n.beta), (1:n.par.beta.interest)])))
    
    beta.hat <- crossprod(t(B1), A)
    
    U <- Score.reduce.beta.perm.reorg[ ,1:n.par.beta.interest] - crossprod(t(Score.reduce.beta.perm.reorg[ ,((n.par.beta.interest + 1):n.beta)]),
                                                                           tcrossprod(ginv(Hess.reduce.beta.perm.reorg[((n.par.beta.interest + 1):n.beta), ((n.par.beta.interest + 1):n.beta)]), Hess.reduce.beta.perm.reorg[(1:n.par.beta.interest), ((n.par.beta.interest + 1):n.beta)] ))
    
    B2 <- matrix(0, n.par.beta.interest, n.par.beta.interest)
    
    for (i in 1:n) {
      B2 <- B2 + U[i, ] %o% U[i, ]
    }
    cov.beta = crossprod(t(B1), crossprod(t(B2), B1))
    
    # Generate the test statistics for each study
    score.beta[[j]] <- beta.hat # restore the score statistics and estimated covariance matrix in lists which are of necessity for some meta-analysis methods
    est.cov[[j]] <- cov.beta
    idx <- col.index.list[[j]] # specify the interested parameters' index for each study
    score.beta.meta[idx] =  score.beta.meta[idx] +  crossprod(ginv(cov.beta), beta.hat)
    est.cov.meta[idx,idx] =  est.cov.meta[idx,idx] + ginv(cov.beta)
  }
  # save the index of those elements that have values greater than zero in score.beta.meta vector
  save.index.pos = which(abs(score.beta.meta) >= 1e-7) # the length of index is the total number of beta parameters which we are used in our meta-analysis
  n.par.save.beta = length(save.index.pos)
  score.beta.meta = score.beta.meta[save.index.pos]
  est.cov.meta = est.cov.meta[save.index.pos,save.index.pos]
  est.cov.inv <- ginv(est.cov.meta)
  if(Method == "FE-MV"){
    score.stat.meta.perm = crossprod( score.beta.meta,crossprod(t(est.cov.inv), score.beta.meta))
  }
  if(Method == "FE-VC"){
    score.stat.meta.perm = crossprod(score.beta.meta) #SKAT-VC
  }
  if(Method == "RE-MV"){
    U.theta = 0
    V.theta = 0
    for( i in 1:stu.num ){
      est.inv = ginv(est.cov[[i]])
      U.theta = U.theta + 1/2 * crossprod(crossprod(t(est.inv), score.beta[[i]])) - 1/2 * tr(est.inv)
      V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
    }
    score.stat.meta.perm = crossprod( score.beta.meta,crossprod(t(est.cov.inv), score.beta.meta)) + U.theta^2/V.theta
  }
  if(Method == "RE-VC"){
    U.theta = 0
    V.theta = 0
    U.tau = 0
    est.inv.sum = matrix(0, nrow = n.par.interest.beta, ncol = n.par.interest.beta)
    for( i in 1:stu.num ){
      idx = col.index.list[[i]]
      est.inv = ginv(est.cov[[i]])# calculate the generalized inverse for each estimate covariance for each study
      U.theta = U.theta +  1/2 * crossprod(crossprod(t(est.inv), score.beta[[i]])) - 1/2 * tr(est.inv)
      # when \tau matrix and W matrix is identity the elements in upper left, bottom right as well as bottom left are the same in RE-VC test
      V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
      est.inv.sum[idx,idx] = est.inv.sum[idx,idx] + est.inv
    }
    V.tau = 1/2 * tr(crossprod(est.inv.sum))
    U.tau = 1/2 * crossprod(score.beta.meta) - 1/2 * tr(est.inv.sum)
    score.stat.meta.perm = crossprod(c(U.tau, U.theta), crossprod(ginv(matrix(c(V.tau,rep(V.theta,3)),ncol = 2)), c(U.tau, U.theta)))
  }
  return(as.numeric(score.stat.meta.perm))  # return the test statistics
}

.Score.test.stat.meta.4Gresampling.c <- function(X.perm.list, X.par.index, n.par.interest.beta, col.index.list, S.beta.list.meta, I.beta.list.meta, Method = "FE-MV"){
  tmp = score_test_stat_meta_resampling_c(X.perm.list, col.index.list, S.beta.list.meta, I.beta.list.meta, X.par.index, n.par.interest.beta)
  est.cov.meta = tmp$est_cov_meta
  score.beta.meta = tmp$score_beta_meta
  est.cov = tmp$est_cov
  score.beta = tmp$score_beta
  # save the index of those elements that have values greater than zero in score.beta.meta vector
  save.index.pos = which(abs(score.beta.meta) >= 1e-7)
  score.beta.meta = score.beta.meta[save.index.pos]
  est.cov.meta = est.cov.meta[save.index.pos,save.index.pos]
  est.cov.inv = ginv(est.cov.meta)
  stu.num = length(X.perm.list)
  if(Method == "FE-MV"){
    score.stat.meta.perm = crossprod( score.beta.meta,crossprod(t(est.cov.inv), score.beta.meta))
  }
  if(Method == "FE-VC"){
    score.stat.meta.perm = crossprod(score.beta.meta) #SKAT-VC
  }
  if(Method == "RE-MV"){
    U.theta = 0
    V.theta = 0
    for( i in 1:stu.num ){
      est.inv = ginv(est.cov[[i]])
      U.theta = U.theta + 1/2 * crossprod(crossprod(t(est.inv), score.beta[[i]])) - 1/2 * tr(est.inv)
      V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
    }
    score.stat.meta.perm = crossprod( score.beta.meta,crossprod(t(est.cov.inv), score.beta.meta)) + U.theta^2/V.theta
  }
  if(Method == "RE-VC"){
    U.theta = 0
    V.theta = 0
    U.tau = 0
    est.inv.sum = matrix(0, nrow = n.par.interest.beta, ncol = n.par.interest.beta)
    for( i in 1:stu.num ){
      idx = col.index.list[[i]]
      est.inv = ginv(est.cov[[i]])# calculate the generalized inverse for each estimate covariance for each study
      U.theta = U.theta +  1/2 * crossprod(crossprod(t(est.inv), score.beta[[i]])) - 1/2 * tr(est.inv)
      # when \tau matrix and W matrix is identity the elements in upper left, bottom right as well as bottom left are the same in RE-VC test
      V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
      # U.tau = U.tau - 1/2 * tr(est.inv)
      est.inv.sum[idx,idx] = est.inv.sum[idx,idx] + est.inv
    }
    V.tau = 1/2 * tr(crossprod(est.inv.sum))
    U.tau = 1/2 * crossprod(score.beta.meta) - 1/2 * tr(est.inv.sum)
    score.stat.meta.perm = crossprod(c(U.tau, U.theta), crossprod(ginv(matrix(c(V.tau,rep(V.theta,3)),ncol = 2)), c(U.tau, U.theta)))
  }
  return(as.numeric(score.stat.meta.perm))
}

.resample.work.one.meta <- function(X.list, X.par.index, n.par.interest.beta, col.index.list, score.stat.meta, S.beta.list.meta, I.beta.list.meta, start.nperm, end.nperm, n.one, one.acc, use.cpp, Method = "FE-MV") {
  n.one.new <- n.one
  one.acc.new <- one.acc
  # adaptive permutation test
  for (k in start.nperm:end.nperm) {
    X.perm.list = X.list
    for(p in 1:length(X.list)){
      idx = sample(1:nrow(X.list[[p]])) # sampling the index and reset the design matrix based on the new index
      X.perm.list[[p]][,X.par.index] = X.list[[p]][idx,X.par.index]
    }
    # get the permutation test statistics
    if(use.cpp){
      score.stat.meta.perm <- try(.Score.test.stat.meta.4Gresampling.c(X.perm.list, X.par.index, n.par.interest.beta, col.index.list, S.beta.list.meta, I.beta.list.meta, Method = Method))
    }
    else{
      score.stat.meta.perm <- try(.Score.test.stat.meta.4Gresampling(X.perm.list, X.par.index, n.par.interest.beta, col.index.list, S.beta.list.meta, I.beta.list.meta, Method = Method))
    }
    if (!("try-error" %in% class(score.stat.meta.perm))) {
      n.one.new <- n.one.new + 1 # if score.stat.meta.perm exists, then n.one.new + 1
      # if the permutation test statistics greater than original test statistics, cnt + 1
      if (score.stat.meta.perm >= score.stat.meta) {
        one.acc.new <- one.acc.new + 1 # if score.stat.meta.perm >= score.stat.metam one.acc.new + 1
      }
    }
  }
  # adaptive adjusting based on the permuation results
  # if the total number of permutation results which are greater than original results too small, enlarge the
  # number of total iterations
  if (one.acc.new < 1) {
    next.end.nperm <- (end.nperm + 1) * 100 - 1
    flag <- 1
  } else if (one.acc.new < 10) {
    next.end.nperm <- (end.nperm + 1) * 10 - 1
    flag <- 1
  }
  
  else {
    next.end.nperm <- (end.nperm + 1) - 1
    flag <- 0
  }
  
  return(list(n.one.new = n.one.new, one.acc.new = one.acc.new, flag = flag, next.end.nperm = next.end.nperm))
}

.Score.test.meta <- function(Y.list, X.list, X.par.index, seed=11, resample=FALSE, n.replicates=NULL, use.cpp = F, Method = "FE-MV"){
  stu.num = length(X.list)
  p.par = length(X.par.index)
  m = ncol(Y.list[[1]])
  p = ncol(X.list[[1]])
  n.par.interest.beta = (m-1)*length(X.par.index)
  
  remove.study = NULL
  for(i in 1:stu.num){
    remove.rows.idx = which(rowSums(Y.list[[i]]) == 0)
    if(length(remove.rows.idx) == nrow(Y.list[[i]])){
      remove.study = append(remove.study, i)
      next
    }
    if(length(remove.rows.idx)>0){
      Y.list[[i]] <- Y.list[[i]][-remove.rows.idx, ,drop=FALSE]
      X.list[[i]] <- X.list[[i]][-remove.rows.idx, , drop=FALSE]
    }
  }
  if(!is.null(remove.study)){
    if(length(remove.study) == stu.num){
      stop("Error: Not proper data")
    }else{
      Y.list <- Y.list[-remove.study]
      X.list <- X.list[-remove.study]
      stu.num = length(X.list)
    }
  }
  
  if(sum(X.par.index == 1)){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }
  # check the parameters of interest to ensure the intercept term is not in it
  if(! Method %in% c("FE-MV","FE-VC","RE-MV","RE-VC")){
    stop("Error: Please Choose a Proper Meta-analysis Method")
  }
  # the taxa of each study should be the same
  if (length(unique(sapply(1:stu.num,function(j) ncol(Y.list[[j]]))))!=1){
    stop("Error: The taxon in each study should be the same")
  }
  # initialize the test statistics
  if(is.null(X.par.index)){
    stop("Error: Please provide the index(es) of covariate(s) of interest")
  }
  
  # initialize for later use
  score.stat.beta = NULL
  score.beta = NULL
  score.pvalue = NA
  score.stat.meta = NA
  df = NA
  est.cov = list()
  S.beta.list.meta = list()
  I.beta.list.meta = list()
  col.index.list = list()
  remove.index = NULL
  # initialize the score statistics and estimate covariance matrix for meta analysis
  score.beta.meta  = rep(0,n.par.interest.beta) ## A
  est.cov.meta = matrix(0, nrow = n.par.interest.beta, ncol = n.par.interest.beta) ## B
  ava.cnt = 0
  stu.num = length(Y.list)
  for(i in 1:stu.num){
    Y = Y.list[[i]]
    col.index = which(colSums(Y)>0) # keep those taxa which have more than one observation in each study
    Y = Y[,col.index, drop = FALSE]
    X = X.list[[i]]
    # nY = rowSums(Y)
    # # remove 03/28/2016
    # nY.index = which(nY==0)
    # if(length(nY.index)>0){
    #   Y = Y[-nY.index, , drop=FALSE]
    #   X = X[-nY.index, , drop=FALSE]
    # }
    if(length(col.index) <= 1){
      remove.index = append(remove.index,i)
      next
    }
    col.index = col.index[-length(col.index)]
    n.beta = (m - 1)*p
    tmp.one = try(.Score.test.stat(Y, X, X.par.index))
    if( "try-error" %in% class(tmp.one) ){
      remove.index = append(remove.index,i)
      next
    }else{
      # number of study which can get score statistics
      ava.cnt = ava.cnt + 1
      # get the index for parameter of interest after score statistics and estimate covariance matrix are reorginzed
      # different across studies because of different column numbers
      idx = kronecker((col.index-1)*p.par, rep(1,p.par)) + c(1:p.par)
      col.index.list[[ava.cnt]] = idx
      score.stat.beta = append(score.stat.beta, tmp.one$score.stat.beta)
      score.beta[[ava.cnt]] = tmp.one$score.beta
      est.cov[[ava.cnt]] = tmp.one$est.cov
      S.beta.list.meta[[ava.cnt]] = tmp.one$S.beta.list
      I.beta.list.meta[[ava.cnt]] = tmp.one$I.beta.list
      score.beta.meta[idx] =  score.beta.meta[idx] +  crossprod(ginv(tmp.one$est.cov), tmp.one$score.beta)
      est.cov.meta[idx,idx] =  est.cov.meta[idx,idx] + ginv(tmp.one$est.cov)
    }
  }
  
  
  if(ava.cnt>0){
    if(length(remove.index) != 0){
      X.list = X.list[-remove.index]
      Y.list = Y.list[-remove.index]
    }
    # save the index of those elements that have values greater than zero in score.beta.meta vector
    save.index.pos = which(abs(score.beta.meta) >= 1e-7) # the length of index is the total number of beta parameters which we are used in our meta-analysis
    n.par.save.beta = length(save.index.pos)
    score.beta.meta = score.beta.meta[save.index.pos]
    est.cov.meta = est.cov.meta[save.index.pos,save.index.pos]
    est.cov.inv = ginv(est.cov.meta)
    if(Method == "FE-MV"){
      score.stat.meta = crossprod( score.beta.meta,crossprod(t(est.cov.inv), score.beta.meta))
      score.pvalue = 1- pchisq(score.stat.meta,df = n.par.save.beta )
      df = n.par.save.beta
    }
    if(Method == "FE-VC"){
      weight.cov.meta = eigen(est.cov.meta)$values #eign.val/sum(eign.val)
      score.stat.meta = crossprod(score.beta.meta)
      score.pvalue = davies(score.stat.meta,weight.cov.meta, h = rep(1,n.par.save.beta), delta = rep(0,n.par.save.beta), sigma = 0, lim = 10000, acc = 0.0001)$Qq
      score.pvalue = ifelse(score.pvalue>0,score.pvalue,0)
      df = n.par.save.beta
    }
    if(Method == "RE-MV"){
      U.theta = 0
      V.theta = 0
      for( i in 1:ava.cnt){
        est.inv = ginv(est.cov[[i]])
        U.theta = U.theta + 1/2 * crossprod(crossprod(t(est.inv), score.beta[[i]])) - 1/2 * tr(est.inv)
        V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
      }
      score.stat.meta = crossprod( score.beta.meta,crossprod(t(est.cov.inv), score.beta.meta)) + U.theta^2/V.theta
    }
    if(Method == "RE-VC"){
      U.theta = 0
      V.theta = 0
      U.tau = 0
      est.inv.sum = matrix(0, nrow = n.par.interest.beta, ncol = n.par.interest.beta)
      for( i in 1:ava.cnt){
        idx = col.index.list[[i]]
        est.inv = ginv(est.cov[[i]])# calculate the generalized inverse for each estimate covariance for each study
        U.theta = U.theta +  1/2 * crossprod(crossprod(t(est.inv), score.beta[[i]])) - 1/2 * tr(est.inv)
        # when \tau matrix and W matrix is identity the elements in upper left, bottom right as well as bottom left are the same in RE-VC test
        V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
        est.inv.sum[idx,idx] = est.inv.sum[idx,idx] + est.inv
      }
      V.tau = 1/2 * tr(crossprod(est.inv.sum))
      U.tau = 1/2 * crossprod(score.beta.meta) - 1/2 * tr(est.inv.sum)
      score.stat.meta = crossprod(c(U.tau, U.theta), crossprod(ginv(matrix(c(V.tau,rep(V.theta,3)),ncol = 2)), c(U.tau, U.theta)))
    }
  }
  beta.meta.results = list(score.stat = score.stat.meta, score.pvalue = score.pvalue, df = df)
  
  # if resample = TRUE then will apply permutation method to get permuted p-value
  # adaptive resampling test
  if(resample){
    
    set.seed(seed)
    if(!is.na(score.stat.meta)){
      
      n.one = 0
      one.acc = 0
      
      start.nperm = 1;
      end.nperm = min(100,n.replicates);
      flag = 1
      while(flag & end.nperm <= n.replicates){
        
        results = .resample.work.one.meta(X.list, X.par.index, n.par.interest.beta, col.index.list, score.stat.meta, S.beta.list.meta, I.beta.list.meta, start.nperm, end.nperm, n.one, one.acc, use.cpp = use.cpp, Method = Method)
        n.one = results$n.one.new
        one.acc = results$one.acc.new
        flag = results$flag
        next.end.nperm = results$next.end.nperm
        
        if(flag){
          start.nperm = end.nperm + 1;
          end.nperm = next.end.nperm;
          
        }
        
        if(start.nperm < n.replicates & end.nperm > n.replicates){
          #warning(paste( "Inaccurate pvalue with", n.replicates, "permutations"))
          results = .resample.work.one.meta(X.list, X.par.index, n.par.interest.beta, col.index.list, score.stat.meta, S.beta.list.meta, I.beta.list.meta, start.nperm, end.nperm, n.one, one.acc, use.cpp = use.cpp, Method = Method)
          n.one = results$n.one.new
          one.acc = results$one.acc.new
          
        }
        
      }
      
      
      #      print(paste("Final number of resamplings: ", n.one) )
      #      if(n.one<end.nperm/2){
      #        print("Number of resamplings too small for one-part test")
      #      }
      
      tmp = (one.acc+1)/(n.one+1) # to avoid n.one be zero # resampling p value
      
      #print(n.one)
      #print(one.acc)
      
    }else{
      
      tmp = NA
    }
    
    beta.meta.results = c(beta.meta.results, score.Rpvalue = tmp)
    
    
  }
  
  return(beta.meta.results)
}


########################################
#                                      #
#           Zero Part Model             #
#                                      #
########################################

.Pi.alpha<-function(m, p, alpha, X.i){
  # calculate the exponential of alpha times X
  Pi.out = rep(NA,m)
  
  for(j in 1:m){
    
    tmp = exp(crossprod(alpha[((j-1)*p+1):(j*p)], X.i))
    if(is.infinite(tmp)){
      Pi.out[j] = 1
    }else{
      Pi.out[j] = tmp/(tmp + 1)
    }
    
  }
  
  # no need for base in GEE method
  return (Pi.out)
}

.fun.score.i.alpha <- function(alpha, data, save.list=FALSE){
  # return  the negative log likelihood of estimated pi for each taxa
  Y = data$Y; Z = data$Z;
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(Z)
  
  vA.list = list()
  Vinv.list = list()
  VY.list = list()
  
  n.alpha = m*p
  # check the dimension of alpha
  if(length(alpha)!=n.alpha){
    
    warning("Dim of initial alpha does not match the dim of covariates")
    
  }else{
    
    Score.alpha.i = matrix(0, n, n.alpha)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      Pi.i = .Pi.alpha(m, p, alpha, Z[i,])
      vA.tmp = Pi.i*(1-Pi.i)
      A.i = .diag2(vA.tmp) # transform into diagonal matrices
      t.D.i = kronecker( A.i, as.matrix(Z[i,], ncol=1) )
      V.i = A.i # independent cor structure
      
      tmp.V.i = ginv(V.i)
      tmp.VY = crossprod(t(tmp.V.i), (Y[i,] - Pi.i))
      Score.alpha.i[i,] = crossprod(t(t.D.i), tmp.VY) # score value
      
      if(save.list){
        vA.list[[i]] = vA.tmp
        Vinv.list[[i]] = tmp.V.i
        VY.list[[i]] = tmp.VY
      }
    }
    
    
  }
  # if save.list = TRUE, this will be used for resampling test
  if(save.list){
    
    return ( list(Score.alpha=Score.alpha.i, vA.list = vA.list, Vinv.list = Vinv.list, VY.list=VY.list) )
    
  }else{
    
    return (Score.alpha.i)
  }
  
}

#fun.hessian.alpha(est.reduce.alpha, data.alpha)
.fun.hessian.alpha <- function(alpha, data){
  
  Y = data$Y; Z = data$Z
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(Z)
  n.alpha = m*p
  
  if(length(alpha)!=n.alpha){
    print("Waring: dim of alpha is not the same as alpha\n")
    
  }else{
    
    Hessian.alpha = matrix(0, nrow=n.alpha, ncol=n.alpha)
    nY = rowSums(Y)
    
    
    for(i in 1:n){
      
      Pi.i = .Pi.alpha(m, p, alpha, Z[i,])
      tmp = Pi.i*(1-Pi.i)
      A.i = .diag2(tmp)
      t.D.i = kronecker( A.i, as.matrix(Z[i,], ncol=1) )
      V.i = A.i # independent cor structure
      
      # the Hessian matrix for GEE model
      Hessian.alpha = Hessian.alpha + crossprod(t(t.D.i), tcrossprod(ginv(V.i), t.D.i))
      
      
    }
    
    return (Hessian.alpha)
    
    
  }
  
  
}

.Score.test.stat.zero <- function(Y0, Z, Z.par.index, cor.stru){
  
  # Z.reduce preserve covariates that are not interested in
  Z.reduce = Z[,-Z.par.index,drop=FALSE]
  n = nrow(Y0)
  m = ncol(Y0)
  p = ncol(Z)
  
  if(is.null(colnames(Y0))){
    colnames(Y0) <- paste0(rep('V',m),c(1:m))
  }
  p.reduce = ncol(Z.reduce)
  outcome = NULL
  id = NULL
  cova.reduce = NULL
  for(i in 1:n){
    
    outcome = c(outcome, Y0[i,])
    index.start = 1
    index.end = p
    
    index.start.reduce = 1
    index.end.reduce = p.reduce
    
    for(j in 1:m){
      
      tmp = rep(0, m*p.reduce)
      tmp[index.start.reduce:index.end.reduce] = Z.reduce[i,]
      cova.reduce = rbind(cova.reduce, tmp )
      index.start.reduce = index.start.reduce + p.reduce
      index.end.reduce = index.end.reduce + p.reduce
      
    }
    
    
    id = c(id, rep(i, m))
  }
  
  #data.full = data.frame(outcome=outcome, cova, id = id, row.names=NULL)
  data.reduce = data.frame(outcome=outcome, cova.reduce, id = id, row.names=NULL)
  #gee.full = geeglm(outcome ~ .  - id - 1, data = data.full, id = factor(id), family="binomial", corstr= "independence")
  ## use geeglm to estimate the parameter not interested in
  gee.reduce = geeglm(outcome ~ . - id - 1, data = data.reduce, id = factor(id), family="binomial", corstr= "independence")
  #wald.test = anova(gee.full, gee.reduce)
  
  
  ########### perform score test
  n.alpha = m * p
  par.interest.index.alpha =  kronecker( ((0:(m-1))*p), rep(1,length(Z.par.index))) + Z.par.index
  n.par.interest.alpha = length(par.interest.index.alpha)
  est.reduce.alpha = rep(NA, n.alpha)
  est.reduce.alpha[par.interest.index.alpha] = 0 # set the est.reduce.aplha corresponding to parameters which are interested in to 0
  est.reduce.alpha[-par.interest.index.alpha] = coef(gee.reduce)
  est.reduce.scale = gee.reduce
  
  data.alpha = list(Y=Y0, Z=Z)
  # estimate the Score statistics for parameter of interest
  tmp = .fun.score.i.alpha(est.reduce.alpha, data.alpha, save.list=TRUE)
  Score.reduce.alpha = tmp$Score.alpha
  # for resampling test
  vA.list = tmp$vA.list
  Vinv.list = tmp$Vinv.list
  VY.list = tmp$VY.list
  
  Hess.reduce.alpha =  .fun.hessian.alpha(est.reduce.alpha, data.alpha)
  # re-organized the score statistics and Hessian matrix according to the index of par.interest.index.alpha
  Score.reduce.reorg = cbind( matrix(Score.reduce.alpha[,par.interest.index.alpha], ncol=n.par.interest.alpha), matrix(Score.reduce.alpha[,-par.interest.index.alpha], ncol=n.alpha - n.par.interest.alpha) )
  Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.alpha[par.interest.index.alpha, par.interest.index.alpha], nrow=n.par.interest.alpha), matrix(Hess.reduce.alpha[par.interest.index.alpha, -par.interest.index.alpha], nrow=n.par.interest.alpha) ),
                            cbind( matrix(Hess.reduce.alpha[-par.interest.index.alpha, par.interest.index.alpha], nrow=n.alpha - n.par.interest.alpha), matrix(Hess.reduce.alpha[-par.interest.index.alpha, -par.interest.index.alpha], nrow= n.alpha - n.par.interest.alpha)))
  
  
  A = colSums(Score.reduce.reorg)[1:n.par.interest.alpha]
  
  B1 <- ginv(Hess.reduce.reorg[(1:n.par.interest.alpha), (1:n.par.interest.alpha)] - crossprod(t(Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha + 1):n.alpha)]),
                                                                                               crossprod(t(ginv(Hess.reduce.reorg[((n.par.interest.alpha + 1):n.alpha), ((n.par.interest.alpha + 1):n.alpha)])), Hess.reduce.reorg[((n.par.interest.alpha + 1):n.alpha), (1:n.par.interest.alpha)])))
  
  alpha.hat <- crossprod(t(B1), A)
  
  
  U <- Score.reduce.reorg[ ,1:n.par.interest.alpha] - crossprod(t(Score.reduce.reorg[ ,((n.par.interest.alpha + 1):n.alpha)]),tcrossprod(ginv(Hess.reduce.reorg[((n.par.interest.alpha + 1):n.alpha),
                                                                                                                                                                ((n.par.interest.alpha + 1):n.alpha)]), Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha + 1):n.alpha)] ))
  
  B2 <- matrix(0, n.par.interest.alpha, n.par.interest.alpha)
  
  for (i in 1:n) {
    B2 <- B2 + U[i, ] %o% U[i, ]
  }
  
  cov.alpha = crossprod(t(B1), crossprod(t(B2), B1))
  
  # save these outcomes for later resampling test
  return(list(score.df.alpha = n.par.interest.alpha,  score.alpha = alpha.hat, est.cov.zero=cov.alpha, vA.list=vA.list, Vinv.list=Vinv.list, VY.list=VY.list )   )
  
}

# Get the permutation statistics by R
.Score.test.stat.zero.meta.4Gresampling <- function(Z.perm.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, vA.list.meta, Vinv.list.meta, VY.list.meta, Method = "FE-MV"){
  
  # the total study number
  stu.num = length(Z.perm.list)
  score.stat.alpha = NULL
  score.alpha = NULL
  est.cov.zero = NULL
  score.alpha.meta  = rep(0,n.par.interest.alpha) ## create a vector which length is the number of parameter of interest to restore the score statistics for our meta analysis
  est.cov.meta = matrix(0, nrow = n.par.interest.alpha, ncol = n.par.interest.alpha) ## ## create a matrix which dimension is the number of parameter of interest to restore the score statistics for our estimate of covariance matrix
  for ( j in 1:stu.num){
    vA.list = vA.list.meta[[j]]
    VY.list = VY.list.meta[[j]]
    Vinv.list = Vinv.list.meta[[j]]
    Z.perm = Z.perm.list[[j]]
    p = ncol(Z.perm.list[[j]])
    m.alpha = length(vA.list[[1]])
    n.alpha = m.alpha*p
    # the index of parameter of interest for each study (different across studies becasue of different column number)
    par.index.alpha =  kronecker( ((0:(m.alpha-1))*p), rep(1,length(Z.par.index))) + Z.par.index # the index of parameter of interest for each study
    # the number of parameter of interest for each study
    n.par.alpha.interest = length(par.index.alpha)
    n = nrow(Z.perm)
    
    # initialize for later use
    Score.reduce.alpha.perm = matrix(0, n, n.alpha )
    Hess.reduce.alpha.perm = matrix(0, n.alpha, n.alpha )
    
    for(i in 1:n){
      
      ###################################################
      #                                                 #
      #         alpha part: resampling Score test        #
      #                                                 #
      ###################################################
      tD.tmp = kronecker(.diag2(vA.list[[i]]), as.matrix(Z.perm[i,], ncol=1))
      
      # the permutated score statistics
      Score.reduce.alpha.perm[i,] = Score.reduce.alpha.perm[i,] + crossprod(t(tD.tmp),VY.list[[i]])
      # the permutated Hessian  matrix
      Hess.reduce.alpha.perm = Hess.reduce.alpha.perm + crossprod(t(tD.tmp), tcrossprod(Vinv.list[[i]], tD.tmp))
      
      
    }
    
    # re-organized the score statistics and Hessian matrix according the index of parameter of interest
    Score.reduce.reorg = cbind( matrix(Score.reduce.alpha.perm[,par.index.alpha], ncol=n.par.alpha.interest), matrix(Score.reduce.alpha.perm[,-par.index.alpha], ncol=n.alpha - n.par.alpha.interest) )
    Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.alpha.perm[par.index.alpha, par.index.alpha], nrow=n.par.alpha.interest), matrix(Hess.reduce.alpha.perm[par.index.alpha, -par.index.alpha], nrow=n.par.alpha.interest) ),
                              cbind( matrix(Hess.reduce.alpha.perm[-par.index.alpha, par.index.alpha], nrow=n.alpha - n.par.alpha.interest), matrix(Hess.reduce.alpha.perm[-par.index.alpha, -par.index.alpha], nrow= n.alpha - n.par.alpha.interest)))
    
    
    A = colSums(Score.reduce.reorg)[1:n.par.alpha.interest]
    
    B1 <- ginv(Hess.reduce.reorg[(1:n.par.alpha.interest), (1:n.par.alpha.interest)] - crossprod(t(Hess.reduce.reorg[(1:n.par.alpha.interest), ((n.par.alpha.interest + 1):n.alpha)]),
                                                                                                 crossprod(t(ginv(Hess.reduce.reorg[((n.par.alpha.interest + 1):n.alpha), ((n.par.alpha.interest + 1):n.alpha)])), Hess.reduce.reorg[((n.par.alpha.interest + 1):n.alpha), (1:n.par.alpha.interest)])))
    
    alpha.hat <- crossprod(t(B1), A)
    
    U <- Score.reduce.reorg[ ,1:n.par.alpha.interest] - crossprod(t(Score.reduce.reorg[ ,((n.par.alpha.interest + 1):n.alpha)]),tcrossprod(ginv(Hess.reduce.reorg[((n.par.alpha.interest + 1):n.alpha),
                                                                                                                                                                  ((n.par.alpha.interest + 1):n.alpha)]), Hess.reduce.reorg[(1:n.par.alpha.interest), ((n.par.alpha.interest + 1):n.alpha)] ))
    
    B2 <- matrix(0, n.par.alpha.interest, n.par.alpha.interest)
    
    for (i in 1:n) {
      B2 <- B2 + U[i, ] %o% U[i, ]
    }
    
    cov.alpha = crossprod(t(B1), crossprod(t(B2), B1))
    score.alpha[[j]] = alpha.hat # restore the score statistics and estimated covariance matrix in lists which are of necessity for some meta-analysis methods
    est.cov.zero[[j]] = cov.alpha
    idx = col.zero.index.list[[j]] # the index of alpha parameter of interest in j-th study
    score.alpha.meta[idx] =  score.alpha.meta[idx] + crossprod(ginv(cov.alpha), alpha.hat) # add according to the index of parameter of interest for each study
    est.cov.meta[idx, idx] =  est.cov.meta[idx, idx] + ginv(cov.alpha)
  }
  # save the index of those elements that have values greater than zero in score.aplha.meta vector
  save.index.zero = which(abs(score.alpha.meta) >= 1e-7)
  n.par.save.alpha = length(save.index.zero) # the length of index is the total number of alpha parameter which we are used in our meta-analysis
  score.alpha.meta =  score.alpha.meta[save.index.zero]
  est.cov.meta =  est.cov.meta[save.index.zero,save.index.zero] # the summation of  estimate covariance for each study
  est.cov.inv = ginv(est.cov.meta)
  if(Method == "FE-MV*"){
    score.stat.alpha.perm = crossprod( score.alpha.meta,crossprod(t(est.cov.inv), score.alpha.meta))
  }
  if(Method == "FE-VC*"){
    score.stat.alpha.perm = crossprod(score.alpha.meta) #SKAT-VC
  }
  if(Method == "RE-MV*"){
    U.theta = 0
    V.theta = 0
    for( i in 1:stu.num ){
      est.inv = ginv(est.cov.zero[[i]])
      U.theta = U.theta + 1/2 * crossprod(crossprod(t(est.inv), score.alpha[[i]])) - 1/2 * tr(est.inv)
      V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
    }
    score.stat.alpha.perm = crossprod( score.alpha.meta,crossprod(t(est.cov.inv), score.alpha.meta)) + U.theta^2/V.theta
  }
  if(Method == "RE-VC*"){
    U.theta = 0
    V.theta = 0
    U.tau = 0
    est.inv.sum = matrix(0, nrow = n.par.interest.alpha, ncol = n.par.interest.alpha)
    for( i in 1:stu.num){
      idx = col.zero.index.list[[i]]
      est.inv = ginv(est.cov.zero[[i]])# calculate the generalized inverse for each estimate covariance for each study
      U.theta = U.theta +  1/2 * crossprod(crossprod(t(est.inv), score.alpha[[i]])) - 1/2 * tr(est.inv)
      # when \tau matrix and W matrix is identity the elements in upper left, bottom right as well as bottom left are the same in RE-VC test
      V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
      est.inv.sum[idx,idx] = est.inv.sum[idx,idx] + est.inv
    }
    V.tau = 1/2 * tr(crossprod(est.inv.sum))
    U.tau = 1/2 * crossprod(score.alpha.meta) - 1/2 * tr(est.inv.sum)
    score.stat.alpha.perm = crossprod(c(U.tau, U.theta), crossprod(ginv(matrix(c(V.tau,rep(V.theta,3)),ncol = 2)), c(U.tau, U.theta)))
  }
  return(as.numeric(score.stat.alpha.perm))
  
}

# Get the permutation statistics by cpp
.Score.test.stat.zero.meta.4Gresampling.c <- function(Z.perm.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, vA.list.meta, Vinv.list.meta, VY.list.meta, Method = "FE-MV"){
  tmp = score_test_stat_zero_meta_resampling_c(Z.perm.list, col.zero.index.list, vA.list.meta, Vinv.list.meta, VY.list.meta, Z.par.index, n.par.interest.alpha)
  est.cov.meta = tmp$est_cov_meta
  score.alpha.meta = tmp$score_alpha_meta
  est.cov.zero = tmp$est_cov
  score.alpha = tmp$score_alpha
  # save the index of those elements that have values greater than zero in score.aplha.meta vector
  save.index.zero = which(abs(score.alpha.meta) >= 1e-7)
  score.alpha.meta = score.alpha.meta[save.index.zero]
  est.cov.meta = est.cov.meta[save.index.zero,save.index.zero]
  est.cov.inv = ginv(est.cov.meta)
  stu.num = length(Z.perm.list)
  if(Method == "FE-MV*"){
    score.stat.alpha.perm = crossprod( score.alpha.meta,crossprod(t(est.cov.inv), score.alpha.meta))
  }
  if(Method == "FE-VC*"){
    score.stat.alpha.perm = crossprod(score.alpha.meta)
  }
  if(Method == "RE-MV*"){
    U.theta = 0
    V.theta = 0
    for( i in 1:stu.num ){
      est.inv = ginv(est.cov.zero[[i]])
      U.theta = U.theta + 1/2 * crossprod(crossprod(t(est.inv), score.alpha[[i]])) - 1/2 * tr(est.inv)
      V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
    }
    score.stat.alpha.perm = crossprod( score.alpha.meta,crossprod(t(est.cov.inv), score.alpha.meta)) + U.theta^2/V.theta
  }
  if(Method == "RE-VC*"){
    U.theta = 0
    V.theta = 0
    U.tau = 0
    est.inv.sum = matrix(0, nrow = n.par.interest.alpha, ncol = n.par.interest.alpha)
    for( i in 1:stu.num){
      idx = col.zero.index.list[[i]]
      est.inv = ginv(est.cov.zero[[i]])# calculate the generalized inverse for each estimate covariance for each study
      U.theta = U.theta +  1/2 * crossprod(crossprod(t(est.inv), score.alpha[[i]])) - 1/2 * tr(est.inv)
      # when \tau matrix and W matrix is identity the elements in upper left, bottom right as well as bottom left are the same in RE-VC test
      V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
      est.inv.sum[idx,idx] = est.inv.sum[idx,idx] + est.inv
    }
    V.tau = 1/2 * tr(crossprod(est.inv.sum))
    U.tau = 1/2 * crossprod(score.alpha.meta) - 1/2 * tr(est.inv.sum)
    score.stat.alpha.perm = crossprod(c(U.tau, U.theta), crossprod(ginv(matrix(c(V.tau,rep(V.theta,3)),ncol = 2)), c(U.tau, U.theta)))
  }
  return(as.numeric(score.stat.alpha.perm))
  
}

# add 10/22/2022 for adaptive resampling
.resample.work.zero.meta <- function(Z.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, score.stat.zero.meta, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, start.nperm, end.nperm, n.zero, zero.acc, use.cpp, Method = "FE-MV"){
  
  stu.num = length(Z.list)
  n.zero.new = n.zero
  zero.acc.new = zero.acc
  
  for(k in start.nperm:end.nperm){
    
    Z.perm.list = Z.list
    for (i in 1:stu.num){
      idx = sample(1:nrow(Z.list[[i]]))
      Z.perm.list[[i]][,Z.par.index] =  Z.perm.list[[i]][idx,Z.par.index]
    }
    if(use.cpp){ # if use.cpp = T, use Rcpp function to calculate this value
      score.stat.alpha.perm = try( .Score.test.stat.zero.meta.4Gresampling.c(Z.perm.list, Z.par.index,n.par.interest.alpha, col.zero.index.list, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, Method = Method) )
    }
    else{
      score.stat.alpha.perm = try( .Score.test.stat.zero.meta.4Gresampling(Z.perm.list, Z.par.index,n.par.interest.alpha, col.zero.index.list, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, Method = Method) )
    }
    
    if(!("try-error" %in% class(score.stat.alpha.perm))){# if score.stat.alpha.perm exists, then n.one.new + 1
      # if the permutation test statistics greater than original test statistics, cnt + 1
      n.zero.new = n.zero.new + 1
      if(score.stat.alpha.perm >= score.stat.zero.meta){
        zero.acc.new = zero.acc.new + 1 # if score.stat.alpha.perm >= score.stat.zero.meta one.acc.new + 1
        
      }
    }
    
    
  }
  
  # adaptive adjust the total number of iterations according to the number of resampling statistics which are more extreme than the original one in each loop
  # if the total number of permutation results which are greater than original results too small, enlarge the
  # number of total iterations
  if(zero.acc.new < 1){
    next.end.nperm = (end.nperm + 1) * 100 - 1;
    flag = 1;
    
  }else if(zero.acc.new<10){
    next.end.nperm = ( end.nperm + 1) * 10 - 1;
    flag = 1;
    
  }
  #   else if(one.acc.new<20){
  #     next.end.nperm = ( end.nperm + 1) * 5 - 1;
  #     flag = 1;
  #
  #   }
  else{
    next.end.nperm = ( end.nperm + 1) - 1;
    flag = 0;
  }
  
  return(list(n.zero.new=n.zero.new, zero.acc.new=zero.acc.new,
              flag=flag, next.end.nperm=next.end.nperm))
  
}

# Z.list: a list of covariates for zero part: first column is always intercept
# Z.par.index: index for the parameter of interest for the Z part

.Score.test.zero.meta <- function(Y.list, Z.list, Z.par.index, seed=11, resample=FALSE, n.replicates=NULL, use.cpp = F, Method = "FE-MV*"){
  stu.num = length(Z.list)
  p.par = length(Z.par.index)
  m = ncol(Y.list[[1]])
  p.zero = ncol(Z.list[[1]])
  Y0.list = Y.list
  # the total number of parameter of interest
  n.par.interest.alpha = m*length(Z.par.index)
  remove.study = NULL
  for(j in 1:stu.num)
  {
    # for each study set those value which are zero to 1
    Y0.list[[j]][Y.list[[j]]==0] = 1
    # for each study set those value which are non-zero to 0
    Y0.list[[j]][Y.list[[j]]>0] = 0
    remove.rows.idx = which(rowSums(Y0.list[[j]]) == 0)
    if(length(remove.rows.idx) == nrow(Y0.list[[j]])){
      remove.study = append(remove.study, i)
      next
    }
    # if(length(remove.rows.idx)>0){
    #   Y.list[[i]] <- Y.list[[i]][-remove.rows.idx, ,drop=FALSE]
    #   X.list[[i]] <- X.list[[i]][-remove.rows.idx, , drop=FALSE]
    # }
  }
  if(!is.null(remove.study)){
    if(length(remove.study) == stu.num){
      stop("Error: Not proper data")
    }else{
      Y0.list <- Y0.list[-remove.study]
      Z.list <- Z.list[-remove.study]
      stu.num = length(Z.list)
    }
  }
  if(! Method %in% c("FE-MV*","FE-VC*","RE-MV*","RE-VC*")){
    stop("Error: Please Choose a Proper Meta-analysis Method")
  }
  # check the parameters of interest to ensure the intercept term is not in it
  if (length(unique(sapply(1:stu.num,function(j) ncol(Y.list[[j]]))))!=1){
    stop("Error: The taxon in each study should be the same")
  }
  if(m<=1){
    stop("Error: Improper dimension for OTU table")
  }
  if(is.null(Z.par.index)){
    stop("Error: Testing parameters for the intercept is not informative. (Alpha part)")
  }
  ava.cnt = 0
  # col.pos.index.lst = lapply(1:stu.num, function(j) )
  ############################# Asymptotic: zero part
  
  
  # if all 0 in one group across across taxa, then output NA
  #if( (ncol(Y0)-length(remove.index))<=1 | sum(Y0[case==1,])==0 | sum(Y0[case==0,])==0){
  
  # initialize for later use
  score.pvalue = NA
  score.stat.meta = NA
  df = NA
  score.stat.alpha = NULL
  remove.index = NULL
  score.alpha = list()
  est.cov.zero = list()
  zero.vA.list.meta = list()
  zero.Vinv.list.meta = list()
  zero.VY.list.meta = list()
  col.zero.index.list = list()
  # initialize the score statistics and estimate covariance matrix for meta analysis
  score.alpha.meta  = rep(0,n.par.interest.alpha) ## A
  est.cov.meta = matrix(0, nrow = n.par.interest.alpha, ncol = n.par.interest.alpha) ## B
  for(i in 1:stu.num){
    Y0 = Y0.list[[i]]
    Z = Z.list[[i]]
    col.zero.index = which(apply(Y0, 2, function(x) length(table(x)) ) > 1) # keep those taxa which have both zero and one (both positive and zero observations)
    Y0 = Y0[, col.zero.index , drop=FALSE] # only save those columns which have both 0 and 1 values
    if(length(col.zero.index)<1)
    {
      remove.index = append(remove.index,i)
      next
    }else{
      # nY0 = rowSums(Y0)
      # ## remove 03/28/2016
      # index.subj.zero = which(nY0>0)
      # if(length(index.subj.zero)== 0){
      #     next
      #}
      #
      # Y0 = Y0[index.subj.zero, , drop=FALSE]
      # Z = Z[index.subj.zero, , drop=FALSE]
      #
      tmp.zero = try( .Score.test.stat.zero(Y0, Z, Z.par.index, "independence") )
    }
    if("try-error" %in% class(tmp.zero)){
      remove.index= append(remove.index,i)
      next
      
    }else{
      # number of study which can get score statistics
      ava.cnt = ava.cnt + 1
      # get the index for parameter of interest after score statistics and estimate covariance matrix are reorginzed
      # different across studies because of different column numbers
      idx = kronecker((col.zero.index-1)*p.par, rep(1,p.par)) + c(1:p.par)
      col.zero.index.list[[ava.cnt]] = idx  # the index of alpha parameter of interest in this study
      # save these outcome in list form for later meta-analysis as well as resampling test
      score.stat.alpha = append(score.stat.alpha, tmp.zero$score.stat.alpha)
      score.alpha[[ava.cnt]] = tmp.zero$score.alpha
      est.cov.zero[[ava.cnt]] = tmp.zero$est.cov.zero
      zero.vA.list.meta[[ava.cnt]] = tmp.zero$vA.list
      zero.Vinv.list.meta[[ava.cnt]] = tmp.zero$Vinv.list
      zero.VY.list.meta[[ava.cnt]] = tmp.zero$VY.list
      score.alpha.meta[idx] =  score.alpha.meta[idx] + crossprod(ginv(tmp.zero$est.cov.zero), tmp.zero$score.alpha) # add according to the index of parameter of interest for each study
      est.cov.meta[idx, idx] =  est.cov.meta[idx, idx] + ginv(tmp.zero$est.cov.zero)
    }
  }
  
  ############################# Asymptotic: combined
  if(ava.cnt>0){
    if(length(remove.index) != 0){
      Z.list = Z.list[-remove.index]
      Y0.list = Y0.list[-remove.index]
    }
    # save the index of those elements that have values greater than zero in score.aplha.meta vector
    save.index.zero = which(abs(score.alpha.meta) >= 1e-7)
    n.par.save.alpha = length(save.index.zero)
    score.alpha.meta = score.alpha.meta[save.index.zero]
    est.cov.meta = est.cov.meta[save.index.zero,save.index.zero]
    est.cov.inv = ginv(est.cov.meta)
    if(Method == "FE-MV*"){
      score.stat.meta = crossprod( score.alpha.meta,crossprod(t(est.cov.inv), score.alpha.meta))
      score.pvalue = 1- pchisq(score.stat.meta,df = n.par.save.alpha )
      df = n.par.save.alpha
    }
    if(Method == "FE-VC*"){
      weight.cov.meta = eigen(est.cov.meta)$values
      score.stat.meta = crossprod(score.alpha.meta)
      score.pvalue = davies(score.stat.meta, weight.cov.meta, h = rep(1,n.par.save.alpha), delta = rep(0,n.par.save.alpha), sigma = 0, lim = 10000, acc = 0.0001)$Qq
      score.pvalue = ifelse(score.pvalue>0,score.pvalue,0)
      df = n.par.save.alpha
    }
    if(Method == "RE-MV*"){
      U.theta = 0
      V.theta = 0
      for( i in 1:ava.cnt){
        est.inv = ginv(est.cov.zero[[i]])
        U.theta = U.theta + 1/2 * crossprod(crossprod(t(est.inv), score.alpha[[i]])) - 1/2 * tr(est.inv)
        V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
      }
      score.stat.meta = crossprod( score.alpha.meta,crossprod(t(est.cov.inv), score.alpha.meta)) + U.theta^2/V.theta
    }
    if(Method == "RE-VC*"){
      U.theta = 0
      V.theta = 0
      U.tau = 0
      est.inv.sum = matrix(0, nrow = n.par.interest.alpha, ncol = n.par.interest.alpha)
      for( i in 1:ava.cnt){
        idx = col.zero.index.list[[i]]
        est.inv = ginv(est.cov.zero[[i]])# calculate the generalized inverse for each estimate covariance for each study
        U.theta = U.theta +  1/2 * crossprod(crossprod(t(est.inv), score.alpha[[i]])) - 1/2 * tr(est.inv)
        # when \tau matrix and W matrix is identity the elements in upper left, bottom right as well as bottom left are the same in RE-VC test
        V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
        est.inv.sum[idx,idx] = est.inv.sum[idx,idx] + est.inv
      }
      V.tau = 1/2 * tr(crossprod(est.inv.sum))
      U.tau = 1/2 * crossprod(score.alpha.meta) - 1/2 * tr(est.inv.sum)
      score.stat.meta = crossprod(c(U.tau, U.theta), crossprod(ginv(matrix(c(V.tau,rep(V.theta,3)),ncol = 2)), c(U.tau, U.theta)))
    }
    zero.results = list(score.stat = score.stat.meta, score.pvalue = score.pvalue, df = df)
  }
  # if resample = TRUE then will apply permutation method to get permuted p-value
  # adaptive resampling test
  if(resample){
    
    #print("simulated stat:")
    set.seed(seed)
    if(!is.na(score.stat.meta)){
      
      n.zero = 0
      zero.acc = 0
      
      start.nperm = 1;
      end.nperm = min(100,n.replicates);
      flag = 1
      while(flag & end.nperm <= n.replicates){
        
        
        results = .resample.work.zero.meta(Z.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, score.stat.meta, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, start.nperm, end.nperm, n.zero, zero.acc, use.cpp = use.cpp, Method = Method)
        n.zero = results$n.zero.new
        zero.acc = results$zero.acc.new
        flag = results$flag
        next.end.nperm = results$next.end.nperm
        
        if(flag){
          start.nperm = end.nperm + 1;
          end.nperm = next.end.nperm;
          
        }
        
        if(start.nperm < n.replicates & end.nperm > n.replicates){
          #warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings"))
          
          results = .resample.work.zero.meta(Z.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, score.stat.meta, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, start.nperm, end.nperm, n.zero, zero.acc, use.cpp = use.cpp, Method = Method)
          
          n.zero = results$n.zero.new
          zero.acc = results$zero.acc.new
          
        }
        
      }
      
      tmp = (zero.acc+1)/(n.zero+1) # to avoid n.one be zero # resampling p value
      
      
    }else{
      
      tmp = NA
    }
    
    zero.results = c(zero.results, score.Rpvalue = tmp)
    
  }
  return(zero.results)
}

## Cauchy Combine
ACAT<-function(Pvals,Weights=NULL){
  #### check if there is NA
  if (sum(is.na(Pvals))>0){
    stop("Cannot have NAs in the p-values!")
  }
  #### check if Pvals are between 0 and 1
  if ((sum(Pvals<0)+sum(Pvals>1))>0){
    stop("P-values must be between 0 and 1!")
  }
  #### check if there are pvals that are either exactly 0 or 1.
  is.zero<-(sum(Pvals==0)>=1)
  is.one<-(sum(Pvals==1)>=1)
  if (is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if (is.zero){
    return(0)
  }
  if (is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }
  
  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(Weights)){
    Weights<-rep(1/length(Pvals),length(Pvals))
  }else if (length(Weights)!=length(Pvals)){
    stop("The length of weights should be the same as that of the p-values")
  }else if (sum(Weights<0)>0){
    stop("All the weights must be positive!")
  }else{
    Weights<-Weights/sum(Weights)
  }
  
  
  #### check if there are very small non-zero p values
  is.small<-(Pvals<1e-16)
  if (sum(is.small)==0){
    cct.stat<-sum(Weights*tan((0.5-Pvals)*pi))
  }else{
    cct.stat<-sum((Weights[is.small]/Pvals[is.small])/pi)
    cct.stat<-cct.stat+sum(Weights[!is.small]*tan((0.5-Pvals[!is.small])*pi))
  }
  #### check if the test statistic is very large.
  if (cct.stat>1e+15){
    pval<-(1/cct.stat)/pi
  }else{
    pval<-1-pcauchy(cct.stat)
  }
  return(pval)
}


Rarefy <- function (otu.tab, depth = min(rowSums(otu.tab))){
  # Rarefaction function: downsample to equal depth
  #
  # Args:
  #		otu.tab: OTU count table, row - n sample, column - q OTU
  #		depth: required sequencing depth
  #
  # Returns:
  # 	otu.tab.rff: Rarefied OTU table
  #		discard: labels of discarded samples
  #
  otu.tab <- as.matrix(otu.tab)
  ind <- (rowSums(otu.tab) < depth)
  sam.discard <- rownames(otu.tab)[ind]
  otu.tab <- otu.tab[!ind, ]
  
  rarefy <- function(x, depth){
    y <- sample(rep(1:length(x), x), depth)
    y.tab <- table(y)
    z <- numeric(length(x))
    z[as.numeric(names(y.tab))] <- y.tab
    z
  }
  otu.tab.rff <- t(apply(otu.tab, 1, rarefy, depth))
  rownames(otu.tab.rff) <- rownames(otu.tab)
  colnames(otu.tab.rff) <- colnames(otu.tab)
  return(list(otu.tab.rff=otu.tab.rff, discard=sam.discard))
}


#' Title
#'
#' @param OTU a list of matrices containing OTU counts with each row corresponding to a sample and each column corresponding to an OTU or taxa. Each matrix's taxas are better to the same. The column name is mandatory.
#' @param X a list of matrices containing covariates for the positive-part test with each column pertaining to one variable (pertains to the covariate of interest or the confounders). The number of elements of X and OTU must be the same. The column number of each matrix in this list must be the same.
#' @param X.index a vector indicate the columns in X for the covariate(s) of interest.
#' @param Tax a matrix defines the taxonomy ranks with each row corresponding to an OTU or a taxa and each column corresponding to a rank (starting from the higher taxonomic level). Row name is mandatory and should be consistent with the column name of the OTU table, Column name should be formatted as "Rank1", "Rank2 ()"... etc
#'        If provided, tests will be performed for lineages based on the taxonomic rank. The output contains P-values for all lineages; a list of significant lineages controlling the false discovery rate (based on resampling p-value if resampling test was performed).
#'        If not provided, one test will be performed with all the OTUs and one p-value will be output.
#' @param Method Meta-analysis method to be used. Including fixed effect methods such as the FE-MV, FE-VC, FE-MV\*, FE-VC\* test and the Cauchy-combined p-value of all the fix effect tests denoted as FE-O test and random effect methods like RE-MV, RE-VC, RE-MV\*, RE-VC\* test and the Cauchy-combined p-value of all the random effect tests denoted as RE-O test.
#' @param min.depth keep samples with depths >= min.depth.
#' @param n.perm perform asymptotic test if n.perm is null, otherwise perform permutation tests using the specified number of resamplings.
#' @param fdr.alpha false discovery rate for multiple tests on the lineages.
#' @param use.cpp Logical value (default F). Whether to use Rcpp or not for resampling test.
#' @param rarefy Logical value (default F). Whether to perform Rarefy for the zero part test and the zero part p value of
#'
#' @return A list with this elements
#'    \item{lineage.pval}{p-values for all lineages. By default ( Method = "FE-MV", n.perm = NULL ), only the asymptotic test will be performed. If using random effect meta-analysis methods ( Method = "RE-VC" or Method = "RE-MV" ), then resampling test must be performed.}
#'    \item{sig.lineage}{a vector of significant lineages.}
#' @export
#'
#' @examples
#' data(data.meta)
#' OTU = data.meta$OTU
#' Tax = data.meta$Tax
#' case = data.meta$covariate
#' Memic(OTU, case, 1, Tax, Method = "FE-MV", min.depth=0, n.perm=NULL, fdr.alpha=0.05)
#' @import MASS
#' @import data.table
#' @import CompQuadForm
#' @import geepack
#' @import brglm2
#' @import psych
#' @importFrom stats coef optim pchisq
#' @importFrom dplyr bind_rows
#' @references
#' Tang ZZ, Chen G, Alekseyenko AV, Li H. (2017) A general framework for association analysis of microbial communities on a taxonomic tree.
#' \emph{Bioinformatics}
#' \doi{10.1093/bioinformatics/btw804}.
#' @references
#' Lee S, Teslovich TM, Boehnke M, Lin X. (2013) General framework for meta-analysis of rare variants in sequencing association studies. Am J Hum Genet.
#' \emph{Am J Hum Genet}
#' \doi{10.1016/j.ajhg.2013.05.010}.
#' @references
#' Benjamini, Yoav, and Yosef Hochberg.(1995) Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.
#' \emph{Journal of the Royal Statistical Society. Series B}
Memic <- function(OTU, X, X.index, Tax=NULL, Method = "FE-MV", min.depth=0, n.perm=NULL, use.cpp = F, fdr.alpha=0.05, rarefy = F){
  n.resample = n.perm
  n.OTU = length(OTU)
  n.X = length(X)
  # drop.col = NULL
  if(Method %in% c("RE-VC", "RE-MV", "RE-VC*", "RE-MV*", "RE-O")){
    if(is.null(n.perm)){
      stop("The p-value for random effect meta-analysis method must be got by resampling test")
    }
  }
  if(n.OTU != n.X)
  {
    stop("The study number of OTU table and Covariate should be the same")
  }
  # if(length(unique(sapply(1:n.OTU,function(j) ncol(OTU[[j]]))))!= 1){
  #   stop("The taxa in each study should be the same")
  # }
  remove.study = NULL
  for(i in 1:n.OTU)
  {
    if(!is.matrix(OTU[[i]])){
      OTU[[i]] = as.matrix(OTU[[i]])
    }
    
    if(!is.matrix(X[[i]])){
      X[[i]] = as.matrix(X[[i]])
    }
    
    if(nrow(OTU[[i]])!=nrow(X[[i]])){
      stop(paste0("Number of samples in the OTU table and the covariate table of study ", i,
                  " should be the same\n"))
    }
    remove.subject = which(rowSums(OTU[[i]])<=min.depth)
    if(length(remove.subject)>0){
      print(paste("Remove",length(remove.subject), "samples with read depth less or equal to", min.depth, "in OTU table", i,"\n"))
      X[[i]] = X[[i]][-remove.subject, ,drop=FALSE]
      OTU[[i]] = OTU[[i]][-remove.subject, ,drop=FALSE]
      if(length(remove.subject) == nrow(OTU[[i]])){
        remove.study = append(remove.study, i)
      }
    }
    # drop.col = union(drop.col,which(colSums(OTU[[i]])==0))
  }
  if(length(remove.study) == n.OTU){
    stop("Please provide proper data")
  }
  if(missing(X.index)){
    X.index = 1:ncol(X[[1]])
  }
  if(!is.null(remove.study)){
    X = X[-remove.study]
    OTU = OTU[-remove.study]
    n.OTU = length(OTU)
  }
  # join all the OTU table
  OTU.comb = as.data.frame(OTU[[1]])
  if(n.OTU>1){
    for (i in 2:n.OTU){
      OTU.comb = bind_rows(OTU.comb,as.data.frame(OTU[[i]]))
    }
  }
  OTU.comb[is.na(OTU.comb)] <- 0 # substitute NA to zero (do not matter because those cloumns will de delete during calculation)
  OTU.comb <- as.matrix(OTU.comb)
  if(!is.null(Tax)){
    # preserve those columns which have taxonomy information
    col.save = intersect(Tax$Rank6,colnames(OTU.comb))
    OTU.comb = OTU.comb[,colnames(OTU.comb) %in% col.save]
    Tax = Tax[Tax$Rank6 %in% col.save,]
  }
  # divide the combined OTU table according to the number of observations for each study
  batch.cnt <- unlist(lapply(OTU, function(x) nrow(x)))
  batch.cnt <- append(1,batch.cnt)
  batch.cnt <- cumsum(batch.cnt)
  count = list()
  for (i in 1:n.OTU) {
    count[[i]] <- OTU.comb[batch.cnt[i]:(batch.cnt[i+1]-1),]
  }
  if(rarefy){
    count.rare = lapply(count,function(X){Rarefy(X)$otu.tab.rff})
  }
  X = lapply(1:n.OTU,function(j) cbind(1, X[[j]])) # add the intercept term
  X.index = X.index + 1
  
  if(is.null(Tax)){ # perform one test using all OTUs
    if(is.null(n.perm)){
      if(Method %in% c('FE-MV', 'FE-VC')){
        tmp = try(.Score.test.meta(count, X, X.index, Method = Method))
        if(!("try-error" %in% class(tmp))){
          pval = c(tmp$score.pvalue)
          names(pval) = paste0("Asymptotic-",Method)
        }else{
          pval = NA
          names(pval) = paste0("Asymptotic-",Method)
        }
      }else if(Method %in% c('FE-MV*', 'FE-VC*')){
        if(rarefy){
          tmp = try(.Score.test.zero.meta(count.rare, X, X.index, Method = Method))
        }else{
          tmp = try(.Score.test.zero.meta(count, X, X.index, Method = Method))
        }
        if(!("try-error" %in% class(tmp))){
          pval = c( tmp$score.pvalue )
          names(pval) = paste0("Asymptotic-",Method)
        }else{
          pval = NA
          names(pval) = paste0("Asymptotic-",Method)
        }
      }else if(Method == 'FE-O'){
        Pval_c = c()
        for(i in c('FE-MV', 'FE-VC','FE-MV*', 'FE-VC*')){
          if(i %in% c('FE-MV', 'FE-VC')){
            tmp = try(.Score.test.meta(count, X, X.index, Method = i))
            if(!("try-error" %in% class(tmp))){
              Pval_c = append(Pval_c, tmp$score.pvalue)
            }
          }else{
            if(rarefy){
              tmp = try(.Score.test.zero.meta(count.rare, X, X.index, Method = i))
              if(!("try-error" %in% class(tmp))){
                Pval_c = append(Pval_c, tmp$score.pvalue)
              }
            }else{
              tmp = try(.Score.test.zero.meta(count, X, X.index, Method = i))
              if(!("try-error" %in% class(tmp))){
                Pval_c = append(Pval_c, tmp$score.pvalue)
              }
            }
          }
        }
        Pval_c = na.omit(Pval_c)
        if(length(Pval_c) == 0){
          pval = NA
          names(pval) = paste0("Asymptotic-",Method)
        }else{
          pval = ACAT(Pval_c)
          names(pval) = paste0("Asymptotic-",Method)
        }
      }
    }else{ # resampling test + asymptotic test
      # (Y.list, X.list, X.par.index, seed=11, resample=FALSE, n.replicates=NULL, Method = "FE-MV", Weight=NULL )
      if(Method %in% c('FE-MV', 'FE-VC', 'RE-MV', 'RE-VC')){
        tmp = try(.Score.test.meta(count, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = Method))
        if(!("try-error" %in% class(tmp))){
          if(Method %in% c("RE-VC", "RE-MV")){
            pval = c(tmp$score.Rpvalue)
            names(pval) = paste0("Resampling-",Method)
          }else{
            pval = c(tmp$score.pvalue, tmp$score.Rpvalue)
            names(pval) = c(paste0("Asymptotic-",Method),paste0("Resampling-",Method))
          }
        }else{
          if(Method %in% c("RE-VC", "RE-MV")){
            pval = NA
            names(pval) = paste0("Resampling-",Method)
          }else{
            pval = c(NA, NA)
            names(pval) = c(paste0("Asymptotic-",Method),paste0("Resampling-",Method))
          }
        }
      }else if(Method %in% c('FE-MV*', 'FE-VC*', 'RE-MV*', 'RE-VC*')){
        if(rarefy){
          tmp = try(.Score.test.zero.meta(count.rare, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = Method))
        }else{
          tmp = try(.Score.test.zero.meta(count, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = Method))
        }
        if(!("try-error" %in% class(tmp))){
          if(Method %in% c("RE-VC*", "RE-MV*")){
            pval = c(tmp$score.Rpvalue)
            names(pval) = paste0("Resampling-",Method)
          }else{
            pval = c(tmp$score.pvalue, tmp$score.Rpvalue)
            names(pval) = c(paste0("Asymptotic-",Method),paste0("Resampling-",Method))
          }
        }else{
          if(Method %in% c("RE-VC*", "RE-MV*")){
            pval = NA
            names(pval) = paste0("Resampling-",Method)
          }else{
            pval = c(NA, NA)
            names(pval) = c(paste0("Asymptotic-",Method),paste0("Resampling-",Method))
          }
        }
      }else if(Method == 'FE-O'){
        Pval_asy_c = c()
        Pval_res_c = c()
        for(i in c('FE-MV', 'FE-VC','FE-MV*', 'FE-VC*')){
          if(i %in% c('FE-MV', 'FE-VC')){
            tmp = try(.Score.test.meta(count, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = i))
            if(!("try-error" %in% class(tmp))){
              Pval_asy_c = append(Pval_asy_c, tmp$score.pvalue)
              Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
            }
          }else{
            if(rarefy){
              tmp = try(.Score.test.zero.meta(count.rare, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = i))
              if(!("try-error" %in% class(tmp))){
                Pval_asy_c = append(Pval_asy_c, tmp$score.pvalue)
                Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
              }
            }else{
              tmp = try(.Score.test.zero.meta(count, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = i))
              if(!("try-error" %in% class(tmp))){
                Pval_asy_c = append(Pval_asy_c, tmp$score.pvalue)
                Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
              }
            }
          }
        }
        Pval_asy_c = na.omit(Pval_asy_c)
        Pval_res_c = na.omit(Pval_res_c)
        if((length(Pval_asy_c) == 0)|(length(Pval_res_c) == 0)){
          pval = c(NA, NA)
          names(pval) = c(paste0("Asymptotic-",Method),paste0("Resampling-",Method))
        }else{
          pval = c(ACAT(Pval_asy_c),ACAT(Pval_res_c))
          names(pval) = c(paste0("Asymptotic-",Method),paste0("Resampling-",Method))
        }
      }else if(Method == 'RE-O'){
        Pval_res_c = c()
        for(i in c('RE-MV', 'RE-VC','RE-MV*', 'RE-VC*')){
          if(i %in% c('RE-MV', 'RE-VC')){
            tmp = try(.Score.test.meta(count, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = i))
            if(!("try-error" %in% class(tmp))){
              Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
            }
          }else{
            if(rarefy){
              tmp = try(.Score.test.zero.meta(count.rare, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = i))
              if(!("try-error" %in% class(tmp))){
                Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
              }
            }else{
              tmp = try(.Score.test.zero.meta(count, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = i))
              if(!("try-error" %in% class(tmp))){
                Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
              }
            }
          }
        }
        Pval_res_c = na.omit(Pval_res_c)
        if(length(Pval_res_c) == 0){
          pval = c(NA)
          names(pval) = c(paste0("Resampling-",Method))
        }else{
          pval = ACAT(Pval_res_c)
          names(pval) = c(paste0("Resampling-",Method))
        }
      }
    }
    return( list(pval=pval) )
  }else{ # perform tests for lineages
    
    if(!is.matrix(Tax)){
      tax = as.matrix(Tax)
    }
    
    # if(length(drop.col)>0){
    #   tax = Tax[-drop.col,,drop = FALSE]
    # }else{
    #   tax = Tax
    # }
    for(i in 1:n.OTU)
    {
      if( sum(!(colnames(count[[i]]) %in% rownames(tax)))>0 ){
        stop(paste0("Error: OTU IDs in OTU table ",i," are not consistent with OTU IDs in Tax table"))
      }
    }
    
    n.rank = ncol(tax)
    # merge the tax and count together for later partition and merge OTU table according to taxonomy information
    W.data.list = lapply(1:n.OTU,function(j) data.table(data.frame(tax, t(count[[j]]))))
    otucols = lapply(1:n.OTU,function(j) names(W.data.list[[j]])[-(1:n.rank)])
    if(rarefy){
      W.data.rare.list = lapply(1:n.OTU,function(j) data.table(data.frame(tax, t(count.rare[[j]]))))
    }
    n.level = n.rank-1
    
    subtree = NULL
    pval = NULL
    
    for(k in 1:n.level){
      
      Rank.low = paste("Rank", n.rank-k,sep="")
      Rank.high = paste("Rank", n.rank-k+1,sep="")
      
      tmp = table(tax[,n.rank-k])
      level.uni = sort( names(tmp)[which(tmp>1)] )
      m.level = length(level.uni)
      # partition and merge OTU table according to taxonomy information
      tt = lapply(1:n.OTU, function(j) W.data.list[[j]][, lapply(.SD , sum, na.rm=TRUE), .SDcols=as.vector(unlist(otucols[j])), by=list( get(Rank.low), get(Rank.high) )])
      tt = lapply(1:n.OTU,function(j) setnames(tt[[j]], 1:2, c(Rank.low, Rank.high)))
      W.tax = as.vector(unlist(tt[[1]][, Rank.low, with=FALSE]))
      W.count = lapply(1:n.OTU,function(j) tt[[j]][, otucols[[j]], with=FALSE])
      if(rarefy){
        # partition and merge OTU table according to taxonomy information
        tt.rare = lapply(1:n.OTU, function(j) W.data.rare.list[[j]][, lapply(.SD , sum, na.rm=TRUE), .SDcols=as.vector(unlist(otucols[j])), by=list( get(Rank.low), get(Rank.high) )])
        tt.rare = lapply(1:n.OTU,function(j) setnames(tt[[j]], 1:2, c(Rank.low, Rank.high)))
        W.rare.tax = as.vector(unlist(tt.rare[[1]][, Rank.low, with=FALSE]))
        W.rare.count = lapply(1:n.OTU,function(j) tt.rare[[j]][, otucols[[j]], with=FALSE])
      }
      
      for(j in 1:m.level){
        
        Y = lapply(1:n.OTU, function(i) t(W.count[[i]][which(W.tax == level.uni[j]), , drop=FALSE]))
        if(rarefy){
          Y.rare = lapply(1:n.OTU, function(i) t(W.rare.count[[i]][which(W.rare.tax == level.uni[j]), , drop=FALSE]))
        }
        #Y = t(W.count[which(W.tax == "f__Veillonellaceae"), , drop=FALSE])
        
        # remove.index = NULL
        # for( i in 1:n.OTU)
        # {
        #   remove.index = union(remove.index,which(colSums(Y[[i]])==0))
        # }
        #
        # if(length(remove.index)==ncol(Y[[1]])){
        #
        #   #print("==skip:0==");
        #   next
        #
        #
        # }else{
        
        # if(length(remove.index)>0){
        #   Y = lapply(1:n.OTU, function(i) Y[[i]][, -remove.index, drop=FALSE])
        # }
        
        
        if(ncol(Y[[1]])==1){
          
          next
          #print("==skip:1==");
          
        }else{
          
          subtree = c(subtree, level.uni[j])
          #print(level.uni[j])
          if(is.null(n.resample)){ # asymptotic test only
            # (Y.list, X.list, X.par.index, seed=11, resample=FALSE, n.replicates=NULL, Method = "FE-MV", Weight=NULL )
            #  run test for each lineage
            if(Method %in% c('FE-MV', 'FE-VC')){
              tmp = try(.Score.test.meta(Y, X, X.index, Method = Method))
              if(!("try-error" %in% class(tmp))){
                pval = cbind(pval, c(tmp$score.pvalue))
              }else{
                pval = cbind(pval, NA)
              }
            }else if(Method %in% c('FE-MV*', 'FE-VC*')){
              if(rarefy){
                tmp = try(.Score.test.zero.meta(Y.rare, X, X.index, Method = Method))
              }else{
                tmp = try(.Score.test.zero.meta(Y, X, X.index, Method = Method))
              }
              if(!("try-error" %in% class(tmp))){
                pval = cbind(pval, tmp$score.pvalue)
              }else{
                pval = cbind(pval, NA)
              }
            }else if(Method == 'FE-O'){
              Pval_c = c()
              for(i in c('FE-MV', 'FE-VC','FE-MV*', 'FE-VC*')){
                if(i %in% c('FE-MV', 'FE-VC')){
                  tmp = try(.Score.test.meta(Y, X, X.index, Method = i))
                  if(!("try-error" %in% class(tmp))){
                    Pval_c = append(Pval_c, tmp$score.pvalue)
                  }
                }else{
                  if(rarefy){
                    tmp = try(.Score.test.zero.meta(Y.rare, X, X.index, Method = i))
                    if(!("try-error" %in% class(tmp))){
                      Pval_c = append(Pval_c, tmp$score.pvalue)
                    }
                  }else{
                    tmp = try(.Score.test.zero.meta(Y, X, X.index, Method = i))
                    if(!("try-error" %in% class(tmp))){
                      Pval_c = append(Pval_c, tmp$score.pvalue)
                    }
                  }
                }
              }
              Pval_c = na.omit(Pval_c)
              if(length(Pval_c) == 0){
                pval = cbind(pval, NA)
              }else{
                pval = cbind(pval, ACAT(Pval_c))
              }
            }
          }
          else{
            # if n.resample in not null, select the significant lineage according to the resampling pvalue
            #  run test for each lineage
            if(Method %in% c('FE-MV', 'FE-VC', 'RE-MV', 'RE-VC')){
              tmp = try(.Score.test.meta(Y, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = Method))
              if(!("try-error" %in% class(tmp))){
                if(Method %in% c("RE-VC", "RE-MV")){
                  pval = cbind(pval, tmp$score.Rpvalue)
                }else{
                  pval = cbind(pval, c(tmp$score.pvalue, tmp$score.Rpvalue) )
                }
              }else{
                if(Method %in% c("RE-VC", "RE-MV")){
                  pval = cbind(pval, NA)
                }else{
                  pval = cbind(pval, c(NA, NA) )
                }
              }
            }else if(Method %in% c('FE-MV*', 'FE-VC*', 'RE-MV*', 'RE-VC*')){
              if(rarefy){
                tmp = try(.Score.test.zero.meta(Y.rare, X, X.index, seed=11, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = Method))
              }else{
                tmp = try(.Score.test.zero.meta(Y, X, X.index, seed=11, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = Method))
              }
              if(!("try-error" %in% class(tmp))){
                if(Method %in% c("RE-VC*", "RE-MV*")){
                  pval = cbind(pval, tmp$score.Rpvalue)
                }else{
                  pval = cbind(pval, c(tmp$score.pvalue, tmp$score.Rpvalue) )
                }
              }else{
                if(Method %in% c("RE-VC*", "RE-MV*")){
                  pval = cbind(pval, NA)
                }else{
                  pval = cbind(pval, c(NA, NA) )
                }
              }
            }else if(Method == 'FE-O'){
              Pval_asy_c = c()
              Pval_res_c = c()
              for(i in c('FE-MV', 'FE-VC','FE-MV*', 'FE-VC*')){
                if(i %in% c('FE-MV', 'FE-VC')){
                  tmp = try(.Score.test.meta(Y, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = i))
                  if(!("try-error" %in% class(tmp))){
                    Pval_asy_c = append(Pval_asy_c, tmp$score.pvalue)
                    Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
                  }
                }else{
                  if(rarefy){
                    tmp = try(.Score.test.zero.meta(Y.rare, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = i))
                    if(!("try-error" %in% class(tmp))){
                      Pval_asy_c = append(Pval_asy_c, tmp$score.pvalue)
                      Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
                    }
                  }else{
                    tmp = try(.Score.test.zero.meta(Y, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = i))
                    if(!("try-error" %in% class(tmp))){
                      Pval_asy_c = append(Pval_asy_c, tmp$score.pvalue)
                      Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
                    }
                  }
                }
              }
              Pval_asy_c = na.omit(Pval_asy_c)
              Pval_res_c = na.omit(Pval_res_c)
              if((length(Pval_asy_c) == 0)|(length(Pval_res_c) == 0)){
                pval = cbind(pval, c(NA, NA))
              }else{
                pval = cbind(pval, c(ACAT(Pval_asy_c),ACAT(Pval_res_c)))
              }
            }else if(Method == 'RE-O'){
              Pval_res_c = c()
              for(i in c('RE-MV', 'RE-VC','RE-MV*', 'RE-VC*')){
                if(i %in% c('RE-MV', 'RE-VC')){
                  tmp = try(.Score.test.meta(Y, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = i))
                  if(!("try-error" %in% class(tmp))){
                    Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
                  }
                }else{
                  if(rarefy){
                    tmp = try(.Score.test.zero.meta(Y.rare, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = i))
                    if(!("try-error" %in% class(tmp))){
                      Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
                    }
                  }else{
                    tmp = try(.Score.test.zero.meta(Y, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = i))
                    if(!("try-error" %in% class(tmp))){
                      Pval_res_c = append(Pval_res_c, tmp$score.Rpvalue)
                    }
                  }
                }
              }
              Pval_res_c = na.omit(Pval_res_c)
              if(length(Pval_res_c) == 0){
                pval = cbind(pval, NA)
              }else{
                pval = cbind(pval, ACAT(Pval_res_c))
              }
            }
          }
          
        }
        
      } # lineage loop
    } # level loop
    
    colnames(pval) = subtree
    if(is.null(n.resample)){
      rownames(pval) = paste0("Asymptotic-",Method)
      score.tmp = pval[1,]
    }else{
      if(Method %in% c("RE-VC", "RE-MV","RE-VC*", "RE-MV*", "RE-O")){
        rownames(pval) = paste0("Resampling-",Method)
        score.tmp = pval[1,]
      }else if(Method %in% c('FE-MV', 'FE-VC','FE-MV*', 'FE-VC*', "FE-O")){
        rownames(pval) = c(paste0("Asymptotic-",Method),paste0("Resampling-",Method))
        score.tmp = pval[2,]
      }
      #print(pval)
    }
    
    # identify significant lineages
    subtree.tmp = subtree
    index.na = which(is.na(score.tmp))
    if(length(index.na)>0){
      # drop those lineages which have NA values
      score.tmp = score.tmp[-index.na]
      subtree.tmp = subtree.tmp[-index.na]
    }
    
    #score.tmp[score.tmp==0] = 1e-4
    m.test = length(score.tmp)
    
    # Benjamini-Hochberg FDR control
    index.p = order(score.tmp)
    p.sort = sort(score.tmp)
    #fdr.alpha = 0.05
    
    # change 2022/12/01
    reject = rep(0, m.test)
    tmp = which(p.sort<=(1:m.test)*fdr.alpha/m.test)
    if(length(tmp)>0){
      index.reject = index.p[1:max(tmp)]
      reject[index.reject] = 1
    }
    
    sig.lineage = subtree.tmp[reject==1]
    
    # return all the p-values as well as significant lineages
    return( list(lineage.pval=pval, sig.lineage=sig.lineage) )
  }
  
}

