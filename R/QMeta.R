# one part model

.F.test <- function(x) {
  # Fisher's p-value combination
  x.stat <- -2 * sum(log(x))
  return(1 - pchisq(x.stat, df = 2 * length(x)))
}

.simes.test <- function(x) {
  return(min(length(x) * x / rank(x)))
}


########################################
#                                      #
#           One Part Model             #
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

  return(-loglik)
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
    est.reduce.beta[-par.interest.index.beta] <- optim(par = beta.ini.reduce, fn = .fun.neg.loglik.beta, gr = .fun.neg.score.beta, data = data.reduce.beta, method = "BFGS")$par


    data.beta <- list(Y = Y, X = X)
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

    B1 <- cbind(diag(n.par.interest.beta), tcrossprod(-Hess.reduce.reorg[(1:n.par.interest.beta), ((n.par.interest.beta + 1):n.beta)], t(ginv(Hess.reduce.reorg[((n.par.interest.beta + 1):n.beta), ((n.par.interest.beta + 1):n.beta)]))))

    B2 <- matrix(0, n.beta, n.beta)

    for (i in 1:n) {
      B2 <- B2 + Score.reduce.reorg[i, ] %o% Score.reduce.reorg[i, ]
    }

    B = crossprod(t(B1), tcrossprod(B2, B1))
    score.stat.beta = crossprod( A, crossprod(t(ginv(B)), A))
  }

  # save those summary statistics for later use
  return(list(score.stat.beta = score.stat.beta, score.beta = A, est.cov = B, S.beta.list = S.beta.list, I.beta.list = I.beta.list))
}

.Score.test.stat.meta.4Gresampling <- function(X.perm.list, X.par.index, n.par.interest.beta, col.index.list, S.beta.list.meta, I.beta.list.meta, Method = "FE-MV") {
  stu.num <- length(X.perm.list) #  the total number of studies for meta analysis

  # initialize those statistics
  score.stat.beta <- NULL
  score.beta <- NULL
  est.cov <- NULL
  score.pvalue.beta <- NULL

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
    A = colSums(Score.reduce.beta.perm.reorg)[1:n.par.beta.interest]

    B1 = cbind(diag(n.par.beta.interest), tcrossprod(-Hess.reduce.beta.perm.reorg[(1:n.par.beta.interest), ((n.par.beta.interest+1):n.beta)] , t(ginv(Hess.reduce.beta.perm.reorg[((n.par.beta.interest+1):n.beta), ((n.par.beta.interest+1):n.beta)]))))

    B2 =  matrix(0, n.beta, n.beta)

    for (i in 1:n) {
      B2 <- B2 + Score.reduce.beta.perm.reorg[i, ] %o% Score.reduce.beta.perm.reorg[i, ]
    }

    B <-  crossprod(t(B1), tcrossprod(B2, B1))

    # Generate the test statistics for each study
    score.stat.beta.perm <- crossprod( A, crossprod(t(ginv(B)), A))
    score.stat.beta <- append(score.stat.beta, score.stat.beta.perm)
    score.beta[[j]] <- A
    est.cov[[j]] <- B
    score.pvalue.beta <- append(score.pvalue.beta, (1 - pchisq(score.stat.beta.perm, n.par.interest.beta)))
  }

  # initialize the score statistics and estimate covariance matrix for meta-analysis
  score.beta.meta <- rep(0, n.par.interest.beta)
  est.cov.meta <- matrix(0, nrow = n.par.interest.beta, ncol = n.par.interest.beta)
  for (i in 1:stu.num)
  {
    idx <- col.index.list[[i]] # specify the interested parameter index for each study
    score.beta.meta[idx] <- score.beta.meta[idx] + score.beta[[i]]
    est.cov.meta[idx, idx] <- est.cov.meta[idx, idx] + est.cov[[i]]
  }
  est.cov.inv <- ginv(est.cov.meta)
  if (Method == "FE-MV") {
    score.stat.meta.perm <- crossprod( score.beta.meta, crossprod(t(est.cov.inv), score.beta.meta)) # Multivariate test
  }
  # if (Method == "SKAT") {
  #   W <- diag(1, nrow = n.par.interest.beta)
  #   score.stat.meta.perm <- score.beta.meta %*% W %*% score.beta.meta  # FE-SKAT-test(FE-SKAT)
  # }
  if (Method == "FE-VC") {
    weight.cov.inv <- eigen(est.cov.inv)$values
    score.stat.meta.perm <- crossprod( score.beta.meta, crossprod(t(est.cov.inv), crossprod( t(est.cov.inv), score.beta.meta))) # FE-VC-test(SKAT-VC)
  }
  # if (Method == "Fisher") {
  #   score.stat.meta.perm <- -2 * sum(log(score.pvalue.beta)) # the Fisher's p-value combination
  # }
  if(Method == "Het-SKAT"){
    #W = diag(1,nrow = n.par.interest.beta)
    score.stat.meta.perm = 0
    for( i in 1:stu.num ){
      est.inv = ginv(est.cov[[i]])
      score.stat.meta.perm = score.stat.meta.perm + tcrossprod(crossprod(score.beta[[i]], est.inv), crossprod(score.beta[[i]], est.inv))
    }
  }
  if(Method == "RE-SKAT"){
    U.tau.b = 0
    a2 = 0
    a3 = 0
    a4 = 0
    for( i in 1:stu.num ){
      est.inv = ginv(est.cov[[i]])
      U.tau.b = U.tau.b + tcrossprod(crossprod(score.beta[[i]],est.inv), crossprod(score.beta[[i]], est.inv))
      a2 = a2 + sum(diag(crossprod(t(est.inv),est.inv)))
      a3 = a3 + sum(diag(crossprod(t(est.inv),est.inv)))
      a4  = a4 + sum(diag(crossprod(t(est.inv),est.inv)))
    }
    a1 = sum(diag(crossprod(t(est.cov.inv),est.cov.inv)))
    U.tau.b = 1/2 * U.tau.b + 1/2 * sum(diag(est.cov.inv))
    U.tau.w = 1/2 * crossprod( score.beta.meta, crossprod(t(est.cov.inv), crossprod( t(est.cov.inv), score.beta.meta)))
    + 1/2 * sum(diag(est.cov.inv))
    score.stat.meta.perm = crossprod(c(U.tau.w, U.tau.b), crossprod(1/2 * matrix(c(a1,a2,a3,a4),ncol = 2), c(U.tau.w, U.tau.b)))
  }

  return(as.numeric(score.stat.meta.perm))  # return the test statistics
}


.Score.test.stat.meta.4Gresampling.c <- function(X.perm.list, X.par.index, n.par.interest.beta, col.index.list, S.beta.list.meta, I.beta.list.meta, Method = "FE-MV"){
  tmp = score_test_stat_meta_resampling_c(X.perm.list, col.index.list, S.beta.list.meta, I.beta.list.meta, X.par.index, n.par.interest.beta)
  est.cov.meta = tmp$est_cov_meta
  score.beta.meta = tmp$score_beta_meta
  est.cov = tmp$est_cov
  score.beta = tmp$score_beta
  save.index.pos = which(abs(score.beta.meta) >= 1e-7)
  score.beta.meta = score.beta.meta[save.index.pos]
  est.cov.meta = est.cov.meta[save.index.pos,save.index.pos]
  est.cov.inv = ginv(est.cov.meta)
  stu.num = length(X.perm.list)
  if (Method == "FE-MV")
  {
    score.stat.meta.perm = crossprod( score.beta.meta, crossprod(t(est.cov.inv), score.beta.meta)) #FE-Burden
  }
  # if (Method == "SKAT")
  # {
  #     W = diag(1,nrow = n.par.interest.beta)
  #   #fesk.p = farebrother(score.stat.fesk,weight, h = rep(1,m-1), delta = rep(0,m-1), maxit = 100000,eps = 10^(-10), mode = 1)$Qq
  #   score.stat.meta.perm = crossprod( score.beta.meta, crossprod(t(W), score.beta.meta))
  # }
  if (Method == "FE-VC")
  {
    weight.cov.inv = eigen(est.cov.inv)$values #eign.val/sum(eign.val)
    score.stat.meta.perm = crossprod( score.beta.meta,  crossprod(t(est.cov.inv), crossprod(t(est.cov.inv), score.beta.meta))) #SKAT-VC
  }
  # if (Method == "Fisher")
  # {
  #   score.stat.meta.perm = -2 * sum(log(score.pvalue.beta))
  # }
  #
  if(Method == "Het-SKAT"){
    #W = diag(1,nrow = n.par.interest.beta)
    score.stat.meta.perm = 0

    for( i in 1:stu.num ){
      est.inv = ginv(est.cov[[i]])
      score.stat.meta.perm = score.stat.meta.perm + tcrossprod(crossprod(score.beta[[i]], est.inv), crossprod(score.beta[[i]], est.inv))
    }
  }
  if(Method == "RE-SKAT"){
    U.tau.b = 0
    a2 = 0
    a3 = 0
    a4 = 0
    for( i in 1:stu.num ){
      est.inv = ginv(est.cov[[i]])
      U.tau.b = U.tau.b + tcrossprod(crossprod(score.beta[[i]],est.inv), crossprod(score.beta[[i]], est.inv))
      a2 = a2 + sum(diag(crossprod(t(est.inv),est.inv)))
      a3 = a3 + sum(diag(crossprod(t(est.inv),est.inv)))
      a4  = a4 + sum(diag(crossprod(t(est.inv),est.inv)))
    }
    a1 = sum(diag(crossprod(t(est.cov.inv),est.cov.inv)))
    U.tau.b = 1/2 * U.tau.b + 1/2 * sum(diag(est.cov.inv))
    U.tau.w = 1/2 * crossprod( score.beta.meta, crossprod(t(est.cov.inv), crossprod( t(est.cov.inv), score.beta.meta)))
    + 1/2 * sum(diag(est.cov.inv))  #SKAT-VC
    score.stat.meta.perm = crossprod(c(U.tau.w, U.tau.b), crossprod(1/2 * matrix(c(a1,a2,a3,a4),ncol = 2), c(U.tau.w, U.tau.b)))
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
  n = nrow(Y.list[[1]])
  p = ncol(X.list[[1]])
  n.OTU = length(Y.list)
  n.par.interest.beta = (m-1)*length(X.par.index)

  if(sum(X.par.index == 1)){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }
  # check the parameters of interest to ensure the intercept term is not in it
  if(! Method %in% c("FE-MV","FE-VC","Het-SKAT","RE-SKAT")){
    stop("Error: Please Choose a Proper Meta-analysis Method")
  }
  # the taxa of each study should be the same
  if (length(unique(sapply(1:n.OTU,function(j) ncol(Y.list[[j]]))))!=1){
    stop("Error: The taxon in each study should be the same")
  }
  # initialize the test statistics
  if(is.null(X.par.index) | n==0){

    score.pvalue = NA
    score.stat.meta = NA
    df = NA

  }else{
    score.stat.beta = NULL
    score.beta = NULL
    score.pvalue.beta = NULL
    score.pvalue = NA
    score.stat.meta = NA
    df = NA
    est.cov = list()
    S.beta.list.meta = list()
    I.beta.list.meta = list()
    col.index.list = list()
    remove.index = c()
    # initialize the score statistics and estimate covariance matrix for meta analysis
    score.beta.meta  = rep(0,n.par.interest.beta) ## A
    est.cov.meta = matrix(0, nrow = n.par.interest.beta, ncol = n.par.interest.beta) ## B
    j = 0
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
      col.index = col.index[-1]
      n.beta = (m - 1)*p
      tmp.one = try( .Score.test.stat(Y, X, X.par.index) )
      if( "try-error" %in% class(tmp.one) ){
        score.stat.beta = score.stat.beta
        score.beta = score.beta
        est.cov = est.cov
        score.pvalue.beta = score.pvalue.beta
        S.beta.list.meta = S.beta.list.meta
        I.beta.list.meta = I.beta.list.meta
      }else{
        j = j+1
        idx = kronecker((col.index-2)*p.par, rep(1,p.par)) + c(1:p.par)
        col.index.list[[j]] = idx
        score.stat.beta = append(score.stat.beta, tmp.one$score.stat.beta)
        score.beta[[j]] = tmp.one$score.beta
        est.cov[[j]] = tmp.one$est.cov
        score.pvalue.beta = append(score.pvalue.beta, (1 - pchisq(tmp.one$score.stat.beta, n.par.interest.beta)))
        S.beta.list.meta[[j]] = tmp.one$S.beta.list
        I.beta.list.meta[[j]] = tmp.one$I.beta.list
        score.beta.meta[idx] =  score.beta.meta[idx] +  tmp.one$score.beta  #FE-Burden
        est.cov.meta[idx,idx] =  est.cov.meta[idx,idx] + tmp.one$est.cov #FE-SKAT
      }
    }

  }
  if(j>0){
    if(length(remove.index) != 0){
      X.list = X.list[-remove.index]
      Y.list = Y.list[-remove.index]
    }
    save.index.pos = which(abs(score.beta.meta) >= 1e-7)
    n.par.save.beta = length(save.index.pos)
    score.beta.meta = score.beta.meta[save.index.pos]
    est.cov.meta = est.cov.meta[save.index.pos,save.index.pos]
    est.cov.inv = ginv(est.cov.meta)
    if (Method == "FE-MV")
    {
      score.stat.meta = crossprod( score.beta.meta,crossprod(t(est.cov.inv), score.beta.meta)) #FE-Burden
      score.pvalue = 1- pchisq(score.stat.meta,df = n.par.save.beta )
      df = n.par.save.beta
    }
    # if (Method == "SKAT")
    # {
    #
    #   W = diag(1,nrow = n.par.interest.beta)
    #   #fesk.p = farebrother(score.stat.fesk,weight, h = rep(1,m-1), delta = rep(0,m-1), maxit = 100000,eps = 10^(-10), mode = 1)$Qq
    #   eigen.cov = eigen(est.cov.meta)
    #   eigen.cov.sqrt = crossprod(t(eigen.cov$vectors), crossprod(diag(sqrt(eigen.cov$values),nrow = length(eigen.cov$values)), solve(eigen.cov$vectors)))
    #   weight.cov = eigen(crossprod(t(eigen.cov.sqrt), crossprod(W, eigen.cov.sqrt)))$values #eign.val/sum(eign.val)
    #   score.stat.meta = crossprod(score.beta.meta, crossprod( W, score.beta.meta))
    #   score.pvalue = davies(score.stat.meta,weight.cov, h = rep(1,n.par.interest.beta), delta = rep(0,n.par.interest.beta), sigma = 0, lim = 10000, acc = 0.0001)$Qq
    #   score.pvalue = ifelse(score.pvalue>0,score.pvalue,1e-7)
    #   df = n.par.interest.beta
    # }
    if (Method == "FE-VC"){
      weight.cov.inv = eigen(est.cov.inv)$values #eign.val/sum(eign.val)
      score.stat.meta = crossprod( score.beta.meta, crossprod(t(est.cov.inv), crossprod( t(est.cov.inv), score.beta.meta))) #SKAT-VC
      score.pvalue = davies(score.stat.meta,weight.cov.inv, h = rep(1,n.par.save.beta), delta = rep(0,n.par.save.beta), sigma = 0, lim = 10000, acc = 0.0001)$Qq
      score.pvalue = ifelse(score.pvalue>0,score.pvalue,1e-7)
      df = n.par.save.beta
    }
    # if (Method == "Fisher")
    # {
    #   score.stat.meta = -2 * sum(log(score.pvalue.beta))
    #   score.pvalue = .F.test(score.pvalue.beta)
    #   df = length(score.pvalue.beta)
    # }
    if(Method == "Het-SKAT"){
      # W = diag(1,nrow = n.par.interest.beta)
      score.stat.meta = 0
      for( i in 1:j ){
        est.inv = ginv(est.cov[[i]])
        tmp  = tcrossprod(crossprod(score.beta[[i]],est.inv), crossprod(score.beta[[i]], est.inv))
        score.stat.meta = score.stat.meta + tcrossprod(crossprod(score.beta[[i]],est.inv), crossprod(score.beta[[i]], est.inv))
      }
    }
    if(Method == "RE-SKAT"){
      U.tau.b = 0
      a2 = 0
      a3 = 0
      a4 = 0
      for( i in 1:j ){
        est.inv = ginv(est.cov[[i]])
        U.tau.b = U.tau.b + tcrossprod(crossprod(score.beta[[i]],est.inv), crossprod(score.beta[[i]], est.inv))
        a2 = a2 + sum(diag(crossprod(t(est.inv),est.inv)))
        a3 = a3 + sum(diag(crossprod(t(est.inv),est.inv)))
        a4  = a4 + sum(diag(crossprod(t(est.inv),est.inv)))
      }
      a1 = sum(diag(crossprod(t(est.cov.inv),est.cov.inv)))
      U.tau.b = 1/2 * U.tau.b + 1/2 * sum(diag(est.cov.inv))
      U.tau.w = 1/2 * crossprod( score.beta.meta, crossprod(t(est.cov.inv), crossprod( t(est.cov.inv), score.beta.meta)))
      + 1/2 * sum(diag(est.cov.inv)) #SKAT-VC
      score.stat.meta = crossprod(c(U.tau.w, U.tau.b), crossprod(1/2 * matrix(c(a1,a2,a3,a4),ncol = 2), c(U.tau.w, U.tau.b)))
    }
  }
  beta.meta.results = list(score.stat = score.stat.meta, score.pvalue = score.pvalue, df = df)

  # if resample = TRUE then will apply permutation method to get permuted p-value
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

      tmp = (one.acc+1)/(n.one+1)

      #print(n.one)
      #print(one.acc)

    }else{

      tmp = NA
    }

    beta.meta.results = c(beta.meta.results, score.Rpvalue = tmp)


  }

  return(beta.meta.results)
}


#' Title
#'
#' @param OTU a list of matrices contains OTU counts with each row corresponding to a sample and each column corresponding to an OTU or taxa. Each matrix's taxas are better to the same. The column name is mandatory.
#' @param X a list of matrices containing covariates for the positive-part test with each column pertaining to one variable (pertains to the covariate of interest or the confounders). The number of elements of X and OTU must be the same. The column number of each matrix in this list must be the same.
#' @param X.index a vector indicate the columns in X for the covariate(s) of interest.
#' @param Tax a matrix defines the taxonomy ranks with each row corresponding to an OTU or a taxa and each column corresponding to a rank (starting from the higher taxonomic level). Row name is mandatory and should be consistent with the column name of the OTU table, Column name should be formatted as "Rank1", "Rank2 ()"... etc
#'        If provided, tests will be performed for lineages based on the taxonomic rank. The output contains P-values for all lineages; a list of significant lineages controlling the false discovery rate (based on resampling p-value if resampling test was performed); p-values of the global tests (Fisher- and Simes-combined the p-values for testing lineages).
#'        If not provided, one test will be performed with all the OTUs and one p-value will be output.
#' @param Method Meta-analysis method to be used. Including fixed effect methods such as the FE-MV test and FE-VC test and random effect methods like Het-SKAT and RE-SKAT test.
#' @param min.depth keep samples with depths >= min.depth.
#' @param n.perm perform asymptotic test if n.perm is null, otherwise perform permutation tests using the specified number of resamplings.
#' @param fdr.alpha false discovery rate for multiple tests on the lineages.
#' @param use.cpp Logical value (default F). Whether to use Rcpp or not for resampling test.
#'
#' @return A list with this elements
#'    \item{lineage.pval}{p-values for all lineages. By default ( Method = "FE-MV", n.perm = NULL ), only the asymptotic test will be performed. If using random effect meta-analysis methods ( Method = "RE-SKAT" or Method = "Het-SKAT" ), then resampling test must be performed.}
#'    \item{sig.lineage}{a vector of significant lineages.}
#' @export
#'
#' @examples
#' data(data.meta)
#' OTU = data.meta$OTU
#' Tax = data.meta$Tax
#' case = data.meta$covariate
#' QCAT_Meta(OTU, case, 1, Tax, Method = "FE-MV", min.depth=0, n.perm=NULL, fdr.alpha=0.05)
#' @import MASS
#' @import data.table
#' @import CompQuadForm
#' @import geepack
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
QCAT_Meta <- function(OTU, X, X.index, Tax=NULL, Method = "FE-MV", min.depth=0, n.perm=NULL, use.cpp = F, fdr.alpha=0.05){
  n.resample = n.perm
  n.OTU = length(OTU)
  n.X = length(X)
  # drop.col = NULL
  if(Method %in% c("RE-SKAT", "Het-SKAT")){
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
  for(i in 1:n.OTU)
  {
    if(!is.matrix(OTU[[i]])){
      OTU[[i]] = as.matrix(OTU[[i]])
    }

    if(!is.matrix(X[[i]])){
      X[[i]] = as.matrix(X[[i]])
    }

    if(nrow(OTU[[i]])!=nrow(X[[i]])){
      stop(paste0("Samples in the OTU table and the covariate table of study ", i,
                  " should be the same\n"))
    }
    remove.subject = which(rowSums(OTU[[i]])<min.depth)
    if(length(remove.subject)>0){
      print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth, "in OTU table", i,"\n"))
      X[[i]] = X[[i]][-remove.subject, ,drop=FALSE]
      OTU[[i]] = OTU[[i]][-remove.subject, ,drop=FALSE]
    }
    # drop.col = union(drop.col,which(colSums(OTU[[i]])==0))
  }

  if(missing(X.index)){
    X.index = 1:ncol(X[[1]])
  }

  OTU.comb = as.data.frame(OTU[[1]])
  for (i in 2:n.OTU){
    OTU.comb = bind_rows(OTU.comb,as.data.frame(OTU[[i]]))
  }
  OTU.comb[is.na(OTU.comb)] <- 0
  OTU.comb <- as.matrix(OTU.comb)
  if(!is.null(Tax)){
    col.save = intersect(Tax$Rank6,colnames(OTU.comb))
    OTU.comb = OTU.comb[,colnames(OTU.comb) %in% col.save]
    Tax = Tax[Tax$Rank6 %in% col.save,]
  }
  batch.cnt <- unlist(lapply(OTU, function(x) nrow(x)))
  batch.cnt <- append(1,batch.cnt)
  batch.cnt <- cumsum(batch.cnt)
  count = list()
  for (i in 1:n.OTU) {
    count[[i]] <- OTU.comb[batch.cnt[i]:(batch.cnt[i+1]-1),]
  }
  X = lapply(1:n.OTU,function(j) cbind(1, X[[j]])) # add the intercept term
  X.index = X.index + 1

  if(is.null(Tax)){ # perform one test using all OTUs
    if(is.null(n.perm)){
      pval = as.matrix(.Score.test.meta(count, X, X.index, Method = Method)$score.pvalue)
      colnames(pval) = paste0("Asymptotic-",Method)
    }else{ # resampling test + asymptotic test
      # (Y.list, X.list, X.par.index, seed=11, resample=FALSE, n.replicates=NULL, Method = "FE-MV", Weight=NULL )
      tmp = .Score.test.meta(count, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = Method)
      if(Method %in% c("RE-SKAT", "Het-SKAT")){
        pval = c(tmp$score.Rpvalue)
        names(pval) = paste0("Resampling-",Method)
      }else{
        pval = c(tmp$score.pvalue, tmp$score.Rpvalue)
        names(pval) = c(paste0("Asymptotic-",Method),paste0("Resampling-",Method))
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
    W.data.list = lapply(1:n.OTU,function(j) data.table(data.frame(tax, t(count[[j]]))))
    otucols = lapply(1:n.OTU,function(j) names(W.data.list[[j]])[-(1:n.rank)])
    n.level = n.rank-1

    subtree = NULL
    pval = NULL

    for(k in 1:n.level){

      Rank.low = paste("Rank", n.rank-k,sep="")
      Rank.high = paste("Rank", n.rank-k+1,sep="")

      tmp = table(tax[,n.rank-k])
      level.uni = sort( names(tmp)[which(tmp>1)] )
      m.level = length(level.uni)

      tt = lapply(1:n.OTU, function(j) W.data.list[[j]][, lapply(.SD , sum, na.rm=TRUE), .SDcols=as.vector(unlist(otucols[j])), by=list( get(Rank.low), get(Rank.high) )])
      tt = lapply(1:n.OTU,function(j) setnames(tt[[j]], 1:2, c(Rank.low, Rank.high)))
      W.tax = as.vector(unlist(tt[[1]][, Rank.low, with=FALSE]))
      W.count = lapply(1:n.OTU,function(j) tt[[j]][, otucols[[j]], with=FALSE])


      for(j in 1:m.level){

        Y = lapply(1:n.OTU, function(i) t(W.count[[i]][which(W.tax == level.uni[j]), , drop=FALSE]))

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
            tmp = .Score.test.meta(Y, X, X.index, Method = Method)
            pval = cbind(pval, c(tmp$score.pvalue))
          }
          else{
            # (Y.list, X.list, X.par.index, seed=11, resample=FALSE, n.replicates=NULL, Method = "FE-MV", Weight=NULL )
            if(Method %in% c("RE-SKAT", "Het-SKAT")){
              tmp = .Score.test.meta(Y, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = Method)
              pval = cbind(pval, c(tmp$score.Rpvalue))
            }else{
              tmp = .Score.test.meta(Y, X, X.index, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = Method)
              pval = cbind(pval, c(tmp$score.pvalue, tmp$score.Rpvalue) )
            }
          }

        }

      }
    }# lineage loop


  }# level loop


  colnames(pval) = subtree
  if(is.null(n.resample)){
    rownames(pval) = paste0("Asymptotic-",Method)
    score.tmp = pval[1,]
  }else{
    if(Method %in% c("RE-SKAT", "Het-SKAT")){
      rownames(pval) = paste0("Resampling-",Method)
      score.tmp = pval[1,]
    }else{
      rownames(pval) = c(paste0("Asymptotic-",Method),paste0("Resampling-",Method))
      score.tmp = pval[2,]
    }
    #print(pval)
  }

  # identify significant lineages
  subtree.tmp = subtree
  index.na = which(is.na(score.tmp))
  if(length(index.na)>0){
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


  return( list(lineage.pval=pval, sig.lineage=sig.lineage) )


}



