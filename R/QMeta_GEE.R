# Zero part model

.F.test <- function(x){
  # Fisher's p-value combination
  x.stat = -2 * sum(log(x))
  return( 1 - pchisq(x.stat, df = 2 * length(x)) )
}

.diag2 <- function(x){
 # transform the numeric into diag matrix
  if(length(x)>1){
    return(diag(x))
  }else{
    return(as.matrix(x))

  }

}

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

  U <- Score.reduce.reorg[ ,1:n.par.interest.alpha] - t(crossprod(t(Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha + 1):n.alpha)]),
                                                                 crossprod(t(ginv(Hess.reduce.reorg[((n.par.interest.alpha + 1):n.alpha), ((n.par.interest.alpha + 1):n.alpha)])), t(Score.reduce.reorg[ ,((n.par.interest.alpha + 1):n.alpha)]))))

  B2 <- matrix(0, n.alpha, n.alpha)

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
  score.pvalue.alpha = NULL
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
    n.par.alpha = length(par.index.alpha)
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
    Score.reduce.reorg = cbind( matrix(Score.reduce.alpha.perm[,par.index.alpha], ncol=n.par.alpha), matrix(Score.reduce.alpha.perm[,-par.index.alpha], ncol=n.alpha - n.par.alpha) )
    Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.alpha.perm[par.index.alpha, par.index.alpha], nrow=n.par.alpha), matrix(Hess.reduce.alpha.perm[par.index.alpha, -par.index.alpha], nrow=n.par.alpha) ),
                              cbind( matrix(Hess.reduce.alpha.perm[-par.index.alpha, par.index.alpha], nrow=n.alpha - n.par.alpha), matrix(Hess.reduce.alpha.perm[-par.index.alpha, -par.index.alpha], nrow= n.alpha - n.par.alpha)))


    A = colSums(Score.reduce.reorg)[1:n.par.interest.alpha]

    B1 <- ginv(Hess.reduce.reorg[(1:n.par.interest.alpha), (1:n.par.interest.alpha)] - crossprod(t(Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha + 1):n.alpha)]),
                                                                                                 crossprod(t(ginv(Hess.reduce.reorg[((n.par.interest.alpha + 1):n.alpha), ((n.par.interest.alpha + 1):n.alpha)])), Hess.reduce.reorg[((n.par.interest.alpha + 1):n.alpha), (1:n.par.interest.alpha)])))

    alpha.hat <- crossprod(t(B1), A)

    U <- Score.reduce.reorg[ ,1:n.par.interest.alpha] - t(crossprod(t(Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha + 1):n.alpha)]),
                                                                    crossprod(t(ginv(Hess.reduce.reorg[((n.par.interest.alpha + 1):n.alpha), ((n.par.interest.alpha + 1):n.alpha)])), t(Score.reduce.reorg[ ,((n.par.interest.alpha + 1):n.alpha)]))))

    B2 <- matrix(0, n.alpha, n.alpha)

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
  if(Method == "FE-MV"){
    score.stat.alpha.perm = crossprod( score.alpha.meta,crossprod(t(est.cov.inv), score.alpha.meta))
  }
  if(Method == "FE-VC"){
    score.stat.alpha.perm = crossprod(score.alpha.meta) #SKAT-VC
  }
  if(Method == "RE-MV"){
    U.theta = 0
    V.theta = 0
    for( i in 1:stu.num ){
      est.inv = ginv(est.cov.zero[[i]])
      U.theta = U.theta + 1/2 * crossprod(crossprod(t(est.inv), score.alpha[[i]])) - 1/2 * tr(est.inv)
      V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
    }
    score.stat.alpha.perm = crossprod( score.alpha.meta,crossprod(t(est.cov.inv), score.alpha.meta)) + U.theta^2/V.theta
  }
  if(Method == "RE-VC"){
    U.theta = 0
    V.theta = 0
    U.tau = 0
    est.inv.sum = 0
    for( i in 1:stu.num){
      est.inv = ginv(est.cov.zero[[i]])# calculate the generalized inverse for each estimate covariance for each study
      U.theta = U.theta +  1/2 * crossprod(crossprod(t(est.inv), score.alpha[[i]])) - 1/2 * tr(est.inv)
      # when \tau matrix and W matrix is identity the elements in upper left, bottom right as well as bottom left are the same in RE-VC test
      V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
      U.tau = U.tau - 1/2 * tr(est.inv)
      est.inv.sum = est.inv.sum + est.inv
    }
    V.tau = 1/2 * tr(crossprod(est.inv.sum))
    U.tau = 1/2 * crossprod(score.alpha.meta) + U.tau
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
  if(Method == "FE-MV"){
    score.stat.alpha.perm = crossprod( score.alpha.meta,crossprod(t(est.cov.inv), score.alpha.meta))
  }
  if(Method == "FE-VC"){
    score.stat.alpha.perm = crossprod(score.alpha.meta)
  }
  if(Method == "RE-MV"){
    U.theta = 0
    V.theta = 0
    for( i in 1:stu.num ){
      est.inv = ginv(est.cov.zero[[i]])
      U.theta = U.theta + 1/2 * crossprod(crossprod(t(est.inv), score.alpha[[i]])) - 1/2 * tr(est.inv)
      V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
    }
    score.stat.alpha.perm = crossprod( score.alpha.meta,crossprod(t(est.cov.inv), score.alpha.meta)) + U.theta^2/V.theta
  }
  if(Method == "RE-VC"){
    U.theta = 0
    V.theta = 0
    U.tau = 0
    est.inv.sum = 0
    for( i in 1:stu.num){
      est.inv = ginv(est.cov.zero[[i]])# calculate the generalized inverse for each estimate covariance for each study
      U.theta = U.theta +  1/2 * crossprod(crossprod(t(est.inv), score.alpha[[i]])) - 1/2 * tr(est.inv)
      # when \tau matrix and W matrix is identity the elements in upper left, bottom right as well as bottom left are the same in RE-VC test
      V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
      U.tau = U.tau - 1/2 * tr(est.inv)
      est.inv.sum = est.inv.sum + est.inv
    }
    V.tau = 1/2 * tr(crossprod(est.inv.sum))
    U.tau = 1/2 * crossprod(score.alpha.meta) + U.tau
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
      tmp = try( .Score.test.stat.zero.meta.4Gresampling.c(Z.perm.list, Z.par.index,n.par.interest.alpha, col.zero.index.list, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, Method = Method) )
    }
    else{
      tmp = try( .Score.test.stat.zero.meta.4Gresampling(Z.perm.list, Z.par.index,n.par.interest.alpha, col.zero.index.list, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, Method = Method) )
    }

    if(!("try-error" %in% class(tmp))){# if score.stat.alpha.perm exists, then n.one.new + 1
      # if the permutation test statistics greater than original test statistics, cnt + 1
      n.zero.new = n.zero.new + 1
      if(tmp >= score.stat.zero.meta){
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

.Score.test.zero.meta <- function(Y.list, Z.list, Z.par.index, seed=11, resample=FALSE, n.replicates=NULL, use.cpp = F, Method = "FE-MV"){
  stu.num = length(Z.list)
  p.par = length(Z.par.index)
  m = ncol(Y.list[[1]])
  n = nrow(Y.list[[1]])
  p.zero = ncol(Z.list[[1]])
  n.OTU = length(Y.list)
  # the total number of parameter of interest
  n.par.interest.alpha = m*length(Z.par.index)
  if(! Method %in% c("FE-MV","FE-VC",'RE-MV',"RE-VC")){
    stop("Error: Please Choose a Proper Meta-analysis Method")
  }
  # check the parameters of interest to ensure the intercept term is not in it
  if (length(unique(sapply(1:n.OTU,function(j) ncol(Y.list[[j]]))))!=1){
    stop("Error: The taxon in each study should be the same")
  }
  if(m<=1){
    stop("Error: Improper dimension for OTU table")
  }
  n.zero = 0
  # col.pos.index.lst = lapply(1:n.OTU, function(j) )
  ############################# Asymptotic: zero part
  Y0.list = Y.list
  for(j in 1:n.OTU)
  {
    # for each study set those value which are zero to 1
    Y0.list[[j]][Y.list[[j]]==0] = 1
    # for each study set those value which are non-zero to 0
    Y0.list[[j]][Y.list[[j]]>0] = 0
  }

  # if all 0 in one group across across taxa, then output NA
  #if( (ncol(Y0)-length(remove.index))<=1 | sum(Y0[case==1,])==0 | sum(Y0[case==0,])==0){
  if(is.null(Z.par.index) || n==0){

    score.pvalue.zero = NA
    score.stat.zero.meta = NA
    df.zero = NA

  }else{
    # initialize for later use
    score.pvalue.zero = NA
    score.stat.zero.meta = NA
    df.zero = NA
    score.stat.alpha = NULL
    score.pvalue.alpha = NULL
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
        n.zero = n.zero + 1
        # get the index for parameter of interest after score statistics and estimate covariance matrix are reorginzed
        # different across studies because of different column numbers
        idx = kronecker((col.zero.index-1)*p.par, rep(1,p.par)) + c(1:p.par)
        col.zero.index.list[[n.zero]] = idx  # the index of alpha parameter of interest in this study
        # save these outcome in list form for later meta-analysis as well as resampling test
        score.stat.alpha = append(score.stat.alpha, tmp.zero$score.stat.alpha)
        score.alpha[[n.zero]] = tmp.zero$score.alpha
        est.cov.zero[[n.zero]] = tmp.zero$est.cov.zero
        score.pvalue.alpha = append(score.pvalue.alpha, (1 - pchisq(tmp.zero$score.stat.alpha, ncol(Y0)*length(Z.par.index))))
        zero.vA.list.meta[[n.zero]] = tmp.zero$vA.list
        zero.Vinv.list.meta[[n.zero]] = tmp.zero$Vinv.list
        zero.VY.list.meta[[n.zero]] = tmp.zero$VY.list
        score.alpha.meta[idx] =  score.alpha.meta[idx] + crossprod(ginv(cov.alpha), alpha.hat) # add according to the index of parameter of interest for each study
        est.cov.meta[idx, idx] =  est.cov.meta[idx, idx] + ginv(cov.alpha)
      }
    }
  }
  ############################# Asymptotic: combined
  if(n.zero>0){
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
    if(Method == "FE-MV"){
      score.stat.meta = crossprod( score.alpha.meta,crossprod(t(est.cov.inv), score.alpha.meta))
      score.pvalue = 1- pchisq(score.stat.meta,df = n.par.save.alpha )
      df = n.par.save.alpha
    }
    if(Method == "FE-VC"){
      weight.cov.meta = eigen(est.cov.meta)$values #eign.val/sum(eign.val)
      score.stat.meta = crossprod(score.alpha.meta) #SKAT-VC
      score.pvalue = davies(score.stat.meta, weight.cov.meta, h = rep(1,n.par.save.alpha), delta = rep(0,n.par.save.alpha), sigma = 0, lim = 10000, acc = 0.0001)$Qq
      score.pvalue = ifelse(score.pvalue>0,score.pvalue,1e-7)
      df = n.par.save.alpha
    }
    if(Method == "RE-MV"){
      U.theta = 0
      V.theta = 0
      for( i in 1:n.zero ){
        est.inv = ginv(est.cov.zero[[i]])
        U.theta = U.theta + 1/2 * crossprod(crossprod(t(est.inv), score.beta[[i]])) - 1/2 * tr(est.inv)
        V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
      }
      score.stat.meta = crossprod( score.beta.meta,crossprod(t(est.cov.inv), score.beta.meta)) + U.theta^2/V.theta
    }
    if(Method == "RE-VC"){
      U.theta = 0
      V.theta = 0
      U.tau = 0
      est.inv.sum = 0
      for( i in 1:n.zero ){
        est.inv = ginv(est.cov.zero[[i]])# calculate the generalized inverse for each estimate covariance for each study
        U.theta = U.theta +  1/2 * crossprod(crossprod(t(est.inv), score.beta[[i]])) - 1/2 * tr(est.inv)
        # when \tau matrix and W matrix is identity the elements in upper left, bottom right as well as bottom left are the same in RE-VC test
        V.theta = V.theta + 1/2 * tr(crossprod(est.inv))
        U.tau = U.tau - 1/2 * tr(est.inv)
        est.inv.sum = est.inv.sum + est.inv
      }
      V.tau = 1/2 * tr(crossprod(est.inv.sum))
      U.tau = 1/2 * crossprod(score.beta.meta) + U.tau
      score.stat.meta = crossprod(c(U.tau, U.theta), crossprod(ginv(matrix(c(V.tau,rep(V.theta,3)),ncol = 2)), c(U.tau, U.theta)))
    }
    zero.results = list(score.stat = score.stat.zero.meta, score.pvalue = score.pvalue.zero, df = df.zero)
  }
  # if resample = TRUE then will apply permutation method to get permuted p-value
  # adaptive resampling test
  if(resample){

    #print("simulated stat:")
    set.seed(seed)
    if(!is.na(score.stat.zero.meta)){

      n.zero = 0
      zero.acc = 0

      start.nperm = 1;
      end.nperm = min(100,n.replicates);
      flag = 1
      while(flag & end.nperm <= n.replicates){


        results = .resample.work.zero.meta(Z.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, score.stat.zero.meta, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, start.nperm, end.nperm, n.zero, zero.acc, use.cpp = use.cpp, Method = Method)
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

          results = .resample.work.zero.meta(Z.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, score.stat.zero.meta, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, start.nperm, end.nperm, n.zero, zero.acc, use.cpp = use.cpp, Method = Method)

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
#
# .Rarefy <- function (otu.tab, depth = min(rowSums(otu.tab))){
#   # Rarefaction function: downsample to equal depth
#   #
#   # Args:
#   #		otu.tab: OTU count table, row - n sample, column - q OTU
#   #		depth: required sequencing depth
#   #
#   # Returns:
#   # 	otu.tab.rff: Rarefied OTU table
#   #		discard: labels of discarded samples
#   #
#   otu.tab <- as.matrix(otu.tab)
#   ind <- (rowSums(otu.tab) < depth)
#   sam.discard <- rownames(otu.tab)[ind]
#   otu.tab <- otu.tab[!ind, ]
#
#   rarefy <- function(x, depth){
#     y <- sample(rep(1:length(x), x), depth)
#     y.tab <- table(y)
#     z <- numeric(length(x))
#     z[as.numeric(names(y.tab))] <- y.tab
#     z
#   }
#   otu.tab.rff <- t(apply(otu.tab, 1, rarefy, depth))
#   rownames(otu.tab.rff) <- rownames(otu.tab)
#   colnames(otu.tab.rff) <- colnames(otu.tab)
#   return(list(otu.tab.rff=otu.tab.rff, discard=sam.discard))
# }

#' Title
#' @inheritParams QCAT_Meta
#' @param Z a list of matrices contains covariates for the zero-part test with each column pertains to one variable (pertains to the covariate of interest or the confounders). The number of elements of Z and OTU must be the same. The column number of each matrix in this list must be the same.
#' @param Z.index a vector indicate the columns in X for the covariate(s) of interest.
#'
#' @return A list with this elements
#'    \item{lineage.pval}{p-values for all lineages. By default ( Method = "FE-MV", n.perm = NULL ), only the asymptotic test will be performed. If the meta-analysis method is random effect methods like RE-VC or RE-MV, then resampling test must be performed.}
#'    \item{sig.lineage}{a vector of significant lineages.}
#' @export
#'
#' @examples
#' data(data.meta)
#' OTU = data.meta$OTU
#' Tax = data.meta$Tax
#' case = data.meta$covariate
#' QCAT_GEE_Meta(OTU, case, 1, Tax, Method = "FE-MV", min.depth=0, n.perm=NULL, fdr.alpha=0.05)
#' @import MASS
#' @import data.table
#' @import CompQuadForm
#' @import geepack
#' @importFrom stats coef optim pchisq
#' @importFrom dplyr bind_rows
#' @references
#' Lee S, Teslovich TM, Boehnke M, Lin X. (2013) General framework for meta-analysis of rare variants in sequencing association studies. Am J Hum Genet.
#' \emph{Am J Hum Genet}
#' \doi{10.1016/j.ajhg.2013.05.010}.
#' @references
#' Benjamini, Yoav, and Yosef Hochberg.(1995) Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.
#' \emph{Journal of the Royal Statistical Society. Series B}
#' @references
#' Zeger, Scott L., and Kung-Yee Liang. (1986) Longitudinal Data Analysis for Discrete and Continuous Outcomes.
#' \emph{Biometrics}
#' \doi{https://doi.org/10.2307/2531248}.
QCAT_GEE_Meta <- function(OTU, Z, Z.index, Tax=NULL, Method = "FE-MV", min.depth=0, n.perm=NULL, use.cpp = F, fdr.alpha=0.05){
  n.OTU = length(OTU)
  n.resample = n.perm
  # if(length(unique(sapply(1:n.OTU,function(j) ncol(OTU[[j]]))))!= 1){
  #   stop("The taxa in each study should be the same")
  # }
  if(Method %in% c("RE-VC", "RE-MV")){
    if(is.null(n.resample)){
      stop("The p-value for random effect meta-analysis method must be got by resampling test")
    }
  }
  n.Z = length(Z)
  if(n.OTU!=n.Z)
  {
    stop("The study number of OTU table and Covariate table for the
       zero part should be the same")
  }

  for(i in 1:n.OTU)
  {
    if(!is.matrix(OTU[[i]])){
      OTU[[i]] = as.matrix(OTU[[i]])
    }

    if(!is.matrix(Z[[i]])){
    Z[[i]] = as.matrix(Z[[i]])
    }

    if(nrow(OTU[[i]])!=nrow(Z[[i]])){
      stop(paste0("Number of samples in the OTU table and the covariate table for the zero part
                  of study ", i, " should be the same"))
    }

    remove.subject = which(rowSums(OTU[[i]])<=min.depth)
    if(length(remove.subject)>0){
      print(paste("Remove",length(remove.subject), "samples with read depth less or equal to", min.depth, "in OTU table", i))
      Z[[i]] = Z[[i]][-remove.subject, ,drop=FALSE]
      OTU[[i]] = OTU[[i]][-remove.subject, ,drop=FALSE]
    }
    # drop.col = union(drop.col,which(colSums(OTU[[i]])==0))
  }

  if(missing(Z.index)){
    Z.index = 1:ncol(Z[[1]])
  }

  # if(length(drop.col)>0){
  #   count = lapply(1:n.OTU, function(j) OTU[[j]][,-drop.col, drop=FALSE])
  # }else{
  #   count = OTU
  # }
  #
  OTU.comb = as.data.frame(OTU[[1]])
  for (i in 2:n.OTU){
    OTU.comb = bind_rows(OTU.comb,as.data.frame(OTU[[i]]))
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
    count[[i]] = OTU.comb[batch.cnt[i]:(batch.cnt[i+1]-1),]
  }

  Z = lapply(1:n.OTU,function(j) cbind(1, Z[[j]])) # add the intercept term
  Z.index = Z.index + 1

  if(is.null(Tax)){ # perform one test using all OTUs
    # count = lapply(count,function(X){.Rarefy(X)$otu.tab.rff})
    if(is.null(n.resample)){ # asymptotic test only
      tmp = .Score.test.zero.meta(count, Z, Z.index, seed=11, resample=FALSE, n.replicates=NULL, Method = Method)
      pval.zero = as.matrix( tmp$score.pvalue )
      colnames(pval.zero) = paste0("Asymptotic-",Method)
    }else{
      tmp = .Score.test.zero.meta(count, Z, Z.index, seed=11, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = Method)
      if(Method %in% c("RE-VC", "RE-MV")){
        pval.zero = c(tmp$score.Rpvalue)
        names(pval.zero) = paste0("Resampling-",Method)
      }else{
        pval.zero = c(tmp$score.pvalue, tmp$score.Rpvalue)
        names(pval.zero) = c(paste0("Asymptotic-",Method),paste0("Resampling-",Method))
      }
    }
    return(list(pval = pval.zero))
  }else{ # perform tests for lineages

    if(!is.matrix(Tax)){
      tax = as.matrix(Tax)
    }

    for(i in 1:n.OTU)
    {
      if( sum(!(colnames(count[[i]]) %in% rownames(tax)))>0 ){
        stop(paste0("Error: OTU IDs in OTU table ",i," are not consistent with OTU IDs in Tax table"))
      }
    }

    n.rank = ncol(tax)
    #merge the tax and count together for later partition and merge OTU table according to taxonomy information
    W.data.list = lapply(1:n.OTU,function(j) data.table(data.frame(tax, t(count[[j]]))))
    otucols = lapply(1:n.OTU,function(j) names(W.data.list[[j]])[-(1:n.rank)])
    n.level = n.rank-1

    subtree = NULL
    pval.zero = NULL

    for(k in 1:n.level){

      #print(k)
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


      for(j in 1:m.level){

        Y = lapply(1:n.OTU, function(i) t(W.count[[i]][which(W.tax == level.uni[j]), , drop=FALSE]))
        # remove.list.idx = NULL
        # for(i in 1:n.OTU)
        # {
        #   remove.subject = which(rowSums(Y[[i]])<=0)
        #   if(length(remove.subject)>0){
        #     if(length(remove.subject) == nrow(Y[[i]])){
        #       remove.list.idx = append(remove.list.idx,i)
        #     }else{
        #       Z[[i]] = Z[[i]][-remove.subject, ,drop=FALSE]
        #       Y[[i]] = Y[[i]][-remove.subject, ,drop=FALSE]
        #     }
        #   }
        # }
        # if(length(remove.list.idx)==length(Y)){
        #   next
        # }
        # if(length(remove.list.idx)>0){
        #   Y = Y[-remove.list.idx]
        #   Z = Z[-remove.list.idx]
        # }
        # Y = lapply(Y,function(X){.Rarefy(X)$otu.tab.rff})
        #Y = t(W.count[which(W.tax == "f__Veillonellaceae"), , drop=FALSE])
        # remove.index = NULL
        # for( i in 1:n.OTU)
        # {
        #   remove.index = union(remove.index,which(colSums(Y[[i]])==0))
        # }
        #
        # if(length(remove.index)==ncol(Y[[1]])){
        #   next
        # }else{
        #
        #   if(length(remove.index)>0){
        #     Y = lapply(1:n.OTU, function(i) Y[[i]][, -remove.index, drop=FALSE])
        #   }
        # }

        if(ncol(Y[[1]])==1){
          next
          #print("==skip:1==");
        }else{

          subtree = c(subtree, level.uni[j])
          if(is.null(n.resample)){
            #  run test for each lineage
            tmp = .Score.test.zero.meta(Y, Z, Z.index, seed=11, resample=FALSE, n.replicates=n.resample, Method = Method)
            pval.zero = cbind(pval.zero, tmp$score.pvalue)
          }
          else{
            # if n.resample in not null, select the significant lineage according to the resampling pvalue
            #  run test for each lineage
            if(Method %in% c("RE-VC", "RE-MV")){
              tmp = .Score.test.zero.meta(Y, Z, Z.index, seed=11, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = Method)
              pval.zero = cbind(pval.zero, tmp$score.Rpvalue)
            }else{
              tmp = .Score.test.zero.meta(Y, Z, Z.index, seed=11, resample=TRUE, n.replicates=n.resample, use.cpp = use.cpp, Method = Method)
              pval.zero = cbind(pval.zero, c(tmp$score.pvalue, tmp$score.Rpvalue) )
            }
          }

        }

      }# lineage loop
    } # level loop
    colnames(pval.zero) = subtree

    if(is.null(n.resample)){

      rownames(pval.zero) =  paste0("Asymptotic-",Method)
      score.zero.tmp = pval.zero[1,]

    }else{
      if(Method %in% c("RE-VC", "RE-MV")){
        rownames(pval.zero) = paste0("Resampling-",Method)
        score.zero.tmp = pval.zero[1,]
      }else{
        rownames(pval.zero) = c(paste0("Asymptotic-",Method),paste0("Resampling-",Method))
        score.zero.tmp = pval.zero[2,]
      }
    }
    # identify significant lineages

    subtree.tmp = subtree
    index.na = which(is.na(score.zero.tmp))
    if(length(index.na)>0){
      # drop those lineages which have NA values
      score.zero.tmp = score.zero.tmp[-index.na]
      subtree.tmp = subtree.tmp[-index.na]
    }

    #score.tmp[score.tmp==0] = 1e-4
    m.test = length(score.zero.tmp)

    # Benjamini-Hochberg FDR control
    index.p = order(score.zero.tmp)
    p.sort = sort(score.zero.tmp)
    #fdr.alpha = 0.05

    # change 12/10/2022
    reject = rep(0, m.test)
    tmp = which(p.sort<=(1:m.test)*fdr.alpha/m.test)
    if(length(tmp)>0){
      index.reject = index.p[1:max(tmp)]
      reject[index.reject] = 1
    }

    sig.lineage = subtree.tmp[reject==1]
    # return all the p-values as well as significant lineages
    return( list(lineage.pval=pval.zero, sig.lineage=sig.lineage) )
  }

}


