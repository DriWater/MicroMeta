# Zero part model

.F.test <- function(x){

  x.stat = -2 * sum(log(x))
  return( 1 - pchisq(x.stat, df = 2 * length(x)) )
}

.simes.test <- function(x){

  return( min(length(x) * x/rank(x)) )

}

.diag2 <- function(x){

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
      A.i = .diag2(vA.tmp)
      t.D.i = kronecker( A.i, as.matrix(Z[i,], ncol=1) )
      V.i = A.i # independent cor structure

      tmp.V.i = ginv(V.i)
      tmp.VY = tmp.V.i %*% (Y[i,] - Pi.i)
      Score.alpha.i[i,] = t.D.i %*% tmp.VY

      if(save.list){
        vA.list[[i]] = vA.tmp
        Vinv.list[[i]] = tmp.V.i
        VY.list[[i]] = tmp.VY
      }
    }


  }
  # save list for later permutation method
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


      Hessian.alpha = Hessian.alpha + crossprod(t(t.D.i), tcrossprod(ginv(V.i), t.D.i))


    }

    return (Hessian.alpha)


  }


}

.Score.test.stat.zero <- function(Y0, Z, Z.par.index, cor.stru){


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
  gee.reduce = geeglm(outcome ~ . - id - 1, data = data.reduce, id = factor(id), family="binomial", corstr= "independence")
  #wald.test = anova(gee.full, gee.reduce)


  ########### perform score test
  n.alpha = m *p
  par.interest.index.alpha =  kronecker( ((0:(m-1))*p), rep(1,length(Z.par.index))) + Z.par.index
  n.par.interest.alpha = length(par.interest.index.alpha)
  est.reduce.alpha = rep(NA, n.alpha)
  est.reduce.alpha[par.interest.index.alpha] = 0
  est.reduce.alpha[-par.interest.index.alpha] = coef(gee.reduce)
  est.reduce.scale = gee.reduce

  data.alpha = list(Y=Y0, Z=Z)

  tmp = .fun.score.i.alpha(est.reduce.alpha, data.alpha, save.list=TRUE)
  Score.reduce.alpha = tmp$Score.alpha
  # for resampling test
  vA.list = tmp$vA.list
  Vinv.list = tmp$Vinv.list
  VY.list = tmp$VY.list

  Hess.reduce.alpha =  .fun.hessian.alpha(est.reduce.alpha, data.alpha)
  # re-organized the score statistics and Hessian matrix
  Score.reduce.reorg = cbind( matrix(Score.reduce.alpha[,par.interest.index.alpha], ncol=n.par.interest.alpha), matrix(Score.reduce.alpha[,-par.interest.index.alpha], ncol=n.alpha - n.par.interest.alpha) )
  Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.alpha[par.interest.index.alpha, par.interest.index.alpha], nrow=n.par.interest.alpha), matrix(Hess.reduce.alpha[par.interest.index.alpha, -par.interest.index.alpha], nrow=n.par.interest.alpha) ),
                            cbind( matrix(Hess.reduce.alpha[-par.interest.index.alpha, par.interest.index.alpha], nrow=n.alpha - n.par.interest.alpha), matrix(Hess.reduce.alpha[-par.interest.index.alpha, -par.interest.index.alpha], nrow= n.alpha - n.par.interest.alpha)))


  A = colSums(Score.reduce.reorg)[1:n.par.interest.alpha]

  B1 = cbind(diag(n.par.interest.alpha), tcrossprod(-Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha+1):n.alpha)] , t(ginv(Hess.reduce.reorg[((n.par.interest.alpha+1):n.alpha), ((n.par.interest.alpha+1):n.alpha)]))))

  B2 =  matrix(0, n.alpha, n.alpha)
  for(i in 1:n){
    B2 = B2 + Score.reduce.reorg[i,] %o% Score.reduce.reorg[i,]
  }

  B = crossprod(t(B1), tcrossprod(B2, B1))
  score.stat.alpha = crossprod( A, crossprod(t(ginv(B)), A))
  score.pvalue.alpha = 1 - pchisq(score.stat.alpha, n.par.interest.alpha)


  return(list(score.df.alpha=n.par.interest.alpha, score.stat.alpha = score.stat.alpha, score.alpha = A, est.cov.zero=B, score.pvalue.alpha=score.pvalue.alpha, vA.list=vA.list, Vinv.list=Vinv.list, VY.list=VY.list )   )

}

.Score.test.stat.zero.meta.4Gresampling <- function(Z.perm.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, vA.list.meta, Vinv.list.meta, VY.list.meta, Method = "FE-MV", W.zero = NULL){

  # W = W.zero
  iter.num = length(Z.perm.list)
  score.stat.alpha = NULL
  score.alpha = NULL
  score.pvalue.alpha = NULL
  est.cov.zero = NULL
  score.alpha.meta  = rep(0,n.par.interest.alpha) ## A
  est.cov.meta = matrix(0, nrow = n.par.interest.alpha, ncol = n.par.interest.alpha) ## B
  for ( j in 1:iter.num){
    vA.list = vA.list.meta[[j]]
    VY.list = VY.list.meta[[j]]
    Vinv.list = Vinv.list.meta[[j]]
    Z.perm = Z.perm.list[[j]]
    p = ncol(Z.perm.list[[j]])
    m.alpha = length(vA.list[[1]])
    n.alpha = m.alpha*p
    par.index.alpha =  kronecker( ((0:(m.alpha-1))*p), rep(1,length(Z.par.index))) + Z.par.index
    n.par.alpha = length(par.index.alpha)
    n = nrow(Z.perm)
    Score.reduce.alpha.perm = matrix(0, n, n.alpha )
    Hess.reduce.alpha.perm = matrix(0, n.alpha, n.alpha )

    for(i in 1:n){

      ###################################################
      #                                                 #
      #         alpha part: resampling Score test        #
      #                                                 #
      ###################################################
      tD.tmp = kronecker(.diag2(vA.list[[i]]), as.matrix(Z.perm[i,], ncol=1))

      Score.reduce.alpha.perm[i,] = Score.reduce.alpha.perm[i,] + crossprod(t(tD.tmp),VY.list[[i]])

      Hess.reduce.alpha.perm = Hess.reduce.alpha.perm + crossprod(t(tD.tmp), tcrossprod(Vinv.list[[i]], tD.tmp))


    }

    # re-organized the score statistics and Hessian matrix
    Score.reduce.reorg = cbind( matrix(Score.reduce.alpha.perm[,par.index.alpha], ncol=n.par.alpha), matrix(Score.reduce.alpha.perm[,-par.index.alpha], ncol=n.alpha - n.par.alpha) )
    Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.alpha.perm[par.index.alpha, par.index.alpha], nrow=n.par.alpha), matrix(Hess.reduce.alpha.perm[par.index.alpha, -par.index.alpha], nrow=n.par.alpha) ),
                              cbind( matrix(Hess.reduce.alpha.perm[-par.index.alpha, par.index.alpha], nrow=n.alpha - n.par.alpha), matrix(Hess.reduce.alpha.perm[-par.index.alpha, -par.index.alpha], nrow= n.alpha - n.par.alpha)))


    A = colSums(Score.reduce.reorg)[1:n.par.alpha]

    B1 = cbind(diag(n.par.alpha), tcrossprod(-Hess.reduce.reorg[(1:n.par.alpha), ((n.par.alpha+1):n.alpha)], t(ginv(Hess.reduce.reorg[((n.par.alpha+1):n.alpha), ((n.par.alpha+1):n.alpha)]))))

    B2 =  matrix(0, n.alpha, n.alpha)
    for(i in 1:n){
      B2 = B2 + Score.reduce.reorg[i,] %o% Score.reduce.reorg[i,]
    }

    B = crossprod(t(B1), tcrossprod(B2, B1))
    score.stat.alpha.perm = crossprod( A, crossprod(t(ginv(B)), A))
    score.stat.alpha = append(score.stat.alpha, score.stat.alpha.perm)
    score.alpha[[j]] = A
    est.cov.zero[[j]] = B
    score.pvalue.alpha = append(score.pvalue.alpha, (1 - pchisq(score.stat.alpha.perm, n.par.interest.alpha)))
    idx = col.zero.index.list[[j]]
    score.alpha.meta[idx] =  score.alpha.meta[idx] +  A
    est.cov.meta[idx, idx] =  est.cov.meta[idx, idx] + B
  }

  save.index.zero = which(abs(score.alpha.meta) >= 1e-7)
  n.par.save.alpha = length(save.index.zero)
  score.alpha.meta =  score.alpha.meta[save.index.zero]
  est.cov.meta =  est.cov.meta[save.index.zero,save.index.zero]
  est.cov.inv = ginv(est.cov.meta)
  weight.cov.zero = NULL
  if (Method == "FE-MV")
  {
    score.stat.alpha.perm = crossprod( score.alpha.meta, crossprod(t(est.cov.inv), score.alpha.meta)) #FE-Burden
  }
  # if (Method == "SKAT")
  # {
  #   if(is.null(W))
  #   {
  #     W = diag(1,nrow = n.par.save.alpha)
  #   }
  #   else{
  #     W = W[save.index.zero,save.index.zero]
  #   }
  #   eigen.cov.zero <- eigen(est.cov.meta)
  #   eigen.cov.zero.sqrt = eigen.cov.zero$vectors %*% diag(sqrt(eigen.cov.zero$values),nrow = length(eigen.cov.zero$values)) %*% solve(eigen.cov.zero$vectors)
  #   weight.cov.zero = eigen(eigen.cov.zero.sqrt %*% W %*% eigen.cov.zero.sqrt)$values
  #   score.stat.alpha.perm = score.alpha.meta %*% W %*% score.alpha.meta
  # }
  if (Method == "FE-VC")
  {
    # weight.cov.zero = eigen(est.cov.inv)$values #eign.val/sum(eign.val)
    #
    score.stat.alpha.perm = crossprod( score.alpha.meta, crossprod(t(est.cov.inv), crossprod( t(est.cov.inv), score.alpha.meta))) #SKAT-VC
  }
  if(Method == "Het-SKAT"){
    #W = diag(1,nrow = n.par.interest.beta)
    score.stat.alpha.perm = 0

    for( i in 1:iter.num ){
      est.inv = ginv(est.cov.zero[[i]])
      score.stat.alpha.perm = score.stat.alpha.perm + tcrossprod(crossprod(score.alpha[[i]], est.inv), crossprod(score.alpha[[i]], est.inv))
    }
  }
  if(Method == "RE-SKAT"){
    U.tau.b = 0
    a2 = 0
    a3 = 0
    a4 = 0
    for( i in 1:iter.num ){
      est.inv = ginv(est.cov.zero[[i]])
      U.tau.b = U.tau.b + tcrossprod(crossprod(score.alpha[[i]],est.inv), crossprod(score.alpha[[i]], est.inv))
      a2 = a2 + sum(diag(crossprod(t(est.inv),est.inv)))
      a3 = a3 + sum(diag(crossprod(t(est.inv),est.inv)))
      a4  = a4 + sum(diag(crossprod(t(est.inv),est.inv)))
    }
    a1 = sum(diag(crossprod(t(est.cov.inv),est.cov.inv)))
    U.tau.b = 1/2 * U.tau.b + 1/2 * sum(diag(est.cov.inv))
    U.tau.w = 1/2 * crossprod( score.alpha.meta, crossprod(t(est.cov.inv), crossprod( t(est.cov.inv), score.alpha.meta)))
    + 1/2 * sum(diag(est.cov.inv))  #SKAT-VC
    score.stat.alpha.perm = crossprod(c(U.tau.w, U.tau.b), crossprod(1/2 * matrix(c(a1,a2,a3,a4),ncol = 2), c(U.tau.w, U.tau.b)))
  }
  # if (Method == "Fisher")
  # {
  #   score.stat.alpha.perm = -2 * sum(log(score.pvalue.alpha))
  # }

  return(list(score.stat.alpha.perm = as.numeric(score.stat.alpha.perm),weight.cov.zero = weight.cov.zero))

}

#Rcpp::List score_test_stat_meta_resampling_c(const Rcpp::List& Z_perm_list, const Rcpp::List& col_zero_index_list, const Rcpp::List& vA_list_meta, const Rcpp::List& Vinv_list_meta,
#                                             const Rcpp::List& VY_list_meta, const arma::vec& Z_par_index,int n_par_interest_alpha)

.Score.test.stat.zero.meta.4Gresampling.c <- function(Z.perm.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, vA.list.meta, Vinv.list.meta, VY.list.meta, Method = "FE-MV", W.zero = NULL){
  tmp = score_test_stat_zero_meta_resampling_c(Z.perm.list, col.zero.index.list, vA.list.meta, Vinv.list.meta, VY.list.meta, Z.par.index, n.par.interest.alpha)
  est.cov.meta = tmp$est_cov_meta
  score.alpha.meta = tmp$score_alpha_meta
  est.cov = tmp$est_cov
  score.alpha = tmp$score_alpha
  save.index.zero = which(abs(score.alpha.meta) >= 1e-7)
  score.alpha.meta = score.alpha.meta[save.index.zero]
  est.cov.meta = est.cov.meta[save.index.zero,save.index.zero]
  est.cov.inv = ginv(est.cov.meta)
  iter.num = length(Z.perm.list)
  if (Method == "FE-MV")
  {
    score.stat.alpha.perm = crossprod(score.alpha.meta, crossprod(t(est.cov.inv), score.alpha.meta)) #FE-Burden
  }
  # if (Method == "SKAT")
  # {
  #   if(is.null(W))
  #   {
  #     W = diag(1,nrow = n.par.save.alpha)
  #   }
  #   else{
  #     W = W[save.index.zero,save.index.zero]
  #   }
  #   eigen.cov.zero <- eigen(est.cov.meta)
  #   eigen.cov.zero.sqrt = eigen.cov.zero$vectors %*% diag(sqrt(eigen.cov.zero$values),nrow = length(eigen.cov.zero$values)) %*% solve(eigen.cov.zero$vectors)
  #   weight.cov.zero = eigen(eigen.cov.zero.sqrt %*% W %*% eigen.cov.zero.sqrt)$values
  #   score.stat.alpha.perm = score.alpha.meta %*% W %*% score.alpha.meta
  # }
  if (Method == "FE-VC")
  {
    #weight.cov.zero = eigen(est.cov.inv)$values #eign.val/sum(eign.val)
    score.stat.alpha.perm = crossprod(score.alpha.meta, crossprod(t(est.cov.inv ), crossprod(t(est.cov.inv ), score.alpha.meta))) #SKAT-VC
  }
  # if (Method == "Fisher")
  # {
  #   score.stat.alpha.perm = -2 * sum(log(score.pvalue.alpha))
  # }
  if(Method == "Het-SKAT"){
    #W = diag(1,nrow = n.par.interest.beta)
    score.stat.alpha.perm = 0

    for( i in 1:iter.num ){
      est.inv = ginv(est.cov[[i]])
      score.stat.alpha.perm = score.stat.alpha.perm + tcrossprod(crossprod(score.alpha[[i]], est.inv), crossprod(score.alpha[[i]], est.inv))
    }
  }
  if(Method == "RE-SKAT"){
    U.tau.b = 0
    a2 = 0
    a3 = 0
    a4 = 0
    for( i in 1:iter.num ){
      est.inv = ginv(est.cov[[i]])
      U.tau.b = U.tau.b + tcrossprod(crossprod(score.alpha[[i]],est.inv), crossprod(score.alpha[[i]], est.inv))
      a2 = a2 + sum(diag(crossprod(t(est.inv),est.inv)))
      a3 = a3 + sum(diag(crossprod(t(est.inv),est.inv)))
      a4  = a4 + sum(diag(crossprod(t(est.inv),est.inv)))
    }
    a1 = sum(diag(crossprod(t(est.cov.inv),est.cov.inv)))
    U.tau.b = 1/2 * U.tau.b + 1/2 * sum(diag(est.cov.inv))
    U.tau.w = 1/2 * crossprod( score.alpha.meta, crossprod(t(est.cov.inv), crossprod( t(est.cov.inv), score.alpha.meta)))
    + 1/2 * sum(diag(est.cov.inv))  #SKAT-VC
    score.stat.alpha.perm = crossprod(c(U.tau.w, U.tau.b), crossprod(1/2 * matrix(c(a1,a2,a3,a4),ncol = 2), c(U.tau.w, U.tau.b)))
  }
  return(list(score.stat.alpha.perm = as.numeric(score.stat.alpha.perm)))

}


# add 10/22/2022 for adaptive resampling
.resample.work.zero.meta <- function(Z.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, score.stat.zero.meta, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, start.nperm, end.nperm, n.zero, zero.acc, Method = "FE-MV", W.zero = NULL){

  iter.num = length(Z.list)
  n.zero.new = n.zero
  zero.acc.new = zero.acc

  for(k in start.nperm:end.nperm){

    perm.index = lapply(Z.list, function(i) sample(1:nrow(i)))

    Z.perm.list = Z.list
    for (i in 1:iter.num){
      Z.perm.list[[i]][,Z.par.index] =  Z.perm.list[[i]][perm.index[[i]],Z.par.index]
    }

    tmp = try( .Score.test.stat.zero.meta.4Gresampling.c(Z.perm.list, Z.par.index,n.par.interest.alpha, col.zero.index.list, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, Method = Method, W.zero = W.zero) )



    if(!("try-error" %in% class(tmp))){

      n.zero.new = n.zero.new + 1
      if(tmp$score.stat.alpha.perm >= score.stat.zero.meta){
        zero.acc.new = zero.acc.new + 1

      }
    }


  }

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

# add 06/01/2016: function Score.test2 to run two-part test in the general setting(zero part: use score GEE test)
# Y: nxm count of microbiomes
# X: covariates for positive part: first column is always intercept
# Z: covariates for zero part: first column is always intercept
# Z.par.index: index for the parameter of interest for the Z part

.Score.test.zero.meta <- function(Y.list, Z.list, Z.par.index, seed=11, resample=FALSE, n.replicates=NULL, Method = "FE-MV", Weight.zero=NULL){
  iter.num = length(Z.list)
  W.zero = Weight.zero
  p.par = length(Z.par.index)
  m = ncol(Y.list[[1]])
  n = nrow(Y.list[[1]])
  p.zero = ncol(Z.list[[1]])
  n.OTU = length(Y.list)
  n.par.interest.alpha = m*length(Z.par.index)
  if(! Method %in% c("FE-MV","FE-VC",'Het-SKAT',"RE-SKAT")){
    stop("Error: Please Choose a Proper Meta-analysis Method")
  }
  if (length(unique(sapply(1:n.OTU,function(j) ncol(Y.list[[j]]))))!=1){
    stop("Error: The taxon in each study should be the same")
  }
  if(!is.null(W.zero))
  {
    if(!is.matrix(W.zero))
    {
      W.zero = as.matrix(W.zero)
      warning("The Weight Matrix for positvie part should be a matrix")
    }
    if (dim(W.zero)[1] != n.par.interest.alpha | dim(W.zero)[2] != n.par.interest.alpha)
    {
      stop("Error: The dimesion of Weight Matrix of SKAT Method for positive part should
           equal to the number of beta parameter of interest")
    }
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
    Y0.list[[j]][Y.list[[j]]==0] = 1
    Y0.list[[j]][Y.list[[j]]>0] = 0
  }

  # if all 0 in one group across across taxa, then output NA
  #if( (ncol(Y0)-length(remove.index))<=1 | sum(Y0[case==1,])==0 | sum(Y0[case==0,])==0){
  if(is.null(Z.par.index) || n==0){

    score.pvalue.zero = NA
    score.stat.zero.meta = NA
    df.zero = NA

  }else{

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
    score.alpha.meta  = rep(0,n.par.interest.alpha) ## A
    est.cov.zero.meta = matrix(0, nrow = n.par.interest.alpha, ncol = n.par.interest.alpha) ## B
    for(i in 1:iter.num){
      Y0 = Y0.list[[i]]
      Z = Z.list[[i]]
      col.zero.index = which(apply(Y0, 2, function(x) length(table(x)) ) > 1)
      Y0 = Y0[, col.zero.index , drop=FALSE]
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

        next

      }else{
        n.zero = n.zero + 1
        # par.index.pos.list[[n.pos]] = kronecker(((0:(m-2))*p), rep(1,length(X1.par.index))) + X1.par.index
        #par.index.zero.list[[n.zero]] = kronecker((col.zero.index - 1)*p.zero, rep(1,length(Z.par.index))) + Z.par.index
        idx = kronecker((col.zero.index-1)*p.par, rep(1,p.par)) + c(1:p.par)
        col.zero.index.list[[n.zero]] = idx
        score.stat.alpha = append(score.stat.alpha, tmp.zero$score.stat.alpha)
        score.alpha[[n.zero]] = tmp.zero$score.alpha
        est.cov.zero[[n.zero]] = tmp.zero$est.cov.zero
        score.pvalue.alpha = append(score.pvalue.alpha, (1 - pchisq(tmp.zero$score.stat.alpha, ncol(Y0)*length(Z.par.index))))
        zero.vA.list.meta[[n.zero]] = tmp.zero$vA.list
        zero.Vinv.list.meta[[n.zero]] = tmp.zero$Vinv.list
        zero.VY.list.meta[[n.zero]] = tmp.zero$VY.list
        score.alpha.meta[idx] =  score.alpha.meta[idx] +  tmp.zero$score.alpha  #FE-Burden
        est.cov.zero.meta[idx,idx] =  est.cov.zero.meta[idx,idx] + tmp.zero$est.cov.zero #FE-SKAT
      }
    }
  }
  ############################# Asymptotic: combined
  if(n.zero>0){
    if(length(remove.index) != 0){
      Z.list = Z.list[-remove.index]
      Y0.list = Y0.list[-remove.index]
    }
    save.index.zero = which(abs(score.alpha.meta) >= 1e-7)
    n.par.save.alpha = length(save.index.zero)
    score.alpha.meta = score.alpha.meta[save.index.zero]
    est.cov.zero.meta = est.cov.zero.meta[save.index.zero,save.index.zero]
    est.cov.zero.inv = ginv(est.cov.zero.meta)
    if(!is.null(W.zero)){
      W.zero = W.zero[save.index.zero,save.index.zero]
    }
    if (Method == "FE-MV")
    {
      score.stat.zero.meta = score.alpha.meta %*% est.cov.zero.inv %*% score.alpha.meta #FE-Burden
      score.pvalue.zero = 1- pchisq(score.stat.zero.meta,df = n.par.save.alpha)
      df.zero = n.par.save.alpha
    }
    #if (Method == "SKAT")
    # {
    #   if(is.null(W.zero))
    #   {
    #     W.zero = diag(1,nrow = n.par.interest.alpha)
    #   }
    #   if(n.zero>0){
    #     eigen.cov.zero <- eigen(est.cov.zero.meta)
    #     eigen.cov.zero.sqrt = eigen.cov.zero$vectors %*% diag(sqrt(eigen.cov.zero$values),nrow = length(eigen.cov.zero$values)) %*% solve(eigen.cov.zero$vectors)
    #     weight.cov.zero = eigen(eigen.cov.zero.sqrt %*% W.zero %*% eigen.cov.zero.sqrt)$values #eign.val/sum(eign.val)
    #     score.stat.zero.meta = score.alpha.meta %*% W.zero %*% score.alpha.meta
    #     score.pvalue.zero = davies(score.stat.zero.meta,weight.cov.zero, h = rep(1,n.par.save.alpha), delta = rep(0,n.par.save.alpha), sigma = 0, lim = 10000, acc = 0.0001)$Qq
    #     score.pvalue.zero = ifelse(score.pvalue.zero>0,score.pvalue.zero,1e-7)
    #     df.zero = n.par.interest.alpha
    #   }else{
    #     score.stat.zero.meta = NA
    #     score.pvalue.zero = NA
    #     df.zero = NA
    #   }
    #   zero.results = list(score.stat = score.stat.zero.meta, score.pvalue = score.pvalue.zero, df = df.zero)
    # }
    if (Method == "FE-VC"){

      weight.cov.zero = eigen(est.cov.zero.inv)$values
      score.stat.zero.meta = score.alpha.meta %*% est.cov.zero.inv %*% est.cov.zero.inv %*% score.alpha.meta #SKAT-VC
      score.pvalue.zero = davies(score.stat.zero.meta, weight.cov.zero, h = rep(1,n.par.save.alpha), delta = rep(0,n.par.save.alpha), sigma = 0, lim = 10000, acc = 0.0001)$Qq
      score.pvalue.zero = ifelse(score.pvalue.zero>0,score.pvalue.zero,1e-7)
      df.zero = n.par.save.alpha
    }
    # if (Method == "Fisher")
    # {
    #
    #   if(n.zero>0){
    #     score.pvalue.zero = .F.test(score.pvalue.alpha)
    #     score.stat.zero.meta = score.pvalue.zero
    #   }else{
    #     score.stat.zero.meta = NA
    #     score.pvalue.zero = NA
    #   }
    #   zero.results = list(score.stat = score.stat.zero.meta, score.pvalue = score.pvalue.zero)
    # }
    # if (Method == "Mix")
    # {
    #   if(n.zero>0){
    #     score.stat.zero.meta = score.alpha.meta %*% ginv(est.cov.zero.meta) %*% score.alpha.meta #FE-Burden
    #     score.pvalue.zero = 1- pchisq(score.stat.zero.meta,df = n.par.save.alpha)
    #     df.zero = n.par.save.alpha
    #   }else{
    #     score.stat.zero.meta = NA
    #     score.pvalue.zero = NA
    #     df.zero = NA
    #   }
    #   zero.results = list(score.stat = score.stat.zero.meta, score.pvalue = score.pvalue.zero, df = df.zero)
    # }
    if(Method == "Het-SKAT"){
      # W = diag(1,nrow = n.par.interest.alpha)
      score.stat.zero.meta = 0
      for( i in 1:n.zero){
        est.inv = ginv(est.cov.zero[[i]])
        tmp  = tcrossprod(crossprod(score.alpha[[i]],est.inv), crossprod(score.alpha[[i]], est.inv))
        score.stat.zero.meta = score.stat.zero.meta + tcrossprod(crossprod(score.alpha[[i]],est.inv), crossprod(score.alpha[[i]], est.inv))
      }
    }
    if(Method == "RE-SKAT"){
      U.tau.b = 0
      a2 = 0
      a3 = 0
      a4 = 0
      for( i in 1:n.zero){
        est.inv = ginv(est.cov.zero[[i]])
        U.tau.b = U.tau.b + tcrossprod(crossprod(score.alpha[[i]],est.inv), crossprod(score.alpha[[i]], est.inv))
        a2 = a2 + sum(diag(crossprod(t(est.inv),est.inv)))
        a3 = a3 + sum(diag(crossprod(t(est.inv),est.inv)))
        a4  = a4 + sum(diag(crossprod(t(est.inv),est.inv)))
      }
      a1 = sum(diag(crossprod(t(est.cov.zero.inv),est.cov.zero.inv)))
      U.tau.b = 1/2 * U.tau.b + 1/2 * sum(diag(est.cov.zero.inv))
      U.tau.w = 1/2 * crossprod( score.alpha.meta, crossprod(t(est.cov.zero.inv), crossprod( t(est.cov.zero.inv), score.alpha.meta)))
      + 1/2 * sum(diag(est.cov.zero.inv)) #SKAT-VC
      score.stat.zero.meta = crossprod(c(U.tau.w, U.tau.b), crossprod(1/2 * matrix(c(a1,a2,a3,a4),ncol = 2), c(U.tau.w, U.tau.b)))
    }
    zero.results = list(score.stat = score.stat.zero.meta, score.pvalue = score.pvalue.zero, df = df.zero)
  }

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


        results = .resample.work.zero.meta(Z.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, score.stat.zero.meta, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, start.nperm, end.nperm, n.zero, zero.acc, Method = Method, W.zero = NULL)
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

          results = .resample.work.zero.meta(Z.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, score.stat.zero.meta, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, start.nperm, end.nperm, n.zero, zero.acc, Method = Method, W.zero = NULL)

          n.zero = results$n.zero.new
          zero.acc = results$zero.acc.new

        }

      }

      tmp = (zero.acc+1)/(n.zero+1)


    }else{

      tmp = NA
    }

    zero.results = c(zero.results, score.Rpvalue = tmp)

  }
  return(zero.results)
}


#' Title
#' @inheritParams QCAT_Meta
#' @param Z a list of matrices contains covariates for the zero-part test with each column pertains to one variable (pertains to the covariate of interest or the confounders). The number of elements of Z and OTU must be the same. The column number of each matrix in this list must be the same.
#' @param Z.index a vector indicate the columns in X for the covariate(s) of interest.
#'
#' @return A list with this elements
#'    \item{lineage.pval}{p-values for all lineages. By default ( Method = "FE-MV", n.perm = NULL ), only the asymptotic test will be performed. If the meta-analysis method is random effect methods like RE-SKAT or Het-SKAT, then resampling test must be performed.}
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
#' \doi{http://www.jstor.org/stable/2346101}
#' @references
#' Zeger, Scott L., and Kung-Yee Liang. (1986) Longitudinal Data Analysis for Discrete and Continuous Outcomes.
#' \emph{Biometrics}
#' \doi{https://doi.org/10.2307/2531248}.
QCAT_GEE_Meta <- function(OTU, Z, Z.index, Tax=NULL, Method = "FE-MV", min.depth=0, n.perm=NULL,  fdr.alpha=0.05){
  n.OTU = length(OTU)
  n.resample = n.perm
  # if(length(unique(sapply(1:n.OTU,function(j) ncol(OTU[[j]]))))!= 1){
  #   stop("The taxa in each study should be the same")
  # }

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
      stop(paste0("Samples in the OTU table and the covariate table for the zero part
                  of study ", i, " should be the same"))
    }

    remove.subject = which(rowSums(OTU[[i]])<min.depth)
    if(length(remove.subject)>0){
      print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth, "in OTU table", i))
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

  Z = lapply(1:n.OTU,function(j) cbind(1, Z[[j]])) # add the intercept term
  Z.index = Z.index + 1

  if(is.null(Tax)){ # perform one test using all OTUs
    if(is.null(n.resample)){ # asymptotic test only
      tmp = .Score.test.zero.meta(count, Z, Z.index, seed=11, resample=FALSE, n.replicates=NULL, Method = Method)
      pval.zero = as.matrix( tmp$score.pvalue )
      colnames(pval.zero) = paste0("Asymptotic-",Method)
    }else{
      tmp = .Score.test.zero.meta(count, Z, Z.index, seed=11, resample=TRUE, n.replicates=NULL, Method = Method)
      if(Method %in% c("RE-SKAT", "Het-SKAT")){
        pval.zero = c(tmp$score.Rpvalue)
        names(pval.zero) = paste0("Resampling-",Method)
      }else{
        pval.zero = c(tmp$score.pvalue, tmp$score.Rpvalue)
        names(pval.zero) = c(paste0("Asymptotic-",Method),paste0("Resampling-",Method))
      }
    }
    return(pval.zero)
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
            tmp = .Score.test.zero.meta(Y, Z, Z.index, seed=11, resample=FALSE, n.replicates=n.resample, Method = Method)
            pval.zero = cbind(pval.zero, tmp$score.pvalue)
          }
          else{

            if(Method %in% c("RE-SKAT", "Het-SKAT")){
              tmp = .Score.test.zero.meta(Y, Z, Z.index, seed=11, resample=TRUE, n.replicates=n.resample, Method = Method)
              pval.zero = cbind(pval.zero, tmp$score.Rpvalue)
            }else{
              tmp = .Score.test.zero.meta(Y, Z, Z.index, seed=11, resample=TRUE, n.replicates=n.resample, Method = Method)
              pval.zero = cbind(pval.zero, c(tmp$score.pvalue, tmp$score.Rpvalue) )
            }
          }

        }

      }# lineage loop
    }

  }# level loop


  colnames(pval.zero) = subtree

  if(is.null(n.resample)){

    rownames(pval.zero) =  paste0("Asymptotic-",Method)
    score.zero.tmp = pval.zero[1,]

  }else{
    if(Method %in% c("RE-SKAT", "Het-SKAT")){
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
    score.zero.tmp = score.zero.tmp[-index.na]
    subtree.tmp = subtree.tmp[-index.na]
  }

  #score.tmp[score.tmp==0] = 1e-4
  m.test = length(score.zero.tmp)

  # Benjamini-Hochberg FDR control
  index.p = order(score.zero.tmp)
  p.sort = sort(score.zero.tmp)
  #fdr.alpha = 0.05

  # change 04/17/2022
  reject = rep(0, m.test)
  tmp = which(p.sort<=(1:m.test)*fdr.alpha/m.test)
  if(length(tmp)>0){
    index.reject = index.p[1:max(tmp)]
    reject[index.reject] = 1
  }

  sig.lineage = subtree.tmp[reject==1]

  return( list(lineage.pval=pval.zero, sig.lineage=sig.lineage) )

}


