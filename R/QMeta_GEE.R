
########################################
#                                      #
#       zero Part Model(GEE)           #
#                                      #
########################################
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

.Pi.alpha <- function(m, p, alpha, X.i) {
  Pi.out <- rep(NA, m)
  # calculate the exponential of alpha times X
  for (j in 1:m) {
    tmp <- exp(alpha[((j - 1) * p + 1):(j * p)] %*% X.i)
    if (is.infinite(tmp)) {
      Pi.out[j] <- 1
    } else {
      Pi.out[j] <- tmp / (tmp + 1)
    }
  }

  # no need for base in GEE method
  return(Pi.out)
}

.fun.score.i.alpha <- function(alpha, data, save.list = FALSE) {
  # return  the negative log likelihood of estimated pi for each taxa
  Y <- data$Y
  Z <- data$Z

  n <- nrow(Y)
  m <- ncol(Y)
  p <- ncol(Z)

  vA.list <- list()
  Vinv.list <- list()
  VY.list <- list()

  n.alpha <- m * p
  # check the dimension of alpha
  if (length(alpha) != n.alpha) {
    warning("Dim of initial alpha does not match the dim of covariates")
  } else {
    Score.alpha.i <- matrix(0, n, n.alpha)
    nY <- rowSums(Y)

    for (i in 1:n) {
      Pi.i <- .Pi.alpha(m, p, alpha, Z[i, ])
      vA.tmp <- Pi.i * (1 - Pi.i)
      A.i <- .diag2(vA.tmp)
      t.D.i <- kronecker(A.i, as.matrix(Z[i, ], ncol = 1))
      V.i <- A.i # independent cor structure

      tmp.V.i <- ginv(V.i)
      tmp.VY <- tmp.V.i %*% (Y[i, ] - Pi.i)
      Score.alpha.i[i, ] <- t.D.i %*% tmp.VY

      if (save.list) {
        vA.list[[i]] <- vA.tmp
        Vinv.list[[i]] <- tmp.V.i
        VY.list[[i]] <- tmp.VY
      }
    }
  }
  # save list for later permutation method
  if (save.list) {
    return(list(Score.alpha = Score.alpha.i, vA.list = vA.list, Vinv.list = Vinv.list, VY.list = VY.list))
  } else {
    return(Score.alpha.i)
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

  B1 = cbind(diag(n.par.interest.alpha), -Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha+1):n.alpha)] %*% ginv(Hess.reduce.reorg[((n.par.interest.alpha+1):n.alpha), ((n.par.interest.alpha+1):n.alpha)]) )

  B2 =  matrix(0, n.alpha, n.alpha)
  for(i in 1:n){
    B2 = B2 + Score.reduce.reorg[i,] %o% Score.reduce.reorg[i,]
  }

  B = B1 %*% B2 %*% t(B1)
  score.stat.alpha = A %*% ginv(B) %*% A
  score.pvalue.alpha = 1 - pchisq(score.stat.alpha, n.par.interest.alpha)


  return(list(score.df.alpha=n.par.interest.alpha, score.stat.alpha = score.stat.alpha, score.alpha = A, est.cov.zero=B, score.pvalue.alpha=score.pvalue.alpha, vA.list=vA.list, Vinv.list=Vinv.list, VY.list=VY.list )   )

}

.Score.test.stat.zero.meta.4Gresampling <- function(Z.perm.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, vA.list.meta, Vinv.list.meta, VY.list.meta, Method = "MV", W.zero = NULL){

  W = W.zero
  iter.num = length(Z.perm.list)
  score.stat.alpha = NULL
  score.alpha = NULL
  score.pvalue.alpha = NULL
  est.cov.zero = NULL

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

      Score.reduce.alpha.perm[i,] = Score.reduce.alpha.perm[i,] + tD.tmp %*% VY.list[[i]]

      Hess.reduce.alpha.perm = Hess.reduce.alpha.perm + tD.tmp %*% Vinv.list[[i]] %*% t(tD.tmp)


    }

    # re-organized the score statistics and Hessian matrix
    Score.reduce.reorg = cbind( matrix(Score.reduce.alpha.perm[,par.index.alpha], ncol=n.par.alpha), matrix(Score.reduce.alpha.perm[,-par.index.alpha], ncol=n.alpha - n.par.alpha) )
    Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.alpha.perm[par.index.alpha, par.index.alpha], nrow=n.par.alpha), matrix(Hess.reduce.alpha.perm[par.index.alpha, -par.index.alpha], nrow=n.par.alpha) ),
                              cbind( matrix(Hess.reduce.alpha.perm[-par.index.alpha, par.index.alpha], nrow=n.alpha - n.par.alpha), matrix(Hess.reduce.alpha.perm[-par.index.alpha, -par.index.alpha], nrow= n.alpha - n.par.alpha)))


    A = colSums(Score.reduce.reorg)[1:n.par.alpha]

    B1 = cbind(diag(n.par.alpha), -Hess.reduce.reorg[(1:n.par.alpha), ((n.par.alpha+1):n.alpha)] %*% ginv(Hess.reduce.reorg[((n.par.alpha+1):n.alpha), ((n.par.alpha+1):n.alpha)]) )

    B2 =  matrix(0, n.alpha, n.alpha)
    for(i in 1:n){
      B2 = B2 + Score.reduce.reorg[i,] %o% Score.reduce.reorg[i,]
    }

    B = B1 %*% B2 %*% t(B1)
    score.stat.alpha.perm = A %*% ginv(B) %*% A
    score.stat.alpha = append(score.stat.alpha, score.stat.alpha.perm)
    score.alpha[[j]] = A
    est.cov.zero[[j]] = B
    score.pvalue.alpha = append(score.pvalue.alpha, (1 - pchisq(score.stat.alpha.perm, n.par.interest.alpha)))

  }
  score.alpha.meta  = rep(0,n.par.interest.alpha) ## A
  est.cov.meta = matrix(0, nrow = n.par.interest.alpha, ncol = n.par.interest.alpha) ## B
  for(i in 1:iter.num)
  {
    idx = col.zero.index.list[[i]]
    score.alpha.meta[idx] =  score.alpha.meta[idx] +  score.alpha[[i]]  #FE-Burden
    est.cov.meta[idx, idx] =  est.cov.meta[idx, idx] + est.cov.zero[[i]] #FE-SKAT
  }
  save.index.zero = which(abs(score.alpha.meta) >= 1e-7)
  n.par.save.alpha = length(save.index.zero)
  score.alpha.meta =  score.alpha.meta[save.index.zero]
  est.cov.meta =  est.cov.meta[save.index.zero,save.index.zero]
  est.cov.inv = ginv(est.cov.meta)
  weight.cov.zero = NULL
  if (Method == "MV")
  {
    score.stat.alpha.perm = score.alpha.meta %*% est.cov.inv %*% score.alpha.meta #FE-Burden
  }
  if (Method == "SKAT")
  {
    if(is.null(W))
    {
      W = diag(1,nrow = n.par.save.alpha)
    }
    else{
      W = W[save.index.zero,save.index.zero]
    }
    eigen.cov.zero <- eigen(est.cov.meta)
    eigen.cov.zero.sqrt = eigen.cov.zero$vectors %*% diag(sqrt(eigen.cov.zero$values),nrow = length(eigen.cov.zero$values)) %*% solve(eigen.cov.zero$vectors)
    weight.cov.zero = eigen(eigen.cov.zero.sqrt %*% W %*% eigen.cov.zero.sqrt)$values
    score.stat.alpha.perm = score.alpha.meta %*% W %*% score.alpha.meta
  }
  if (Method == "FE-VC")
  {
    weight.cov.zero = eigen(est.cov.inv)$values #eign.val/sum(eign.val)
    score.stat.alpha.perm = score.alpha.meta %*% est.cov.inv %*% est.cov.inv %*% score.alpha.meta #SKAT-VC
  }
  if (Method == "Fisher")
  {
    score.stat.alpha.perm = -2 * sum(log(score.pvalue.alpha))
  }

  return(list(score.stat.alpha.perm = as.numeric(score.stat.alpha.perm),weight.cov.zero = weight.cov.zero))

}

# add 10/22/2022 for adaptive resampling
.resample.work.zero.meta <- function(Z.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, score.stat.zero.meta, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, start.nperm, end.nperm, n.zero, zero.acc, Method = "MV", W.zero = NULL){

  iter.num = length(Z.list)
  n.zero.new = n.zero
  zero.acc.new = zero.acc

  for(k in start.nperm:end.nperm){

    perm.index = lapply(Z.list, function(i) sample(1:nrow(i)))

    Z.perm.list = Z.list
    for (i in 1:iter.num){
      Z.perm.list[[i]][,Z.par.index] =  Z.perm.list[[i]][perm.index[[i]],Z.par.index]
    }

    tmp = try( .Score.test.stat.zero.meta.4Gresampling(Z.perm.list, Z.par.index,n.par.interest.alpha, col.zero.index.list, zero.vA.list.meta, zero.Vinv.list.meta, zero.VY.list.meta, Method = Method, W.zero = W.zero) )



    if(class(tmp) != "try-error"){

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

.Score.test.zero.meta <- function(Y.list, Z.list, Z.par.index, seed=11, resample=FALSE, n.replicates=NULL, Method = "MV", Weight.zero=NULL){
  iter.num = length(X.list)
  W.zero = Weight.zero
  m = ncol(Y.list[[1]])
  p.zero = ncol(Z.list[[1]])
  n.OTU = length(Y.list)
  n.par.interest.alpha = m*length(Z.par.index)
  if(! Method %in% c('Fisher',"MV","FE-VC",'SKAT','Mix')){
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

    score.stat.alpha = NULL
    score.alpha = NULL
    est.cov.zero = NULL
    score.pvalue.alpha = NULL
    zero.vA.list.meta = NULL
    zero.Vinv.list.meta = NULL
    zero.VY.list.meta = NULL
    col.zero.index.list = NULL

  }else{
    score.stat.alpha = NULL
    score.alpha = NULL
    est.cov.zero = NULL
    score.pvalue.alpha = NULL
    zero.vA.list.meta = list()
    zero.Vinv.list.meta = list()
    zero.VY.list.meta = list()
    col.zero.index.list = list()
    for(i in 1:iter.num){
      Y0 = Y0.list[[i]]
      Z = Z.list[[i]]
      col.zero.index = which(apply(Y0, 2, function(x) length(table(x)) ) > 1)
      Y0 = Y0[, col.zero.index , drop=FALSE]
      if(ncol(Y0)<=1)
      {
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
      if(class(tmp.zero) == "try-error"){

        next

      }else{
        n.zero = n.zero + 1
        # par.index.pos.list[[n.pos]] = kronecker(((0:(m-2))*p), rep(1,length(X1.par.index))) + X1.par.index
        #par.index.zero.list[[n.zero]] = kronecker((col.zero.index - 1)*p.zero, rep(1,length(Z.par.index))) + Z.par.index
        col.zero.index.list[[n.zero]] = col.zero.index
        score.stat.alpha = append(score.stat.alpha, tmp.zero$score.stat.alpha)
        score.alpha[[n.zero]] = tmp.zero$score.alpha
        est.cov.zero[[n.zero]] = tmp.zero$est.cov.zero
        score.pvalue.alpha = append(score.pvalue.alpha, (1 - pchisq(tmp.zero$score.stat.alpha, ncol(Y0)*length(Z.par.index))))
        zero.vA.list.meta[[n.zero]] = tmp.zero$vA.list
        zero.Vinv.list.meta[[n.zero]] = tmp.zero$Vinv.list
        zero.VY.list.meta[[n.zero]] = tmp.zero$VY.list
      }
    }
  }
  ############################# Asymptotic: combined
  if(n.zero==0){
    score.pvalue = NA
  }else{
    score.alpha.meta  = rep(0,n.par.interest.alpha) ## A
    est.cov.zero.meta = matrix(0, nrow = n.par.interest.alpha, ncol = n.par.interest.alpha) ## B
    if(n.zero > 0)
    {
      for(i in 1:n.zero)
      {
        idx = col.zero.index.list[[i]]
        score.alpha.meta[idx] =  score.alpha.meta[idx] +  score.alpha[[i]]
        est.cov.zero.meta[idx, idx] =  est.cov.zero.meta[idx,idx] + est.cov.zero[[i]]
      }
    }
    save.index.zero = which(abs(score.alpha.meta) >= 1e-7)
    n.par.save.alpha = length(save.index.zero)
    score.alpha.meta = score.alpha.meta[save.index.zero]
    est.cov.zero.meta = est.cov.zero.meta[save.index.zero,save.index.zero]
    if(!is.null(W.zero)){
      W.zero = W.zero[save.index.zero,save.index.zero]
    }
    if (Method == "MV")
    {
      if(n.zero>0){
        est.cov.zero.inv = ginv(est.cov.zero.meta)
        score.stat.zero.meta = score.alpha.meta %*% est.cov.zero.inv %*% score.alpha.meta #FE-Burden
        score.zero.pvalue = 1- pchisq(score.stat.zero.meta,df = n.par.save.alpha)
        df.zero = n.par.save.alpha
      }else{
        score.stat.zero.meta = NA
        score.zero.pvalue = NA
        df.zero = NA
      }
      zero.results = list(score.stat = score.stat.zero.meta, score.pvalue = score.zero.pvalue, df = df.zero)
    }
    if (Method == "SKAT")
    {
      if(is.null(W.zero))
      {
        W.zero = diag(1,nrow = n.par.interest.alpha)
      }
      if(n.zero>0){
        eigen.cov.zero <- eigen(est.cov.zero.meta)
        eigen.cov.zero.sqrt = eigen.cov.zero$vectors %*% diag(sqrt(eigen.cov.zero$values),nrow = length(eigen.cov.zero$values)) %*% solve(eigen.cov.zero$vectors)
        weight.cov.zero = eigen(eigen.cov.zero.sqrt %*% W.zero %*% eigen.cov.zero.sqrt)$values #eign.val/sum(eign.val)
        score.stat.zero.meta = score.alpha.meta %*% W.zero %*% score.alpha.meta
        score.zero.pvalue = davies(score.stat.zero.meta,weight.cov.zero, h = rep(1,n.par.save.alpha), delta = rep(0,n.par.save.alpha), sigma = 0, lim = 10000, acc = 0.0001)$Qq
        score.zero.pvalue = ifelse(score.zero.pvalue>0,score.zero.pvalue,1e-7)
        df.zero = n.par.interest.alpha
      }else{
        score.stat.zero.meta = NA
        score.zero.pvalue = NA
        df.zero = NA
      }
      zero.results = list(score.stat = score.stat.zero.meta, score.pvalue = score.zero.pvalue, df = df.zero)
    }
    if (Method == "FE-VC"){
      if(n.zero>0){
        est.cov.zero.inv = ginv(est.cov.zero.meta)
        weight.cov.zero = eigen(est.cov.zero.inv)$values
        score.stat.zero.meta = score.alpha.meta %*% est.cov.zero.inv %*% est.cov.zero.inv %*% score.alpha.meta #SKAT-VC
        score.zero.pvalue = davies(score.stat.zero.meta, weight.cov.zero, h = rep(1,n.par.save.alpha), delta = rep(0,n.par.save.alpha), sigma = 0, lim = 10000, acc = 0.0001)$Qq
        score.zero.pvalue = ifelse(score.zero.pvalue>0,score.zero.pvalue,1e-7)
        df.zero = n.par.save.alpha
      }else{
        score.stat.zero.meta = NA
        score.zero.pvalue = NA
        df.zero = NA
      }
      zero.results = list(score.stat = score.stat.zero.meta, score.pvalue = score.zero.pvalue, df = df.zero)
    }
    if (Method == "Fisher")
    {

      if(n.zero>0){
        score.zero.pvalue = .F.test(score.pvalue.alpha)
        score.stat.zero.meta = score.zero.pvalue
      }else{
        score.stat.zero.meta = NA
        score.zero.pvalue = NA
      }
      zero.results = list(score.stat = score.stat.zero.meta, score.pvalue = score.zero.pvalue)
    }
    if (Method == "Mix")
    {
      if(n.zero>0){
        score.stat.zero.meta = score.alpha.meta %*% ginv(est.cov.zero.meta) %*% score.alpha.meta #FE-Burden
        score.zero.pvalue = 1- pchisq(score.stat.zero.meta,df = n.par.save.alpha)
        df.zero = n.par.save.alpha
      }else{
        score.stat.zero.meta = NA
        score.zero.pvalue = NA
        df.zero = NA
      }
      zero.results = list(score.stat = score.stat.zero.meta, score.pvalue = score.zero.pvalue, df = df.zero)
    }
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

      score.Rpvalue.zero = (zero.acc+1)/(n.zero+1)
      zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)


    }else{

      zero.results = c(zero.results, score.Rpvalue = NA)
    }



  }
  return(zero.results = zero.results)
}
