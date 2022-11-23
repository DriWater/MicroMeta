
########################################
#                                      #
#           zero Part Model             #
#                                      #
########################################

.Pi.alpha <- function(m, p, alpha, X.i) {
  Pi.out <- rep(NA, m)

  for (j in 1:m) {
    tmp <- exp(alpha[((j - 1) * p + 1):(j * p)] %*% X.i)
    if (is.infinite(tmp)) {
      Pi.out[j] <- 1
    } else {
      Pi.out[j] <- tmp / (tmp + 1)
    }
  }


  return(Pi.out)
}

.fun.score.i.alpha <- function(alpha, data, save.list = FALSE) {
  Y <- data$Y
  Z <- data$Z

  n <- nrow(Y)
  m <- ncol(Y)
  p <- ncol(Z)

  vA.list <- list()
  Vinv.list <- list()
  VY.list <- list()

  n.alpha <- m * p

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

  if (save.list) {
    return(list(Score.alpha = Score.alpha.i, vA.list = vA.list, Vinv.list = Vinv.list, VY.list = VY.list))
  } else {
    return(Score.alpha.i)
  }
}


.Score.test.stat.zero <- function(Y0, Z, Z.par.index, cor.stru) {
  Z.reduce <- Z[, -Z.par.index, drop = FALSE]
  n <- nrow(Y0)
  m <- ncol(Y0)
  p <- ncol(Z)
  p.reduce <- ncol(Z.reduce)
  outcome <- NULL
  id <- NULL
  cova.reduce <- NULL
  for (i in 1:n) {
    outcome <- c(outcome, Y0[i, ])
    index.start <- 1
    index.end <- p

    index.start.reduce <- 1
    index.end.reduce <- p.reduce

    for (j in 1:m) {
      tmp <- rep(0, m * p.reduce)
      tmp[index.start.reduce:index.end.reduce] <- Z.reduce[i, ]
      cova.reduce <- rbind(cova.reduce, tmp)
      index.start.reduce <- index.start.reduce + p.reduce
      index.end.reduce <- index.end.reduce + p.reduce
    }


    id <- c(id, rep(i, m))
  }

  # data.full = data.frame(outcome=outcome, cova, id = id, row.names=NULL)
  data.reduce <- data.frame(outcome = outcome, cova.reduce, id = id, row.names = NULL)
  # gee.full = geeglm(outcome ~ .  - id - 1, data = data.full, id = factor(id), family="binomial", corstr= "independence")
  gee.reduce <- geeglm(outcome ~ . - id - 1, data = data.reduce, id = factor(id), family = "binomial", corstr = "independence")
  # wald.test = anova(gee.full, gee.reduce)


  ########### perform score test
  n.alpha <- m * p
  par.interest.index.alpha <- kronecker(((0:(m - 1)) * p), rep(1, length(Z.par.index))) + Z.par.index
  n.par.interest.alpha <- length(par.interest.index.alpha)
  est.reduce.alpha <- rep(NA, n.alpha)
  est.reduce.alpha[par.interest.index.alpha] <- 0
  est.reduce.alpha[-par.interest.index.alpha] <- coef(gee.reduce)
  est.reduce.scale <- gee.reduce

  data.alpha <- list(Y = Y0, Z = Z)

  tmp <- .fun.score.i.alpha(est.reduce.alpha, data.alpha, save.list = TRUE)
  Score.reduce.alpha <- tmp$Score.alpha
  # for resampling test
  vA.list <- tmp$vA.list
  Vinv.list <- tmp$Vinv.list
  VY.list <- tmp$VY.list

  Hess.reduce.alpha <- .fun.hessian.alpha(est.reduce.alpha, data.alpha)
  # re-organized the score statistics and Hessian matrix
  Score.reduce.reorg <- cbind(matrix(Score.reduce.alpha[, par.interest.index.alpha], ncol = n.par.interest.alpha), matrix(Score.reduce.alpha[, -par.interest.index.alpha], ncol = n.alpha - n.par.interest.alpha))
  Hess.reduce.reorg <- rbind(
    cbind(matrix(Hess.reduce.alpha[par.interest.index.alpha, par.interest.index.alpha], nrow = n.par.interest.alpha), matrix(Hess.reduce.alpha[par.interest.index.alpha, -par.interest.index.alpha], nrow = n.par.interest.alpha)),
    cbind(matrix(Hess.reduce.alpha[-par.interest.index.alpha, par.interest.index.alpha], nrow = n.alpha - n.par.interest.alpha), matrix(Hess.reduce.alpha[-par.interest.index.alpha, -par.interest.index.alpha], nrow = n.alpha - n.par.interest.alpha))
  )


  A <- colSums(Score.reduce.reorg)[1:n.par.interest.alpha]

  B1 <- cbind(diag(n.par.interest.alpha), -Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha + 1):n.alpha)] %*% ginv(Hess.reduce.reorg[((n.par.interest.alpha + 1):n.alpha), ((n.par.interest.alpha + 1):n.alpha)]))

  B2 <- matrix(0, n.alpha, n.alpha)
  for (i in 1:n) {
    B2 <- B2 + Score.reduce.reorg[i, ] %o% Score.reduce.reorg[i, ]
  }

  B <- B1 %*% B2 %*% t(B1)
  score.stat.alpha <- A %*% ginv(B) %*% A
  score.pvalue.alpha <- 1 - pchisq(score.stat.alpha, n.par.interest.alpha)


  return(list(score.df.alpha = n.par.interest.alpha, score.stat.alpha = score.stat.alpha, score.alpha = A, est.cov.zero = B, score.pvalue.alpha = score.pvalue.alpha, vA.list = vA.list, Vinv.list = Vinv.list, VY.list = VY.list))
}

.Score.test.stat.zero.meta.4Gresampling <- function(Z.perm.list, Z.par.index, n.par.interest.alpha, col.zero.index.list, vA.list.meta, Vinv.list.meta, VY.list.meta, Method = "MV", W.zero = NULL) {
  W <- W.zero
  iter.num <- length(Z.perm.list)
  score.stat.alpha <- NULL
  score.alpha <- NULL
  score.pvalue.alpha <- NULL
  est.cov.zero <- NULL

  for (j in 1:iter.num) {
    vA.list <- vA.list.meta[[j]]
    VY.list <- VY.list.meta[[j]]
    Vinv.list <- Vinv.list.meta[[j]]
    Z.perm <- Z.perm.list[[j]]
    p <- ncol(Z.perm.list[[j]])
    m.alpha <- length(vA.list[[1]])
    n.alpha <- m.alpha * p
    par.index.alpha <- kronecker(((0:(m.alpha - 1)) * p), rep(1, length(Z.par.index))) + Z.par.index
    n.par.alpha <- length(par.index.alpha)
    n <- nrow(Z.perm)
    Score.reduce.alpha.perm <- matrix(0, n, n.alpha)
    Hess.reduce.alpha.perm <- matrix(0, n.alpha, n.alpha)
    for (i in 1:n) {
      ###################################################
      #                                                 #
      #         alpha part: resampling Score test        #
      #                                                 #
      ###################################################
      tD.tmp <- kronecker(.diag2(vA.list[[i]]), as.matrix(Z.perm[i, ], ncol = 1))

      Score.reduce.alpha.perm[i, ] <- Score.reduce.alpha.perm[i, ] + tD.tmp %*% VY.list[[i]]

      Hess.reduce.alpha.perm <- Hess.reduce.alpha.perm + tD.tmp %*% Vinv.list[[i]] %*% t(tD.tmp)
    }

    # re-organized the score statistics and Hessian matrix
    Score.reduce.reorg <- cbind(matrix(Score.reduce.alpha.perm[, par.index.alpha], ncol = n.par.alpha), matrix(Score.reduce.alpha.perm[, -par.index.alpha], ncol = n.alpha - n.par.alpha))
    Hess.reduce.reorg <- rbind(
      cbind(matrix(Hess.reduce.alpha.perm[par.index.alpha, par.index.alpha], nrow = n.par.alpha), matrix(Hess.reduce.alpha.perm[par.index.alpha, -par.index.alpha], nrow = n.par.alpha)),
      cbind(matrix(Hess.reduce.alpha.perm[-par.index.alpha, par.index.alpha], nrow = n.alpha - n.par.alpha), matrix(Hess.reduce.alpha.perm[-par.index.alpha, -par.index.alpha], nrow = n.alpha - n.par.alpha))
    )


    A <- colSums(Score.reduce.reorg)[1:n.par.alpha]

    B1 <- cbind(diag(n.par.alpha), -Hess.reduce.reorg[(1:n.par.alpha), ((n.par.alpha + 1):n.alpha)] %*% ginv(Hess.reduce.reorg[((n.par.alpha + 1):n.alpha), ((n.par.alpha + 1):n.alpha)]))

    B2 <- matrix(0, n.alpha, n.alpha)
    for (i in 1:n) {
      B2 <- B2 + Score.reduce.reorg[i, ] %o% Score.reduce.reorg[i, ]
    }

    B <- B1 %*% B2 %*% t(B1)
    score.stat.alpha.perm <- A %*% ginv(B) %*% A
    score.stat.alpha <- append(score.stat.alpha, score.stat.alpha.perm)
    score.alpha[[j]] <- A
    est.cov.zero[[j]] <- B
    score.pvalue.alpha <- append(score.pvalue.alpha, (1 - pchisq(score.stat.alpha.perm, n.par.interest.alpha)))
  }
  score.alpha.meta <- rep(0, n.par.interest.alpha) ## A
  est.cov.meta <- matrix(0, nrow = n.par.interest.alpha, ncol = n.par.interest.alpha) ## B
  for (i in 1:iter.num)
  {
    idx <- col.zero.index.list[[i]]
    score.alpha.meta[idx] <- score.alpha.meta[idx] + score.alpha[[i]] # FE-Burden
    est.cov.meta[idx, idx] <- est.cov.meta[idx, idx] + est.cov.zero[[i]] # FE-SKAT
  }
  save.index.zero <- which(abs(score.alpha.meta) >= 1e-7)
  n.par.save.alpha <- length(save.index.zero)
  score.alpha.meta <- score.alpha.meta[save.index.zero]
  est.cov.meta <- est.cov.meta[save.index.zero, save.index.zero]
  est.cov.inv <- ginv(est.cov.meta)
  if (Method == "MV") {
    score.stat.alpha.perm <- score.alpha.meta %*% est.cov.inv %*% score.alpha.meta # FE-Burden
  }
  if (Method == "SKAT") {
    if (is.null(W)) {
      W <- diag(1, nrow = n.par.save.alpha)
    } else {
      W <- W[save.index.zero, save.index.zero]
    }
    score.stat.alpha.perm <- score.alpha.meta %*% W %*% score.alpha.meta
  }
  if (Method == "FE-VC") {
    weight.cov.inv <- eigen(est.cov.inv)$values # eign.val/sum(eign.val)
    score.stat.alpha.perm <- score.alpha.meta %*% est.cov.inv %*% est.cov.inv %*% score.alpha.meta # SKAT-VC
  }
  if (Method == "Fisher") {
    score.stat.alpha.perm <- -2 * sum(log(score.pvalue.alpha))
  }

  return(as.numeric(score.stat.alpha.perm))
}
