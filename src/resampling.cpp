#include <Rcpp.h>
using namespace Rcpp;


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::vec colSums(arma::mat tmp){
  int cols = tmp.n_cols;
  arma::vec res(cols, arma::fill::none);
  for (size_t i = 0; i< cols; i++){
    res(i) = sum(tmp.col(i));
  }
  return (res);
}

// [[Rcpp::export]]
Rcpp::List score_test_stat_meta_resampling_c(const Rcpp::List& X_perm_list, const Rcpp::List& col_index_list, const Rcpp::List& S_beta_list_meta, const Rcpp::List& I_beta_list_meta,
                                             const arma::vec& X_par_index,int n_par_interest_beta){
  int total_num = X_perm_list.length();
  arma::vec score_stat_beta = arma::vec(n_par_interest_beta,arma::fill::zeros);
  arma::vec score_stat_beta_vec = arma::vec(total_num);
  arma::mat est_cov_meta = arma::zeros(n_par_interest_beta, n_par_interest_beta);
  // Rcpp::List score_beta = Rcpp::List::create();
  // Rcpp::List est_cov = Rcpp::List::create();
  for(int i = 0; i < total_num; i++){
    Rcpp::List S_beta_list = S_beta_list_meta[i];
    Rcpp::List I_beta_list = I_beta_list_meta[i];
    arma::mat X_perm = X_perm_list[i];
    int p = X_perm.n_cols;
    int n = X_perm.n_rows;
    int m_beta = S_beta_list.length();
    int n_beta = m_beta * p;
    arma::mat index_mat  = arma::kron(arma::linspace(0,(m_beta-1) * p, m_beta),arma::ones(X_par_index.n_elem)) + arma::kron(arma::ones(m_beta), X_par_index);
    arma::uvec par_interest_index_beta = arma::conv_to< arma::uvec>::from(arma::vectorise(index_mat)-1);
    arma::vec par_disinterest_index = arma::linspace(1,n_beta,n_beta);
    par_disinterest_index.elem(par_interest_index_beta) = arma::zeros(par_interest_index_beta.n_elem);
    par_disinterest_index = arma::nonzeros(par_disinterest_index) - 1;
    arma::uvec par_disinterest_index_beta = arma::conv_to< arma::uvec >::from(par_disinterest_index);
    int n_par_beta_interest = par_interest_index_beta.n_elem;
    arma::mat Score_reduce_beta_perm = arma::zeros(n, n_beta);
    arma::mat Hess_reduce_beta_perm = arma::zeros(n_beta, n_beta);
    for(int j = 0; j < n; j++){
      arma::vec vec1 = S_beta_list[j];
      arma::vec vec2 = X_perm.row(j);
      arma::colvec tmp = Score_reduce_beta_perm.col(j) + arma::kron(vec1, vec2);
      Score_reduce_beta_perm.col(j) = tmp;
      arma::mat mat1 = I_beta_list[j];
      Hess_reduce_beta_perm = Hess_reduce_beta_perm + arma::kron(mat1, arma::cross(X_perm.row(i),X_perm.row(i)));
    }
    arma::mat Score_reduce_beta_perm_reorg = arma::join_cols(Score_reduce_beta_perm.cols(par_interest_index_beta),Score_reduce_beta_perm.cols(par_disinterest_index_beta));
    arma::mat Hess_reduce_beta_perm_reorg = arma::join_rows(arma::join_cols(Hess_reduce_beta_perm.submat(par_interest_index_beta,par_interest_index_beta), Hess_reduce_beta_perm.submat(par_interest_index_beta,par_disinterest_index_beta)),
                                                            arma::join_cols(Hess_reduce_beta_perm.submat(par_disinterest_index_beta,par_interest_index_beta), Hess_reduce_beta_perm.submat(par_disinterest_index_beta,par_disinterest_index_beta)));
    arma::uvec idx = arma::conv_to< arma::uvec >::from(arma::linspace(0, (n_par_beta_interest-1), n_par_beta_interest));
    arma::colvec A = colSums(Score_reduce_beta_perm_reorg)(idx);
    // cbind(diag(n.par.interest.beta), -Hess.reduce.beta.perm.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.beta.perm.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
    arma::mat tmp;
    tmp.eye(n_par_beta_interest, n_par_beta_interest);
    arma::mat B1 = arma::join_cols(tmp,-Hess_reduce_beta_perm_reorg(arma::span(0,(n_par_beta_interest-1)),arma::span(0,(n_par_beta_interest-1))) *
      arma::pinv(Hess_reduce_beta_perm_reorg(arma::span(n_par_beta_interest,(n_beta-1)),arma::span(n_par_beta_interest,(n_beta-1)))));
    arma::mat B2 = arma::zeros(n_beta, n_beta);
    for(int j = 0; j < n; j++){
      arma::vec tmp1 = Score_reduce_beta_perm_reorg.row(j);
      arma::mat tmp = arma::kron(tmp1,tmp1);
      tmp.reshape(Score_reduce_beta_perm_reorg.row(j).n_elem,Score_reduce_beta_perm_reorg.row(j).n_elem);
      B2 = B2 + tmp;
    }
    arma::mat B = B1 * B2 * B1.t();
    arma::mat t = A.t() * arma::pinv(B) * A;
    double score_stat_beta_perm = t(0,0);
    score_stat_beta_vec(i) = score_stat_beta_perm;
    // score_beta[i] = A;
    // est_cov[i] = B;
    arma::vec mat_idx  = col_index_list[i];
    arma::uvec index = arma::conv_to< arma::uvec >::from(mat_idx);
    est_cov_meta.submat(index,index)  = est_cov_meta.submat(index,index) + B;
    score_stat_beta(index) = score_stat_beta(index) + A;
  }
  arma::mat est_cov_inv = arma::pinv(est_cov_meta);
  return Rcpp::List::create(Rcpp::Named("score_statistics") = score_stat_beta_vec,
                            Rcpp::Named("score_stat_beta") = score_stat_beta,
                            // Rcpp::Named("est_cov_meta")  = est_cov_meta,
                            Rcpp::Named("est_cov_inv") = est_cov_inv);
}

// .Score.test.stat.meta.4Gresampling <- function(X.perm.list, X.par.index, n.par.interest.beta, col.index.list, S.beta.list.meta, I.beta.list.meta, Method = "MV", W = NULL){
//
// stu.num = length(X.perm.list)
//   score.stat.beta = NULL
//   score.beta = NULL
//   est.cov = NULL
//   score.pvalue.beta = NULL
//
//   for(j in c(1:stu.num)){
//     S.beta.list = S.beta.list.meta[[j]]
//     I.beta.list = I.beta.list.meta[[j]]
//     X.perm = X.perm.list[[j]]
//     p = ncol(X.perm)
//     n = nrow(X.perm)
//     m.beta = length(S.beta.list[[1]])
// #n.par.interest.beta = m.beta
//     n.beta = m.beta*p
// #n.beta = m.beta*2
//     par.interest.index.beta =  kronecker( ((0:(m.beta-1))*p), rep(1,length(X.par.index))) + X.par.index
// #par.interest.index.beta = (1:m.beta )*2
//     n.par.beta.interest = length(par.interest.index.beta)
//     Score.reduce.beta.perm = matrix(0, n, n.beta )
//     Hess.reduce.beta.perm = matrix(0, n.beta, n.beta )
//     for(i in 1:n){
//
// ###################################################
// #                                                 #
// #         Beta part: resampling Score test        #
// #                                                 #
// ###################################################
//       Score.reduce.beta.perm[i,] = Score.reduce.beta.perm[i,] + kronecker(matrix(S.beta.list[[i]], ncol=1),  matrix(X.perm[i,], ncol=1))
//
//         Hess.reduce.beta.perm = Hess.reduce.beta.perm + kronecker(I.beta.list[[i]], (  X.perm[i,] %o% X.perm[i,] ) )
// #     if(sum(is.na(Hess.reduce.beta.perm))>0){
// #       print(i); break;
// #
// #     }
//     }
// ###################################################
// #                                                 #
// #         Beta part: resampling Score test        #
// #                                                 #
// ###################################################
//     Score.reduce.beta.perm.reorg = cbind( matrix(Score.reduce.beta.perm[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta.perm[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
//       Hess.reduce.beta.perm.reorg = rbind(cbind( matrix(Hess.reduce.beta.perm[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta.perm[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ),
//                                           cbind( matrix(Hess.reduce.beta.perm[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta.perm[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
//
//
// # re-organized the score statistics and Hessian matrix according to parameter of interest
//       A = colSums(Score.reduce.beta.perm.reorg)[1:n.par.interest.beta]
//
//       B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.beta.perm.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.beta.perm.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
//
//         B2 =  matrix(0, n.beta, n.beta)
//
//         for(i in 1:n){
//           B2 = B2 + Score.reduce.beta.perm.reorg[i,] %o% Score.reduce.beta.perm.reorg[i,]
//         }
//
//         B = B1 %*% B2 %*% t(B1)
//           score.stat.beta.perm = A %*% ginv(B) %*% A
//           score.stat.beta = append(score.stat.beta, score.stat.beta.perm)
//           score.beta[[j]] = A
//           est.cov[[j]] = B
//           score.pvalue.beta = append(score.pvalue.beta, (1 - pchisq(score.stat.beta.perm, n.par.interest.beta)))
//
//   }
//   score.beta.meta  = rep(0,n.par.interest.beta) ## A
//     est.cov.meta = matrix(0, nrow = n.par.interest.beta, ncol = n.par.interest.beta) ## B
//     for(i in 1:stu.num)
//     {
//       idx = col.index.list[[i]]
//       score.beta.meta[idx] =  score.beta.meta[idx] +  score.beta[[i]]  #FE-Burden
//       est.cov.meta[idx,idx] =  est.cov.meta[idx,idx] + est.cov[[i]] #FE-SKAT
//     }
//     est.cov.inv = ginv(est.cov.meta)
//       if (Method == "MV")
//       {
//         score.stat.meta.perm = score.beta.meta %*% est.cov.inv %*% score.beta.meta #FE-Burden
//       }
//       if (Method == "SKAT")
//       {
//         if(is.null(W))
//         {
//           W = diag(1,nrow = n.par.interest.beta)
//         }
// #fesk.p = farebrother(score.stat.fesk,weight, h = rep(1,m-1), delta = rep(0,m-1), maxit = 100000,eps = 10^(-10), mode = 1)$Qq
//         score.stat.meta.perm = score.beta.meta %*% W %*% score.beta.meta
//       }
//       if (Method == "VC")
//       {
//         weight.cov.inv = eigen(est.cov.inv)$values #eign.val/sum(eign.val)
//         score.stat.meta.perm = score.beta.meta %*% est.cov.inv %*% est.cov.inv %*% score.beta.meta #SKAT-VC
//       }
//       if (Method == "Fisher")
//       {
//         score.stat.meta.perm = -2 * sum(log(score.pvalue.beta))
//       }
//
//
//       return(as.numeric(score.stat.meta.perm))
//
// }





