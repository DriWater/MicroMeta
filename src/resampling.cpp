#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::vec colSums_c(arma::mat tmp){
  // return a column vector that is a column summation of matrix
  int cols = tmp.n_cols;
  arma::vec res(cols, arma::fill::none);
  for (int i = 0; i< cols; i++){
    res(i) = sum(tmp.col(i));
  }
  return (res);
}

// [[Rcpp::export]]
Rcpp::List score_test_stat_meta_resampling_c(const Rcpp::List& X_perm_list, const Rcpp::List& col_index_list, const Rcpp::List& S_beta_list_meta, const Rcpp::List& I_beta_list_meta,
                                             const arma::vec& X_par_index,int n_par_interest_beta){
  // the total study number
  int total_num = X_perm_list.length();
  // initialize for later use
  arma::vec score_stat_beta = arma::vec(n_par_interest_beta,arma::fill::zeros);
  // arma::vec score_stat_beta_vec = arma::vec(total_num);
  arma::mat est_cov_meta(n_par_interest_beta,n_par_interest_beta);
  // create null list to save the results for each study
  Rcpp::List score_beta(total_num);
  Rcpp::List est_cov(total_num);
  est_cov_meta.zeros();
  for(int i = 0; i < total_num; i++){
    Rcpp::List S_beta_list = S_beta_list_meta[i];
    Rcpp::List I_beta_list = I_beta_list_meta[i];
    arma::mat X_perm = X_perm_list[i];
    int p = X_perm.n_cols;
    int n = X_perm.n_rows;
    arma::vec vec_tmp = S_beta_list[0];
    int m_beta = vec_tmp.n_elem;
    int n_beta = m_beta * p;
    //  the index of parameters of interest
    arma::vec index_mat = arma::kron(arma::linspace(0,(m_beta-1) * p, m_beta),arma::ones(X_par_index.n_elem)) + arma::kron(arma::ones(m_beta), X_par_index);
    arma::uvec par_interest_index_beta = arma::conv_to< arma::uvec>::from(index_mat-1);
    arma::vec par_disinterest_index = arma::linspace(1,n_beta,n_beta);
    par_disinterest_index.elem(par_interest_index_beta) = arma::zeros(par_interest_index_beta.n_elem); // set the positions of parameters of interest to zero
    // the negative index of parameters of interest
    par_disinterest_index = arma::nonzeros(par_disinterest_index);  // extract the non-zero values (parameters of not interest)
    arma::uvec par_disinterest_index_beta = arma::conv_to< arma::uvec >::from(par_disinterest_index-1);
    int n_par_beta_interest = par_interest_index_beta.n_elem;
    arma::mat Score_reduce_beta_perm(n, n_beta);
    Score_reduce_beta_perm.zeros();
    arma::mat Hess_reduce_beta_perm(n_beta, n_beta);
    Hess_reduce_beta_perm.zeros();
    for(int j = 0; j < n; j++){
      //  Beta part: resampling Score test
      arma::rowvec vec1 = S_beta_list[j];
      arma::rowvec vec2 = arma::conv_to< arma::rowvec >::from(X_perm.row(j));
      Score_reduce_beta_perm.row(j) += arma::kron(vec1, vec2);
      arma::mat mat1 = I_beta_list[j];
      // arma::mat mat2 = X_perm.row(0).t() * X_perm.row(0);
      Hess_reduce_beta_perm += arma::kron(mat1, X_perm.row(j).t() * X_perm.row(j));
      // arma::mat mat2 = arma::kron(mat1, X_perm.row(0).t() * X_perm.row(0));
    }
    arma::mat Score_reduce_beta_perm_reorg = arma::join_rows(Score_reduce_beta_perm.cols(par_interest_index_beta),Score_reduce_beta_perm.cols(par_disinterest_index_beta));
    arma::mat Hess_reduce_beta_perm_reorg = arma::join_cols(arma::join_rows(Hess_reduce_beta_perm.submat(par_interest_index_beta,par_interest_index_beta), Hess_reduce_beta_perm.submat(par_interest_index_beta,par_disinterest_index_beta)),
                                                            arma::join_rows(Hess_reduce_beta_perm.submat(par_disinterest_index_beta,par_interest_index_beta), Hess_reduce_beta_perm.submat(par_disinterest_index_beta,par_disinterest_index_beta)));
    arma::uvec idx = arma::conv_to< arma::uvec >::from(arma::linspace(0, (n_par_beta_interest-1), n_par_beta_interest));
    // re-organized the score statistics and estimate covariance matrix based on the parameters of interest
    arma::colvec A = colSums_c(Score_reduce_beta_perm_reorg)(idx);
    // arma::mat tmp1(n_par_beta_interest, n_par_beta_interest);
    // tmp1.eye();
    // ginv(Hess.reduce.reorg[(1:n.par.interest.beta), (1:n.par.interest.beta)] - crossprod(t(Hess.reduce.reorg[(1:n.par.interest.beta), ((n.par.interest.beta + 1):n.beta)]),
    //                         crossprod(t(ginv(Hess.reduce.reorg[((n.par.interest.beta + 1):n.beta), ((n.par.interest.beta + 1):n.beta)])), Hess.reduce.reorg[((n.par.interest.beta + 1):n.beta), (1:n.par.interest.beta)])))

    arma::mat B1 = arma::pinv(Hess_reduce_beta_perm_reorg(arma::span(0,(n_par_beta_interest-1)),arma::span(0,(n_par_beta_interest-1))) - Hess_reduce_beta_perm_reorg(arma::span(0,(n_par_beta_interest-1)),arma::span(n_par_beta_interest,(n_beta-1))) *
                              arma::pinv(Hess_reduce_beta_perm_reorg(arma::span(n_par_beta_interest,(n_beta-1)), arma::span(n_par_beta_interest,(n_beta-1)))) * Hess_reduce_beta_perm_reorg(arma::span(n_par_beta_interest,(n_beta-1)),arma::span(0,(n_par_beta_interest-1))));
    // arma::mat B1 = arma::join_rows(tmp1,-Hess_reduce_beta_perm_reorg(arma::span(0,(n_par_beta_interest-1)),arma::span(n_par_beta_interest,(n_beta-1))) *
    //   arma::pinv(Hess_reduce_beta_perm_reorg(arma::span(n_par_beta_interest,(n_beta-1)),arma::span(n_par_beta_interest,(n_beta-1)))));
    arma::colvec beta_hat = B1 * A;
    // U <- Score.reduce.reorg[ ,1:n.par.interest.beta] - t(crossprod(t(Hess.reduce.reorg[(1:n.par.interest.beta), ((n.par.interest.beta + 1):n.beta)]),
    //                        crossprod(t(ginv(Hess.reduce.reorg[((n.par.interest.beta + 1):n.beta), ((n.par.interest.beta + 1):n.beta)])), t(Score.reduce.reorg[ ,((n.par.interest.beta + 1):n.beta)]))))
    arma::mat U = Score_reduce_beta_perm_reorg.cols(arma::span(0,(n_par_beta_interest-1))) - Score_reduce_beta_perm_reorg.cols(arma::span(n_par_beta_interest,(n_beta-1))) *
      arma::pinv(Hess_reduce_beta_perm_reorg(arma::span(n_par_beta_interest,(n_beta-1)), arma::span(n_par_beta_interest,(n_beta-1)))).t() * Hess_reduce_beta_perm_reorg(arma::span(0,(n_par_beta_interest-1)),arma::span(n_par_beta_interest,(n_beta-1))).t();
    arma::mat B2(n_beta, n_beta);
    B2.zeros();
    for(int j = 0; j < n; j++){
       B2 += U.row(j).t() * U.row(j) ;
    }
    arma::mat cov_beta = B1 * B2 * B1;
    // arma::mat t = A.t() * arma::pinv(B) * A;
    // double score_stat_beta_perm = t(0,0);
    score_beta[i] = beta_hat; // restore the score statistics and estimated covariance matrix in lists which are of necessity for some meta-analysis methods
    est_cov[i] = cov_beta;
    arma::vec mat_idx  = col_index_list[i];
    arma::uvec index = arma::conv_to< arma::uvec >::from(mat_idx - 1);
    est_cov_meta.submat(index,index)  += arma::pinv(cov_beta);
    score_stat_beta(index) = score_stat_beta(index) + arma::pinv(cov_beta) * beta_hat;
    // score_stat_beta_vec(i) = score_stat_beta_perm;
  }
  return Rcpp::List::create(Rcpp::Named("score_beta_meta") = score_stat_beta,
                            Rcpp::Named("score_beta") = score_beta,
                            Rcpp::Named("est_cov_meta") = est_cov_meta,
                            Rcpp::Named("est_cov") = est_cov);
}


// [[Rcpp::export]]
Rcpp::List score_test_stat_zero_meta_resampling_c(const Rcpp::List& Z_perm_list, const Rcpp::List& col_zero_index_list, const Rcpp::List& vA_list_meta, const Rcpp::List& Vinv_list_meta,
                                                  const Rcpp::List& VY_list_meta, const arma::vec& Z_par_index,int n_par_interest_alpha){
  // the total study number
  int total_num = Z_perm_list.length();
  // initialize for later use
  arma::vec score_stat_alpha = arma::vec(n_par_interest_alpha,arma::fill::zeros);
  arma::vec score_stat_alpha_vec = arma::vec(total_num);
  arma::mat est_cov_meta(n_par_interest_alpha,n_par_interest_alpha);
  est_cov_meta.zeros();
  // create null list to save the results for each study
  Rcpp::List score_alpha(total_num);
  Rcpp::List est_cov(total_num);
  for(int i = 0; i < total_num; i++){
    Rcpp::List vA_list = vA_list_meta[i];
    Rcpp::List VY_list = VY_list_meta[i];
    Rcpp::List Vinv_list = Vinv_list_meta[i];
    arma::mat Z_perm = Z_perm_list[i];
    int p = Z_perm.n_cols;
    int n = Z_perm.n_rows;
    arma::vec vec_tmp = vA_list[0];
    int m_alpha = vec_tmp.n_elem;
    int n_alpha = m_alpha * p;
    //  the index of parameters of interest
    arma::mat index_mat  = arma::kron(arma::linspace(0,(m_alpha-1) * p, m_alpha),arma::ones(Z_par_index.n_elem)) + arma::kron(arma::ones(m_alpha), Z_par_index);
    arma::uvec par_interest_index_alpha = arma::conv_to< arma::uvec>::from(arma::vectorise(index_mat)-1);
    // the negative index of parameters of interest
    arma::vec par_disinterest_index = arma::linspace(1,n_alpha,n_alpha);
    par_disinterest_index.elem(par_interest_index_alpha) = arma::zeros(par_interest_index_alpha.n_elem); // set the positions of parameters of interest to zero
    par_disinterest_index = arma::nonzeros(par_disinterest_index) - 1; // extract the non-zero values (parameters of not interest)
    arma::uvec par_disinterest_index_alpha = arma::conv_to< arma::uvec >::from(par_disinterest_index);
    int n_par_alpha_interest = par_interest_index_alpha.n_elem;
    arma::mat Score_reduce_alpha_perm(n, n_alpha);
    Score_reduce_alpha_perm.zeros();
    arma::mat Hess_reduce_alpha_perm(n_alpha, n_alpha);
    Hess_reduce_alpha_perm.zeros();
    for(int j = 0; j < n; j++){
      // alpha part: resampling Score test
      arma::vec vec1 = vA_list[j];
      arma::mat tmp1 = arma::diagmat(vec1);
      arma::rowvec vec2 = Z_perm.row(j);
      arma::mat tD_tmp = arma::kron(tmp1, vec2.t());
      arma::colvec tmp2 = VY_list[j];
      arma::colvec vec3 = tD_tmp * tmp2;
      // the permutated score statistics
      Score_reduce_alpha_perm.row(j) += vec3.t();
      arma::mat tmp3 = Vinv_list[j];
      // the permutated Hessian  matrix
      Hess_reduce_alpha_perm = Hess_reduce_alpha_perm + tD_tmp * tmp3 * tD_tmp.t();
    }
    arma::mat Score_reduce_alpha_perm_reorg = arma::join_rows(Score_reduce_alpha_perm.cols(par_interest_index_alpha),Score_reduce_alpha_perm.cols(par_disinterest_index_alpha));
    arma::mat Hess_reduce_alpha_perm_reorg = arma::join_cols(arma::join_rows(Hess_reduce_alpha_perm.submat(par_interest_index_alpha,par_interest_index_alpha), Hess_reduce_alpha_perm.submat(par_interest_index_alpha,par_disinterest_index_alpha)),
                                                             arma::join_rows(Hess_reduce_alpha_perm.submat(par_disinterest_index_alpha,par_interest_index_alpha), Hess_reduce_alpha_perm.submat(par_disinterest_index_alpha,par_disinterest_index_alpha)));
    arma::uvec idx = arma::conv_to< arma::uvec >::from(arma::linspace(0, (n_par_alpha_interest-1), n_par_alpha_interest));
    // re-organized the score statistics and estimate covariance matrix based on the parameters of interest
    arma::colvec A = colSums_c(Score_reduce_alpha_perm_reorg)(idx);
    arma::mat B1 = arma::pinv(Hess_reduce_alpha_perm_reorg(arma::span(0,(n_par_alpha_interest-1)),arma::span(0,(n_par_alpha_interest-1))) - Hess_reduce_alpha_perm_reorg(arma::span(0,(n_par_alpha_interest-1)),arma::span(n_par_alpha_interest,(n_alpha-1))) *
      arma::pinv(Hess_reduce_alpha_perm_reorg(arma::span(n_par_alpha_interest,(n_alpha-1)), arma::span(n_par_alpha_interest,(n_alpha-1)))) * Hess_reduce_alpha_perm_reorg(arma::span(n_par_alpha_interest,(n_alpha-1)),arma::span(0,(n_par_alpha_interest-1))));
    // arma::mat B1 = arma::join_rows(tmp1,-Hess_reduce_beta_perm_reorg(arma::span(0,(n_par_beta_interest-1)),arma::span(n_par_beta_interest,(n_beta-1))) *
    //   arma::pinv(Hess_reduce_beta_perm_reorg(arma::span(n_par_beta_interest,(n_beta-1)),arma::span(n_par_beta_interest,(n_beta-1)))));
    arma::colvec beta_hat = B1 * A;
    // U <- Score.reduce.reorg[ ,1:n.par.interest.beta] - t(crossprod(t(Hess.reduce.reorg[(1:n.par.interest.beta), ((n.par.interest.beta + 1):n.beta)]),
    //                        crossprod(t(ginv(Hess.reduce.reorg[((n.par.interest.beta + 1):n.beta), ((n.par.interest.beta + 1):n.beta)])), t(Score.reduce.reorg[ ,((n.par.interest.beta + 1):n.beta)]))))
    arma::mat U = Score_reduce_alpha_perm_reorg.cols(arma::span(0,(n_par_alpha_interest-1))) - Score_reduce_alpha_perm_reorg.cols(arma::span(n_par_alpha_interest,(n_alpha-1))) *
      arma::pinv(Hess_reduce_alpha_perm_reorg(arma::span(n_par_alpha_interest,(n_alpha-1)), arma::span(n_par_alpha_interest,(n_alpha-1)))).t() * Hess_reduce_alpha_perm_reorg(arma::span(0,(n_par_alpha_interest-1)),arma::span(n_par_alpha_interest,(n_alpha-1))).t();
    arma::mat B2(n_alpha, n_alpha);
    B2.zeros();
    for(int j = 0; j < n; j++){
      B2 += U.row(j).t() * U.row(j) ;
    }
    arma::mat cov_beta = B1 * B2 * B1;
    score_alpha[i] = beta_hat; // restore the score statistics and estimated covariance matrix in lists which are of necessity for some meta-analysis methods
    est_cov[i] = cov_beta;
    arma::vec mat_idx  = col_zero_index_list[i];
    arma::uvec index = arma::conv_to< arma::uvec >::from(mat_idx - 1);
    est_cov_meta.submat(index,index)  += arma::pinv(cov_beta);
    score_stat_alpha(index) = score_stat_alpha(index) + arma::pinv(cov_beta) * beta_hat;
  }
  // arma::mat est_cov_inv = arma::pinv(est_cov_meta);
  return Rcpp::List::create(Rcpp::Named("score_alpha_meta") = score_stat_alpha,
                            Rcpp::Named("score_alpha") = score_alpha,
                            Rcpp::Named("est_cov_meta") = est_cov_meta,
                            Rcpp::Named("est_cov") = est_cov);
}

