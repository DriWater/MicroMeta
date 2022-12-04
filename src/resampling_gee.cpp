#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::vec colSums_c(arma::mat tmp){
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
  arma::mat est_cov_meta(n_par_interest_beta,n_par_interest_beta);
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
    arma::vec index_mat  = arma::kron(arma::linspace(0,(m_beta-1) * p, m_beta),arma::ones(X_par_index.n_elem)) + arma::kron(arma::ones(m_beta), X_par_index);
    arma::uvec par_interest_index_beta = arma::conv_to< arma::uvec>::from(index_mat-1);
    arma::vec par_disinterest_index = arma::linspace(1,n_beta,n_beta);
    par_disinterest_index.elem(par_interest_index_beta) = arma::zeros(par_interest_index_beta.n_elem);
    par_disinterest_index = arma::nonzeros(par_disinterest_index);
    arma::uvec par_disinterest_index_beta = arma::conv_to< arma::uvec >::from(par_disinterest_index-1);
    int n_par_beta_interest = par_interest_index_beta.n_elem;
    arma::mat Score_reduce_beta_perm(n, n_beta);
    Score_reduce_beta_perm.zeros();
    arma::mat Hess_reduce_beta_perm(n_beta, n_beta);
    Hess_reduce_beta_perm.zeros();
    for(int j = 0; j < n; j++){
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
    arma::colvec A = colSums_c(Score_reduce_beta_perm_reorg)(idx);
    arma::mat tmp1(n_par_beta_interest, n_par_beta_interest);
    tmp1.eye();
    arma::mat B1 = arma::join_rows(tmp1,-Hess_reduce_beta_perm_reorg(arma::span(0,(n_par_beta_interest-1)),arma::span(n_par_beta_interest,(n_beta-1))) *
      arma::pinv(Hess_reduce_beta_perm_reorg(arma::span(n_par_beta_interest,(n_beta-1)),arma::span(n_par_beta_interest,(n_beta-1)))));
    arma::mat B2(n_beta, n_beta);
    B2.zeros();
    for(int j = 0; j < n; j++){
       B2 += Score_reduce_beta_perm_reorg.row(j).t() * Score_reduce_beta_perm_reorg.row(j) ;
    }
    arma::mat B = B1 * B2 * B1.t();
    arma::mat t = A.t() * arma::pinv(B) * A;
    double score_stat_beta_perm = t(0,0);
    score_beta[i] = A;
    est_cov[i] = B;
    arma::vec mat_idx  = col_index_list[i];
    arma::uvec index = arma::conv_to< arma::uvec >::from(mat_idx - 1);
    est_cov_meta.submat(index,index)  +=  B;
    score_stat_beta(index) = score_stat_beta(index) + A;
    score_stat_beta_vec(i) = score_stat_beta_perm;
  }
  return Rcpp::List::create(Rcpp::Named("score_statistics") = score_stat_beta_vec,
                            Rcpp::Named("score_beta_meta") = score_stat_beta,
                            Rcpp::Named("score_beta") = score_beta,
                            Rcpp::Named("est_cov_meta") = est_cov_meta,
                            Rcpp::Named("est_cov") = est_cov);
}

