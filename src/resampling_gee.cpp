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
Rcpp::List score_test_stat_meta_resampling_c(const Rcpp::List& Z_perm_list, const Rcpp::List& col_zero_index_list, const Rcpp::List& vA_list_meta, const Rcpp::List& Vinv_list_meta,
                                             const Rcpp::List& VY_list_meta, const arma::vec& Z_par_index,int n_par_interest_alpha){
  int total_num = Z_perm_list.length();
  arma::vec score_stat_alpha = arma::vec(n_par_interest_alpha,arma::fill::zeros);
  arma::vec score_stat_alpha_vec = arma::vec(total_num);
  arma::mat est_cov_meta(n_par_interest_alpha,n_par_interest_alpha);
  est_cov_meta.zeros();
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
    arma::mat index_mat  = arma::kron(arma::linspace(0,(m_alpha-1) * p, m_alpha),arma::ones(Z_par_index.n_elem)) + arma::kron(arma::ones(m_alpha), Z_par_index);
    arma::uvec par_interest_index_alpha = arma::conv_to< arma::uvec>::from(arma::vectorise(index_mat)-1);
    arma::vec par_disinterest_index = arma::linspace(1,n_alpha,n_alpha);
    par_disinterest_index.elem(par_interest_index_alpha) = arma::zeros(par_interest_index_alpha.n_elem);
    par_disinterest_index = arma::nonzeros(par_disinterest_index) - 1;
    arma::uvec par_disinterest_index_alpha = arma::conv_to< arma::uvec >::from(par_disinterest_index);
    int n_par_alpha_interest = par_interest_index_alpha.n_elem;
    arma::mat Score_reduce_alpha_perm(n, n_alpha);
    Score_reduce_alpha_perm.zeros();
    arma::mat Hess_reduce_alpha_perm(n_alpha, n_alpha);
    Hess_reduce_alpha_perm.zeros();
    for(int j = 0; j < n; j++){
      arma::vec vec1 = vA_list[j];
      arma::mat tmp1 = arma::diagmat(vec1);
      arma::rowvec vec2 = Z_perm.row(j);
      arma::mat tD_tmp = arma::kron(tmp1, vec2.t());
      arma::colvec tmp2 = VY_list[j];
      arma::colvec vec3 = tD_tmp * tmp2;
      Score_reduce_alpha_perm.row(j) += vec3.t();
      arma::mat tmp3 = Vinv_list[j];
      Hess_reduce_alpha_perm = Hess_reduce_alpha_perm + tD_tmp * tmp3 * tD_tmp.t();
    }
    arma::mat Score_reduce_alpha_perm_reorg = arma::join_rows(Score_reduce_alpha_perm.cols(par_interest_index_alpha),Score_reduce_alpha_perm.cols(par_disinterest_index_alpha));
    arma::mat Hess_reduce_alpha_perm_reorg = arma::join_cols(arma::join_rows(Hess_reduce_alpha_perm.submat(par_interest_index_alpha,par_interest_index_alpha), Hess_reduce_alpha_perm.submat(par_interest_index_alpha,par_disinterest_index_alpha)),
                                                             arma::join_rows(Hess_reduce_alpha_perm.submat(par_disinterest_index_alpha,par_interest_index_alpha), Hess_reduce_alpha_perm.submat(par_disinterest_index_alpha,par_disinterest_index_alpha)));
    arma::uvec idx = arma::conv_to< arma::uvec >::from(arma::linspace(0, (n_par_alpha_interest-1), n_par_alpha_interest));
    arma::colvec A = colSums_c(Score_reduce_alpha_perm_reorg)(idx);
    arma::mat tmp4(n_par_alpha_interest, n_par_alpha_interest);
    tmp4.eye();
    arma::mat B1 = arma::join_rows(tmp4,-Hess_reduce_alpha_perm_reorg(arma::span(0,(n_par_alpha_interest-1)),arma::span(0,(n_par_alpha_interest-1))) *
      arma::pinv(Hess_reduce_alpha_perm_reorg(arma::span(n_par_alpha_interest,(n_alpha-1)),arma::span(n_par_alpha_interest,(n_alpha-1)))));
    arma::mat B2(n_alpha, n_alpha);
    B2.zeros();
    for(int j = 0; j < n; j++){
      B2 = B2 + Score_reduce_alpha_perm_reorg.row(j).t() * Score_reduce_alpha_perm_reorg.row(j) ;
    }
    arma::mat B = B1 * B2 * B1.t();
    arma::mat t = A.t() * arma::pinv(B) * A;
    double score_stat_alpha_perm = t(0,0);
    score_stat_alpha_vec(i) = score_stat_alpha_perm;
    score_alpha[i] = A;
    est_cov[i] = B;
    arma::vec mat_idx  = col_zero_index_list[i];
    arma::uvec index = arma::conv_to< arma::uvec >::from(mat_idx-1);
    est_cov_meta.submat(index,index)  += B;
    score_stat_alpha(index) += A;
  }
  // arma::mat est_cov_inv = arma::pinv(est_cov_meta);
  return Rcpp::List::create(Rcpp::Named("score_statistics") = score_stat_alpha_vec,
                            Rcpp::Named("score_alpha_meta") = score_stat_alpha,
                            Rcpp::Named("score_alpha") = score_alpha,
                            Rcpp::Named("est_cov_meta") = est_cov_meta,
                            Rcpp::Named("est_cov") = est_cov);
}
