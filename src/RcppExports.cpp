// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// update_loadings_bin
arma::mat update_loadings_bin(const arma::mat& F_T, arma::mat& L, const arma::mat& Y_T, const arma::mat& N_T, const std::vector<int> update_indices, unsigned int num_iter, const bool line_search, const double alpha, const double beta, const double ccd_iter_tol);
RcppExport SEXP _plash_update_loadings_bin(SEXP F_TSEXP, SEXP LSEXP, SEXP Y_TSEXP, SEXP N_TSEXP, SEXP update_indicesSEXP, SEXP num_iterSEXP, SEXP line_searchSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ccd_iter_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type F_T(F_TSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y_T(Y_TSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type N_T(N_TSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type update_indices(update_indicesSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< const bool >::type line_search(line_searchSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type ccd_iter_tol(ccd_iter_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(update_loadings_bin(F_T, L, Y_T, N_T, update_indices, num_iter, line_search, alpha, beta, ccd_iter_tol));
    return rcpp_result_gen;
END_RCPP
}
// update_factors_bin
arma::mat update_factors_bin(const arma::mat& L_T, arma::mat& FF, const arma::mat& Y, const arma::mat& N, const std::vector<int> update_indices, unsigned int num_iter, const bool line_search, const double alpha, const double beta, const double ccd_iter_tol);
RcppExport SEXP _plash_update_factors_bin(SEXP L_TSEXP, SEXP FFSEXP, SEXP YSEXP, SEXP NSEXP, SEXP update_indicesSEXP, SEXP num_iterSEXP, SEXP line_searchSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ccd_iter_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L_T(L_TSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type FF(FFSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type update_indices(update_indicesSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< const bool >::type line_search(line_searchSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type ccd_iter_tol(ccd_iter_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(update_factors_bin(L_T, FF, Y, N, update_indices, num_iter, line_search, alpha, beta, ccd_iter_tol));
    return rcpp_result_gen;
END_RCPP
}
// update_loadings_missing_bin
arma::mat update_loadings_missing_bin(const arma::mat& F_T, arma::mat& L, const arma::mat& Y_T, const arma::mat& N_T, Rcpp::List nonmissing_index_list, const std::vector<int> update_indices, unsigned int num_iter, const bool line_search, const double alpha, const double beta, const double ccd_iter_tol);
RcppExport SEXP _plash_update_loadings_missing_bin(SEXP F_TSEXP, SEXP LSEXP, SEXP Y_TSEXP, SEXP N_TSEXP, SEXP nonmissing_index_listSEXP, SEXP update_indicesSEXP, SEXP num_iterSEXP, SEXP line_searchSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ccd_iter_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type F_T(F_TSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y_T(Y_TSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type N_T(N_TSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type nonmissing_index_list(nonmissing_index_listSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type update_indices(update_indicesSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< const bool >::type line_search(line_searchSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type ccd_iter_tol(ccd_iter_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(update_loadings_missing_bin(F_T, L, Y_T, N_T, nonmissing_index_list, update_indices, num_iter, line_search, alpha, beta, ccd_iter_tol));
    return rcpp_result_gen;
END_RCPP
}
// update_factors_missing_bin
arma::mat update_factors_missing_bin(const arma::mat& L_T, arma::mat& FF, const arma::mat& Y, const arma::mat& N, Rcpp::List nonmissing_index_list, const std::vector<int> update_indices, unsigned int num_iter, const bool line_search, const double alpha, const double beta, const double ccd_iter_tol);
RcppExport SEXP _plash_update_factors_missing_bin(SEXP L_TSEXP, SEXP FFSEXP, SEXP YSEXP, SEXP NSEXP, SEXP nonmissing_index_listSEXP, SEXP update_indicesSEXP, SEXP num_iterSEXP, SEXP line_searchSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ccd_iter_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L_T(L_TSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type FF(FFSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type N(NSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type nonmissing_index_list(nonmissing_index_listSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type update_indices(update_indicesSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< const bool >::type line_search(line_searchSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type ccd_iter_tol(ccd_iter_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(update_factors_missing_bin(L_T, FF, Y, N, nonmissing_index_list, update_indices, num_iter, line_search, alpha, beta, ccd_iter_tol));
    return rcpp_result_gen;
END_RCPP
}
// solve_nb_reg_cpp
arma::vec solve_nb_reg_cpp(const arma::mat X, const arma::mat X_sqrd, const arma::vec y, const arma::vec phi, arma::vec b, const std::vector<int> update_indices, unsigned int num_iter, const bool line_search, const double alpha, const double beta, const double ccd_iter_tol);
RcppExport SEXP _plash_solve_nb_reg_cpp(SEXP XSEXP, SEXP X_sqrdSEXP, SEXP ySEXP, SEXP phiSEXP, SEXP bSEXP, SEXP update_indicesSEXP, SEXP num_iterSEXP, SEXP line_searchSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ccd_iter_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type X_sqrd(X_sqrdSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type update_indices(update_indicesSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< const bool >::type line_search(line_searchSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type ccd_iter_tol(ccd_iter_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_nb_reg_cpp(X, X_sqrd, y, phi, b, update_indices, num_iter, line_search, alpha, beta, ccd_iter_tol));
    return rcpp_result_gen;
END_RCPP
}
// update_loadings
void update_loadings(const arma::mat& F_T, arma::mat& L, const arma::mat& Y_T, const arma::mat& deriv_const_mat, const std::vector<int> update_indices, unsigned int num_iter, const bool line_search, const double alpha, const double beta, const double ccd_iter_tol);
RcppExport SEXP _plash_update_loadings(SEXP F_TSEXP, SEXP LSEXP, SEXP Y_TSEXP, SEXP deriv_const_matSEXP, SEXP update_indicesSEXP, SEXP num_iterSEXP, SEXP line_searchSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ccd_iter_tolSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type F_T(F_TSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y_T(Y_TSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type deriv_const_mat(deriv_const_matSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type update_indices(update_indicesSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< const bool >::type line_search(line_searchSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type ccd_iter_tol(ccd_iter_tolSEXP);
    update_loadings(F_T, L, Y_T, deriv_const_mat, update_indices, num_iter, line_search, alpha, beta, ccd_iter_tol);
    return R_NilValue;
END_RCPP
}
// update_factors
double update_factors(const arma::mat& L_T, arma::mat& FF, const arma::mat& Y, const arma::mat& deriv_const_mat, const std::vector<int> update_indices, unsigned int num_iter, const bool line_search, const double alpha, const double beta, const double ccd_iter_tol);
RcppExport SEXP _plash_update_factors(SEXP L_TSEXP, SEXP FFSEXP, SEXP YSEXP, SEXP deriv_const_matSEXP, SEXP update_indicesSEXP, SEXP num_iterSEXP, SEXP line_searchSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ccd_iter_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L_T(L_TSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type FF(FFSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type deriv_const_mat(deriv_const_matSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type update_indices(update_indicesSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< const bool >::type line_search(line_searchSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type ccd_iter_tol(ccd_iter_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(update_factors(L_T, FF, Y, deriv_const_mat, update_indices, num_iter, line_search, alpha, beta, ccd_iter_tol));
    return rcpp_result_gen;
END_RCPP
}
// update_loadings_sp
void update_loadings_sp(const arma::mat& F_T, arma::mat& L, const arma::sp_mat& Y_T, const arma::mat& deriv_const_mat, const std::vector<int> update_indices, unsigned int num_iter, const bool line_search, const double alpha, const double beta, const double ccd_iter_tol);
RcppExport SEXP _plash_update_loadings_sp(SEXP F_TSEXP, SEXP LSEXP, SEXP Y_TSEXP, SEXP deriv_const_matSEXP, SEXP update_indicesSEXP, SEXP num_iterSEXP, SEXP line_searchSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ccd_iter_tolSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type F_T(F_TSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Y_T(Y_TSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type deriv_const_mat(deriv_const_matSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type update_indices(update_indicesSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< const bool >::type line_search(line_searchSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type ccd_iter_tol(ccd_iter_tolSEXP);
    update_loadings_sp(F_T, L, Y_T, deriv_const_mat, update_indices, num_iter, line_search, alpha, beta, ccd_iter_tol);
    return R_NilValue;
END_RCPP
}
// update_factors_sp
double update_factors_sp(const arma::mat& L_T, arma::mat& FF, const arma::sp_mat& Y, const arma::mat& deriv_const_mat, const std::vector<int> update_indices, unsigned int num_iter, const bool line_search, const double alpha, const double beta, const double ccd_iter_tol);
RcppExport SEXP _plash_update_factors_sp(SEXP L_TSEXP, SEXP FFSEXP, SEXP YSEXP, SEXP deriv_const_matSEXP, SEXP update_indicesSEXP, SEXP num_iterSEXP, SEXP line_searchSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP ccd_iter_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L_T(L_TSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type FF(FFSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type deriv_const_mat(deriv_const_matSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type update_indices(update_indicesSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type num_iter(num_iterSEXP);
    Rcpp::traits::input_parameter< const bool >::type line_search(line_searchSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type ccd_iter_tol(ccd_iter_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(update_factors_sp(L_T, FF, Y, deriv_const_mat, update_indices, num_iter, line_search, alpha, beta, ccd_iter_tol));
    return rcpp_result_gen;
END_RCPP
}
// big_exp_crossprod
double big_exp_crossprod(const arma::mat& L, const arma::mat& F, const int n, const int p);
RcppExport SEXP _plash_big_exp_crossprod(SEXP LSEXP, SEXP FSEXP, SEXP nSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type F(FSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(big_exp_crossprod(L, F, n, p));
    return rcpp_result_gen;
END_RCPP
}
// big_elementwise_mult_crossprod
double big_elementwise_mult_crossprod(const arma::mat& L, const arma::mat& F, const arma::vec& nonzero_y, const std::vector<int> nonzero_y_i_idx, const std::vector<int> nonzero_y_j_idx, const int num_nonzero_y);
RcppExport SEXP _plash_big_elementwise_mult_crossprod(SEXP LSEXP, SEXP FSEXP, SEXP nonzero_ySEXP, SEXP nonzero_y_i_idxSEXP, SEXP nonzero_y_j_idxSEXP, SEXP num_nonzero_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type F(FSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type nonzero_y(nonzero_ySEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type nonzero_y_i_idx(nonzero_y_i_idxSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type nonzero_y_j_idx(nonzero_y_j_idxSEXP);
    Rcpp::traits::input_parameter< const int >::type num_nonzero_y(num_nonzero_ySEXP);
    rcpp_result_gen = Rcpp::wrap(big_elementwise_mult_crossprod(L, F, nonzero_y, nonzero_y_i_idx, nonzero_y_j_idx, num_nonzero_y));
    return rcpp_result_gen;
END_RCPP
}
// deriv_product
arma::mat deriv_product(const arma::mat& L, const arma::mat& F);
RcppExport SEXP _plash_deriv_product(SEXP LSEXP, SEXP FSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type F(FSEXP);
    rcpp_result_gen = Rcpp::wrap(deriv_product(L, F));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_plash_update_loadings_bin", (DL_FUNC) &_plash_update_loadings_bin, 10},
    {"_plash_update_factors_bin", (DL_FUNC) &_plash_update_factors_bin, 10},
    {"_plash_update_loadings_missing_bin", (DL_FUNC) &_plash_update_loadings_missing_bin, 11},
    {"_plash_update_factors_missing_bin", (DL_FUNC) &_plash_update_factors_missing_bin, 11},
    {"_plash_solve_nb_reg_cpp", (DL_FUNC) &_plash_solve_nb_reg_cpp, 11},
    {"_plash_update_loadings", (DL_FUNC) &_plash_update_loadings, 10},
    {"_plash_update_factors", (DL_FUNC) &_plash_update_factors, 10},
    {"_plash_update_loadings_sp", (DL_FUNC) &_plash_update_loadings_sp, 10},
    {"_plash_update_factors_sp", (DL_FUNC) &_plash_update_factors_sp, 10},
    {"_plash_big_exp_crossprod", (DL_FUNC) &_plash_big_exp_crossprod, 4},
    {"_plash_big_elementwise_mult_crossprod", (DL_FUNC) &_plash_big_elementwise_mult_crossprod, 6},
    {"_plash_deriv_product", (DL_FUNC) &_plash_deriv_product, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_plash(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
