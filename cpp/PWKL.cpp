#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
int vecmaxInd(NumericVector x) {
  double m = max(na_omit(x));
  int j = NA_REAL;
  for( int i = 0; i<x.size(); ++i){
    if( x[i] == m) {
      j = i;
    }
  }
  return j;
}

// [[Rcpp::export]]
NumericVector PWKL_Score_cpp(NumericMatrix theta, IntegerVector unasked, NumericVector lpYit_pre){
  int C = theta.nrow();
  int Yit = vecmaxInd(lpYit_pre);
  int remaining = unasked.size();
  NumericVector Score (remaining);
  double pXijq1Yi = NA_REAL;
  double pXijq1Yit = NA_REAL;
  double pXijq0Yi = NA_REAL;
  double pXijq0Yit = NA_REAL;
  NumericMatrix Dj( C , remaining );
  for(unsigned int c = 0; c < C; ++c) {
    for(unsigned int j = 0; j < remaining; ++j) {
      pXijq1Yi = theta(c,unasked[j]-1);
      pXijq1Yit = theta(Yit,unasked[j]-1);
      pXijq0Yi = 1-pXijq1Yi;
      pXijq0Yit = 1-pXijq1Yit;
      Dj(c, j) = ( log(pXijq0Yit) - log(pXijq0Yi) ) * pXijq0Yit + ( log(pXijq1Yit) - log(pXijq1Yi) ) * pXijq1Yit;
    }
  }
  for(unsigned int j = 0; j < remaining; ++j) {
      for(unsigned int c = 0; c < C; ++c) {
        Score[j] = Score[j] + Dj(c,j)*exp(lpYit_pre[c]);
    }
  }
  return Score;
}

// [[Rcpp::export]]
IntegerVector PWKL_cpp(NumericMatrix X, NumericVector Pi, NumericMatrix theta, int st, int i){
  int C = theta.nrow();
  int p = theta.ncol();
  int j = NA_INTEGER;
  IntegerVector stp = IntegerVector::create(st,p);
  int st2 = min(stp);
  IntegerVector unasked = seq(1, p);
  IntegerVector asked(st2);
  int remaining = NA_INTEGER;
  NumericVector lpYit_pre = log(Pi);
  NumericVector lpYit_pre_ (C, NA_REAL);
  int jjj = NA_INTEGER;
  for(unsigned int ite = 0; ite < st2; ++ite) {
    NumericVector Score = PWKL_Score_cpp(theta, unasked, lpYit_pre);
    double maxScore = max(Score);
    jjj = 0;
    remaining = unasked.size();
    IntegerVector unasked_ (remaining-1, NA_INTEGER);
    for(unsigned int jj = 0; jj < remaining; ++jj) {
      unasked_[jjj] = unasked[jj];
      jjj = jjj + 1;
      if(Score[jj] == maxScore){
        j = unasked[jj];
        jjj = jjj - 1;
      }
    }
    asked[ite] = j;
    unasked = unasked_;
    if(NumericVector::is_na(X(i-1,j-1))) {
      lpYit_pre_ = lpYit_pre;
    }else if(X(i-1,j-1)==1){
      for(unsigned int c = 0; c < C; ++c) {
        lpYit_pre_[c] = lpYit_pre[c] + log( theta(c,j-1) );
      }
    }else if(X(i-1,j-1)==0){
      for(unsigned int c = 0; c < C; ++c) {
        lpYit_pre_[c] = lpYit_pre[c] + log( 1-theta(c,j-1) );
      }
    }
    lpYit_pre = lpYit_pre_ - log(sum(exp(lpYit_pre_)));
  }
  return asked;
}