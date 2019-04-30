#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;
// This is a simple function using Rcpp that creates an R list
// containing a character vector and a numeric vector.
//
// Learn more about how to use Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   http://gallery.rcpp.org/
//

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
int findMaxMean(arma::vec num, int binNum, int lengthNum){
  // num is bins, binNum is bin's number,
  // arma::vec values = arma::zeros<arma::vec>(binNum);

  float sum;
  int i = 0;
  int big;
  int max;

  while ( i <= (lengthNum-binNum))
  {
    sum = 0;
    for (int j = i; j < (i+binNum); ++j)
    {
      sum = sum + num[j];
    }
    int mean = sum/binNum;

    if ( mean > max )
    {
      max = mean;
      big = i;
    }

    i++;
  }

  // for (int k = big; k < (big + binNum); ++k)
  // {
  //   values(k-big) = num[k];
  // }
  return (big+1);
}
