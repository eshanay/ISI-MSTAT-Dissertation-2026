#include <RcppArmadillo.h>

using namespace Rcpp;

// takes in as argument the dimension (N) of the Wigner matrix X_N to be simulated and
// the power (p) of the matrix X_N to be computed from which the principal submatrix of dimension (k)
// specified as a function of N is extracted for computing the sum of the entries from the upper triangular block
double wigner_submatrix_sum_cpp(int N, double p) {
  
  // compute the submatrix dimension as a function of N
  int k = (int) std::sqrt((double)N);
  
  // Safety checks to ensure k is valid as a dimension of a matrix
  if (k < 1) k = 1;
  if (k > N) k = N;
  
  // Generating a Gaussian Wigner Matrix (scaled GOE)
  // starting with a matrix of i.i.d. standard normal entries
  arma::mat X = arma::randn<arma::mat>(N, N);
  
  // symmetrizing and scaling
  arma::mat W = (X + X.t()) / std::sqrt(2.0 * N);
  
  // Computing Matrix Power using Eigen-decomposition
  // W is symmetric, so we use efficient symmetric decomposition
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, W);
  
  // Raising eigenvalues to power p
  arma::vec eigval_pow = arma::pow(eigval, p);
  
  // Reconstructing W^p = V * D^p * V^T
  arma::mat W_power = eigvec * arma::diagmat(eigval_pow) * eigvec.t();
  
  // Extracting Principal k x k Submatrix
  arma::mat submatrix = W_power.submat(0, 0, k - 1, k - 1);
  
  // Sum all entries in the submatrix
  double total_sum = arma::accu(submatrix);
  
  return total_sum;
}
