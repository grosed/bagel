#include <iostream>
#include <cmath>
#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <string>

// [[Rcpp::depends(RcppEigen)]]
#include <Eigen/Dense>

using namespace Rcpp;

// perform matrix operations
// double KL_approx(const Eigen::MatrixXd& Sigma1, const Eigen::MatrixXd& Sigma2,
//                 const Eigen::VectorXd& mu1, const Eigen::VectorXd& mu2) {
//     Eigen::MatrixXd Sigma1_inv = Sigma1.inverse();
//     double trace_part = (Sigma1_inv * Sigma2 - Eigen::MatrixXd::Identity(Sigma1.rows(), Sigma1.cols())).trace();
//     Eigen::VectorXd mu_diff = mu1 - mu2;
//     double quadratic_form = (mu_diff.transpose() * Sigma1_inv * mu_diff)(0, 0);
//     double det_part = (Sigma2 * Sigma1_inv).determinant();
// 
//     double result = trace_part + quadratic_form - std:: log(det_part);
//     return result;
// }

double KL_approx(const Eigen::MatrixXd& Sigma1, const Eigen::MatrixXd& Sigma2,
                 const Eigen::VectorXd& mu1, const Eigen::VectorXd& mu2) {
  Eigen::MatrixXd Sigma2_inv = Sigma2.inverse();
  double trace_part = (Sigma2_inv * Sigma1 - Eigen::MatrixXd::Identity(Sigma1.rows(), Sigma1.cols())).trace();
  Eigen::VectorXd mu_diff = mu2 - mu1;
  double quadratic_form = (mu_diff.transpose() * Sigma2_inv * mu_diff)(0, 0);
  double det_part = (Sigma2_inv*Sigma1).determinant();
  
  double result = trace_part + quadratic_form - std:: log(det_part);
  return result;
}

//[[Rcpp::export]]
int  KL_min_index(const std::vector<double>& w, Rcpp::List Sigma,Rcpp::List mu,char model_type) {
  //  std::vector<double> errors;
  
  int N = w.size();
  Eigen::Matrix<double, 2, 4> A;
  
  // Pre-transform all mu and Sigma
  std::vector<Eigen::Vector2d> mu2d(N);
  std::vector<Eigen::Matrix2d> Sigma2d(N);
  if(model_type == 'd'){
    for (int i = 0; i < N; ++i) {
      A << 1, i+1, 1, i+1,
           0, 1, 0, 1;
      Eigen::VectorXd mu_i = Rcpp::as<Eigen::VectorXd>(mu[i]);       // 4x1
      Eigen::MatrixXd Sigma_i = Rcpp::as<Eigen::MatrixXd>(Sigma[i]); // 4x4
      mu2d[i] = A * mu_i;
      Sigma2d[i] = A * Sigma_i * A.transpose();
    }}else if(model_type == 'c'){
      for (int i = 0; i < N; ++i) {
        A << 0, 0, 1, i+1,
             0, 1, 0, 1;
        Eigen::VectorXd mu_i = Rcpp::as<Eigen::VectorXd>(mu[i]);       // 4x1
        Eigen::MatrixXd Sigma_i = Rcpp::as<Eigen::MatrixXd>(Sigma[i]); // 4x4
        mu2d[i] = A * mu_i;
        Sigma2d[i] = A * Sigma_i * A.transpose();
      }}
  

  
  // Find minimum weighted KL divergence
  double min_err = std::numeric_limits<double>::max(); //std::vector<double> 
  int min_index = -1;
  for (int i = 0; i < w.size()-1; ++i) {
    double err = w[i]*KL_approx(Sigma2d[i],Sigma2d[i+1],
                           mu2d[i],mu2d[i+1]);
     if (err < min_err) {
        min_err = err;
        min_index = i;
    }}
    
  return min_index + 1;
}


