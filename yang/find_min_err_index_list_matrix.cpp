#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#include <RcppEigen.h>
#include <string>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wignored-attributes"
// #include <Eigen/Dense>
// #pragma GCC diagnostic pop
using namespace Rcpp;

double err_diff(double w, double mu_1, double mu_2, double sigma_1, double sigma_2) {
//  double a = (mu_1 - mu_2) / (sigma_2 - sigma_1);
//  double b = sqrt(a * a + log(sigma_2 / sigma_1) / (sigma_2 - sigma_1));
  double a = 1/sigma_2-1/sigma_1;
  double b = 2*(mu_1/sigma_1-mu_2/sigma_2);
  double c = pow(mu_2,2)/sigma_2-pow(mu_1,2)/sigma_1+log(sigma_2/sigma_1);
  // double a = mu_1*sigma_2-mu_2*sigma_1;
  // double b = sqrt(sigma_1*sigma_2)*sqrt(pow((mu_1-mu_2),2)+(sigma_2-sigma_1)*log(sigma_2/sigma_1));
  // double c_1=(a-b)/(sigma_2-sigma_1);
  // double c_2=(a+b)/(sigma_2-sigma_1);
  double delta = sqrt(pow(b,2)-4*a*c);
  double c_1 = (-b-delta)/2/a;
  double c_2 = (-b+delta)/2/a;
  if (c_1 > c_2) {
  std::swap(c_1, c_2);  // Ensures c_1 is the smaller root and c_2 is the larger root
}
  std::vector<double> terms(4);
  terms[0] = R::pnorm(c_2,mu_1,sqrt(sigma_1), TRUE, FALSE);
  terms[1] = R::pnorm(c_1,mu_1,sqrt(sigma_1), TRUE, FALSE);
  terms[2] = R::pnorm(c_1,mu_2,sqrt(sigma_2), TRUE, FALSE);
  terms[3] = R::pnorm(c_2,mu_2,sqrt(sigma_2), TRUE, FALSE);
//  terms[0] = a * sqrt(sigma_1) + b * sqrt(sigma_2);
//  terms[1] = a * sqrt(sigma_1) - b * sqrt(sigma_2);
//  terms[2] = a * sqrt(sigma_2) - b * sqrt(sigma_1);
//  terms[3] = a * sqrt(sigma_2) + b * sqrt(sigma_1);
//  for (int i = 0; i < 4; ++i) {
//    terms[i] = R::pnorm(terms[i],0.0,1.0, TRUE, FALSE); // Using CDF instead of PDF
//  }
//  double result = w * (terms[0]-terms[1]+terms[2]-terms[3]);
  double result = w * (terms[0]-terms[1]+terms[2]-terms[3]);
  return result;
}
//[[Rcpp::export]]

List find_min_err_index(const std::vector<double>& w, Rcpp::List mu,Rcpp::List Sigma) {
  int N = w.size();
  Eigen::Matrix<double, 1, 2> A;
  // Pre-transform all mu and Sigma
  std::vector<double> mu2d(N);
  std::vector<double> Sigma2d(N);
    for (int i = 0; i < N; ++i) {
      A << 1, 1;
      Eigen::VectorXd mu_i = Rcpp::as<Eigen::VectorXd>(mu[i]);       // 2x1
      Eigen::MatrixXd Sigma_i = Rcpp::as<Eigen::MatrixXd>(Sigma[i]); // 2x2
      mu2d[i] = (A * mu_i)(0);
      Sigma2d[i] = (A * Sigma_i * A.transpose())(0,0);
    }
    // Find minimum weighted KL divergence
  std::vector<double> errors;
  double min_err = std::numeric_limits<double>::max(); //std::vector<double>
  int min_index = -1;
  double min_ratio = 0.0;
  for (int i = 0; i < w.size()-1; ++i) {
    double err =err_diff(w[i], mu2d[i], mu2d[i+1], Sigma2d[i], Sigma2d[i+1]);
    errors.push_back(err);
    double ratio = w[i]/w[i+1];
    if (err < min_err) {
      min_err = err;
      min_index = i;
      min_ratio = ratio;
    }}
  return List::create(Named("index") = min_index + 1, Named("ratio") = min_ratio);
//  return errors;
}
// List find_min_err_index(const std::vector<double>& w, const std::vector<double>& mu,
//                        const std::vector<double>& sigma) {
//     double min_err = std::numeric_limits<double>::max(); //std::vector<double>
//     int min_index = -1;
//     double min_ratio = 0.0;
// //    std::vector<double> errors;
//     for (int i = 0; i < w.size()-1; ++i) {
//       double err = err_diff(w[i], mu[i], mu[i+1], sigma[i], sigma[i+1]);
//       double ratio = w[i]/w[i+1];
//             if (err < min_err) {
//               min_err = err;
//               min_index = i;
//               min_ratio = ratio;
//             }
//             }
// //      errors.push_back(err);}
//       return List::create(Named("index") = min_index + 1, Named("ratio") = min_ratio);
// }
// 
