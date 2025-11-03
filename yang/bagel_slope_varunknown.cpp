#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <algorithm>

struct UpdateOut {
  Eigen::VectorXd mu_post;
  Eigen::MatrixXd Sigma_post;
  double alpha_prior;
  double beta_prior;
  double w;
};

inline double student_t_pdf(double x, double df, double mu, double sigma) {
  // density of Student-t(df, mu, sigma)
  double z = (x - mu) / sigma;
  double c = std::tgamma((df + 1.0) / 2.0) /
             (std::sqrt(df * M_PI) * std::tgamma(df / 2.0) * sigma);
  return c * std::pow(1.0 + (z * z) / df, -(df + 1.0) / 2.0);
}
// Build seasonal design h(t):
//  - L == 2 + (S-1): h = [1, t, onehot(S-1)]
//  - L == 4 + (S-1): h = [1, t, 1, t, onehot(S-1)]
inline Eigen::VectorXd make_h_seasonal(std::size_t L, int t, int season_period) {
  const int S = std::max(1, season_period);
  const int K = S - 1;
  const int idx = (t - 1) % S; // 0..S-1

  if (S == 1) {
    if (L == 2) { Eigen::VectorXd h(2); h << 1.0, double(t); return h; }
    if (L == 4) { Eigen::VectorXd h(4); h << 1.0, double(t), 1.0, double(t); return h; }
    throw std::runtime_error("make_h_seasonal: bad L for S=1.");
  }

  if (L == std::size_t(2 + K)) {
    Eigen::VectorXd h = Eigen::VectorXd::Zero(2 + K);
    h(0) = 1.0; h(1) = double(t);
    if (idx > 0) h(2 + (idx - 1)) = 1.0; // category 1..S-1
    return h;
  }
  if (L == std::size_t(4 + K)) {
    Eigen::VectorXd h = Eigen::VectorXd::Zero(4 + K);
    h(0) = 1.0; h(1) = double(t);
    h(2) = 1.0; h(3) = double(t);
    if (idx > 0) h(4 + (idx - 1)) = 1.0;
    return h;
  }
  throw std::runtime_error("make_h_seasonal: L must be 2+K or 4+K.");
}

// Known-variance update (Gaussian predictive; covariance form)
// Q = sigma2 + h' Σ h, A = Σ h / Q, Σ' = Σ - A A' Q, μ' = μ + A e
inline UpdateOut update_post_parameter_unknown(
    const Eigen::VectorXd& mu_prior,
    const Eigen::MatrixXd& Lambda_prior,
    double alpha_prior,
    double beta_prior,
    double obs,
    double w_prev,
    int t,
    int season_period)
{
  const std::size_t L = std::size_t(mu_prior.size());
  if (Lambda_prior.rows() != int(L) || Lambda_prior.cols() != int(L))
    throw std::runtime_error("Lambda_prior shape mismatch.");

  // design vector
  const Eigen::VectorXd h = make_h_seasonal(L, t, season_period);

  // posterior updates
  const double alpha_post = alpha_prior + 0.5;
  const double A = 1.0 + (h.transpose() * Lambda_prior * h)(0, 0);
  const double B = obs - (h.transpose() * mu_prior)(0, 0);
  const Eigen::VectorXd C = (Lambda_prior * h) / A;

  Eigen::MatrixXd Lambda_post = Lambda_prior - (C * C.transpose()) * A;
  Eigen::VectorXd mu_post     = mu_prior + C * B;
  double beta_post            = beta_prior + 0.5 * (B * B) / A;

  // marginal likelihood (Student-t)
  const double df = 2.0 * alpha_prior;
  const double mean = (h.transpose() * mu_prior)(0, 0);
  const double scale = std::sqrt(beta_prior * A / alpha_prior);
  const double like = student_t_pdf(obs, df, mean, scale);
  const double w = w_prev * like;

  return {mu_post, Lambda_post, alpha_post, beta_post, w};
}
// 4×4 block prior update 
inline std::pair<Eigen::Vector4d, Eigen::Matrix4d> update_prior_4x4(const Eigen::Matrix4d& Sig0,
                 const Eigen::Vector4d& mu0,
                 const Eigen::Vector2d& mu_post,
                 const Eigen::Matrix2d& Sigma_post){
  const Eigen::Matrix2d Sig_bb = Sig0.topLeftCorner<2,2>();
  const Eigen::Matrix2d Sig_bg = Sig0.topRightCorner<2,2>();
  const Eigen::Matrix2d Sig_gb = Sig0.bottomLeftCorner<2,2>();
  const Eigen::Matrix2d Sig_gg = Sig0.bottomRightCorner<2,2>();

  const Eigen::Vector2d mu_beta  = mu0.head<2>();
  const Eigen::Vector2d mu_gamma = mu0.tail<2>();

  const Eigen::Matrix2d Sig_bb_inv = Sig_bb.inverse();

  const Eigen::Vector2d mu_gamma_updated = mu_gamma + Sig_gb * Sig_bb_inv * (mu_post - mu_beta);

  const Eigen::Matrix2d Sig_gg_updated = Sig_gg + Sig_gb * Sig_bb_inv * (Sigma_post - Sig_bb) * Sig_bb_inv * Sig_bg;
  const Eigen::Matrix2d Sig_bg_updated = Sig_gb * Sig_bb_inv * Sigma_post;
  const Eigen::Matrix2d Sig_gb_updated = Sigma_post * Sig_bb_inv * Sig_bg;

  Eigen::Vector4d mu_updated;
  mu_updated << mu_post, mu_gamma_updated;

  Eigen::Matrix4d Sig_updated;
  Sig_updated.topLeftCorner<2,2>()     = Sigma_post;
  Sig_updated.topRightCorner<2,2>()    = Sig_bg_updated;
  Sig_updated.bottomLeftCorner<2,2>()  = Sig_gb_updated;
  Sig_updated.bottomRightCorner<2,2>() = Sig_gg_updated;

  return {mu_updated, Sig_updated};
}

// Inject 4×4 update into (4+K) state and retain seasonal blocks from an existing prior
inline std::pair<Eigen::VectorXd, Eigen::MatrixXd> rearrange_prior_seasonal(const std::pair<Eigen::Vector4d, Eigen::Matrix4d>& updated4,
                         int season_period,
                         const Eigen::VectorXd& mu_nochange,
                         const Eigen::MatrixXd& Sigma_nochange){
  const int S = std::max(1, season_period);
  const int K = S - 1;
  const int L = 4 + K;

  if (mu_nochange.size() < 2 + K) throw std::runtime_error("mu_nochange too short.");
  if (Sigma_nochange.rows() != mu_nochange.size() ||
      Sigma_nochange.cols() != mu_nochange.size())
    throw std::runtime_error("Sigma_nochange shape mismatch.");

  Eigen::VectorXd mu0 = Eigen::VectorXd::Zero(L);
  Eigen::MatrixXd Sig0 = Eigen::MatrixXd::Zero(L, L);

  mu0.head<4>() = updated4.first;
  Sig0.topLeftCorner<4,4>() = updated4.second;

  if (K > 0) {
    // keep seasonal mean/cov blocks from the "nochange" prior
    mu0.segment(4, K) = mu_nochange.segment(2, K); // note: mu_nochange expected to be length (2+K) or (4+K); we take seasonal means relative to first 2 states
    // seasonal x beta(1:2)
    Sig0.block(4, 0, K, 2) = Sigma_nochange.block(2, 0, K, 2);
    // beta(1:2) x seasonal
    Sig0.block(0, 4, 2, K) = Sigma_nochange.block(0, 2, 2, K);
    // seasonal x seasonal
    Sig0.block(4, 4, K, K) = Sigma_nochange.block(2, 2, K, K);
  }
  return {mu0, Sig0};
}

