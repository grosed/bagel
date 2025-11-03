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
  double w;
};

inline double normal_pdf(double x, double mean, double var) {
  const double two_pi = 2.0 * M_PI;
  return std::exp(-0.5 * (x-mean)*(x-mean) / var) / std::sqrt(two_pi * var);
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
inline UpdateOut update_post_parameter_known(const Eigen::VectorXd& mu_prior,
                                             const Eigen::MatrixXd& Sigma_prior,
                                             double sigma2,
                                             double obs,
                                             double w_prev,
                                             int t,
                                             int season_period)
{
  const std::size_t L = std::size_t(mu_prior.size());
  if (Sigma_prior.rows() != int(L) || Sigma_prior.cols() != int(L))
    throw std::runtime_error("Sigma_prior shape mismatch.");

  const Eigen::VectorXd h = make_h_seasonal(L, t, season_period);
  const double Q = sigma2 + (h.transpose() * Sigma_prior * h)(0,0);
  const double e = obs - (h.transpose() * mu_prior)(0,0);

  const Eigen::VectorXd A = (Sigma_prior * h) / Q;
  const Eigen::MatrixXd Sigma_post = Sigma_prior - (A * A.transpose()) * Q;
  const Eigen::VectorXd mu_post    = mu_prior + A * e;

  const double var_pred  = (h.transpose() * Sigma_prior * h)(0,0) + sigma2;
  const double mean_pred = (h.transpose() * mu_prior)(0,0);
  const double like      = normal_pdf(obs, mean_pred, var_pred);

  return {mu_post, Sigma_post, w_prev * like};
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

// -----------------------------------------
// OBCD loop: known-variance + (optional) seasonality
//  - season_period = 1  -> no seasonality
//  - season_period > 1  -> categorical season (S-1 one-hot)
// -----------------------------------------
struct OBCDResultKV {
  std::vector<double> uncertainty_dis;
  std::vector<double> uncertainty_conti;
  int stoppingtime; // t or +inf
  std::vector<double> change_trace; // size TT (1..TT), NA as NaN
  std::vector<int>    prune_dis, prune_conti;         // pruned candidate indices (time)
  std::vector<int>    candidate_dis, candidate_conti; // all candidate times kept in each stream
  std::vector<double> error_ratio_dis, error_ratio_conti;
  // debugging/inspection
  std::vector<Eigen::VectorXd> mu_dis, mu_conti;
  std::vector<Eigen::MatrixXd> Sigma_dis, Sigma_conti;
};

inline OBCDResultKV OBCD_slope_varknown_seasonal(const std::vector<double>& x,
                                                 double p_dis,
                                                 double p_conti,
                                                 double decay,
                                                 const Eigen::Vector2d& mu0_base,  // length 2
                                                 double delta_1, double delta_2,    // base trend variances
                                                 double delta_3, double delta_4,    // aux for 4x4 priors (if prior=="prior2")
                                                 double sigma2,
                                                 double thresh,
                                                 const std::string& error = "benchmark",
                                                 int Kmax = 100,
                                                 const std::string& prior = "prior1",
                                                 int season_period = 1,
                                                 const Eigen::VectorXd& season_mu = Eigen::VectorXd(),
                                                 const Eigen::VectorXd& season_Sig_diag = Eigen::VectorXd()){
  const int TT = int(x.size());
  if (TT <= 0) throw std::runtime_error("x must be non-empty.");

  const int S = std::max(1, season_period);
  const int K = S - 1;

  // Build initial (2+K)-dim prior: mu0 = [μ0_base, season_mu], Σ0=diag(delta_1, delta_2, season_Sig_diag)
  Eigen::VectorXd mu0 = (K>0) ? (Eigen::VectorXd(2+K)) : (Eigen::VectorXd(2));
  mu0.head<2>() = mu0_base;
  if (K > 0) {
    if (season_mu.size() != K) throw std::runtime_error("season_mu length must be S-1.");
    mu0.tail(K) = season_mu;
  }
  Eigen::MatrixXd Sigma0 = Eigen::MatrixXd::Zero(mu0.size(), mu0.size());
  Sigma0(0,0) = delta_1;
  Sigma0(1,1) = delta_2;
  if (K > 0) {
    if (season_Sig_diag.size() != K) throw std::runtime_error("season_Sig_diag length must be S-1.");
    for (int i=0;i<K;++i) Sigma0(2+i, 2+i) = season_Sig_diag(i);
  }

  const double p = p_dis + p_conti;

  // t=1 update
  UpdateOut s1 = update_post_parameter_known(mu0, Sigma0, sigma2, x[0], 1.0 - p, /*t=*/1, S);

  std::vector<Eigen::VectorXd> mu_dis_vec{ s1.mu_post };
  std::vector<Eigen::MatrixXd> Sig_dis_vec{ s1.Sigma_post };
  std::vector<double> w_dis_vec{ s1.w };

  std::vector<Eigen::VectorXd> mu_cont_vec{ s1.mu_post };
  std::vector<Eigen::MatrixXd> Sig_cont_vec{ s1.Sigma_post };
  std::vector<double> w_cont_vec{ s1.w };

  std::vector<int> candidate_dis{1}, candidate_cont{1};
  std::vector<int> prune_dis, prune_cont;
  std::vector<double> err_ratio_dis, err_ratio_cont;

  std::vector<double> be_change(TT, std::numeric_limits<double>::quiet_NaN());

  // decay weights trackers
  std::vector<double> g_dis_prev{1.0,1.0}, g_cont_prev{1.0,1.0};
  std::vector<double> g_dis{1.0,1.0}, g_cont{1.0,1.0};

  auto make_decay = [&](int t)->std::vector<double>{
    // produce q for indices (t-1):1 normalized (length t-1)
    std::vector<double> q;
    q.reserve(t-1);
    double sumq = 0.0;
    for (int k=t-1;k>=1;--k) { double v = std::pow(decay, double(k)); q.push_back(v); sumq += v; }
    for (double& v : q) v /= sumq;
    // prepend 1 for the "root"
    std::vector<double> g; g.reserve(1 + q.size()); g.push_back(1.0);
    g.insert(g.end(), q.begin(), q.end());
    return g;
  };

  for (int t=2; t<=TT; ++t) {
    const double y = x[t-1];

    // ---- current updating: DISCRETE
    const int len_dis = int(w_dis_vec.size());
    g_dis_prev = g_dis;
    g_dis = (t==2) ? std::vector<double>{1.0,1.0} : make_decay(t);

    std::vector<Eigen::VectorXd> mu_dis_new; mu_dis_new.reserve(len_dis);
    std::vector<Eigen::MatrixXd> Sig_dis_new; Sig_dis_new.reserve(len_dis);
    std::vector<double>           w_dis_new;  w_dis_new.reserve(len_dis);

    for (int n=0;n<len_dis;++n) {
      const double w_prev = w_dis_vec[n] * (g_dis[n] / g_dis_prev[n]);
      UpdateOut u = update_post_parameter_known(mu_dis_vec[n], Sig_dis_vec[n], sigma2, y, w_prev, t, S);
      mu_dis_new.push_back(u.mu_post);
      Sig_dis_new.push_back(u.Sigma_post);
      w_dis_new.push_back(u.w);
    }
    // append and then drop olds
    mu_dis_vec.insert(mu_dis_vec.end(), mu_dis_new.begin(), mu_dis_new.end());
    Sig_dis_vec.insert(Sig_dis_vec.end(), Sig_dis_new.begin(), Sig_dis_new.end());
    w_dis_vec.insert(w_dis_vec.end(), w_dis_new.begin(), w_dis_new.end());
    mu_dis_vec.erase(mu_dis_vec.begin(), mu_dis_vec.begin() + len_dis);
    Sig_dis_vec.erase(Sig_dis_vec.begin(), Sig_dis_vec.begin() + len_dis);
    w_dis_vec.erase(w_dis_vec.begin(), w_dis_vec.begin() + len_dis);

    // ---- current updating: CONTINUOUS
    const int len_cont = int(w_cont_vec.size());
    g_cont_prev = g_cont;
    g_cont = (t==2) ? std::vector<double>{1.0,1.0} : make_decay(t);

    std::vector<Eigen::VectorXd> mu_cont_new; mu_cont_new.reserve(len_cont);
    std::vector<Eigen::MatrixXd> Sig_cont_new; Sig_cont_new.reserve(len_cont);
    std::vector<double>           w_cont_new;  w_cont_new.reserve(len_cont);

    for (int n=0;n<len_cont;++n) {
      const double w_prev = w_cont_vec[n] * (g_cont[n] / g_cont_prev[n]);
      UpdateOut u = update_post_parameter_known(mu_cont_vec[n], Sig_cont_vec[n], sigma2, y, w_prev, t, S);
      mu_cont_new.push_back(u.mu_post);
      Sig_cont_new.push_back(u.Sigma_post);
      w_cont_new.push_back(u.w);
    }
    mu_cont_vec.insert(mu_cont_vec.end(), mu_cont_new.begin(), mu_cont_new.end());
    Sig_cont_vec.insert(Sig_cont_vec.end(), Sig_cont_new.begin(), Sig_cont_new.end());
    w_cont_vec.insert(w_cont_vec.end(), w_cont_new.begin(), w_cont_new.end());
    mu_cont_vec.erase(mu_cont_vec.begin(), mu_cont_vec.begin() + len_cont);
    Sig_cont_vec.erase(Sig_cont_vec.begin(), Sig_cont_vec.begin() + len_cont);
    w_cont_vec.erase(w_cont_vec.begin(), w_cont_vec.begin() + len_cont);

    // ---- new leaf: DISCRETE
    // prior1:
    {
      Eigen::Vector4d mu0_dis4;
      if (prior == "prior1") {
        mu0_dis4 << mu0_base(0), mu0_base(1), -mu0_base(1) * double(t), 0.0;
      } else { // prior2
        mu0_dis4 << mu0_base(0), mu0_base(1), 0.0, 0.0;
      }
      Eigen::Matrix4d Sig0_dis4 = Eigen::Matrix4d::Zero();
      if (prior == "prior1") {
        // matrix(data=c( d1,0,-d1,0, 0,d2,0,-d2, -d1,0,2*d1+t^2*d2,-t*d2, 0,-d2,-t*d2,2*d2), byrow=TRUE)
        Sig0_dis4 <<  delta_1, 0.0,     -delta_1,          0.0,
                      0.0,     delta_2,  0.0,             -delta_2,
                     -delta_1, 0.0,     (2*delta_1 + t*1.0*t*delta_2),  -t*delta_2,
                      0.0,    -delta_2, -t*delta_2,        2*delta_2;
      } else {
        Sig0_dis4 <<  delta_1, 0, 0, 0,
                      0, delta_2, 0, 0,
                      0, 0, delta_3, 0,
                      0, 0, 0, delta_4;
      }

      // prior update using the *current root candidate* (index 0)
      const Eigen::Vector2d mu_post2 = mu_dis_vec[0].head<2>();
      const Eigen::Matrix2d Sig_post2 = Sig_dis_vec[0].topLeftCorner<2,2>();
      auto upd = update_prior_4x4(Sig0_dis4, mu0_dis4, mu_post2, Sig_post2);

      // expand to include season (if any)
      // Use "nochange" mu/Sigma as the current root candidate prior (length 2+K)
      const auto [mu0_dis_full, Sig0_dis_full] =
        rearrange_prior_seasonal(upd, S, /*mu_nochange=*/mu_dis_vec[0].head(2 + K), /*Sigma_nochange=*/Sig_dis_vec[0].topLeftCorner(2 + K, 2 + K));

      const double w_prev_newleaf = w_dis_vec[0] * p_dis * ( g_dis.back() / g_dis_prev[0] );
      UpdateOut u = update_post_parameter_known(mu0_dis_full, Sig0_dis_full, sigma2, y, w_prev_newleaf, t, S);
      mu_dis_vec.push_back(u.mu_post);
      Sig_dis_vec.push_back(u.Sigma_post);
      w_dis_vec.push_back(u.w);
    }

    // ---- new leaf: CONTINUOUS
    {
      Eigen::Vector4d mu0_cont4;
      mu0_cont4 << mu0_base(0), mu0_base(1), 0.0, 0.0;
      Eigen::Matrix4d Sig0_cont4 = Eigen::Matrix4d::Zero();
      // matrix(data=c(d1,0,0,0, 0,d2,0,0, 0,0,(t-1)^2*d4, -(t-1)*d4, 0,0,-(t-1)*d4, d4), byrow=TRUE)
      Sig0_cont4 <<  delta_1, 0, 0, 0,
                     0, delta_2, 0, 0,
                     0, 0, (t-1.0)*(t-1.0)*delta_4, -(t-1.0)*delta_4,
                     0, 0, -(t-1.0)*delta_4,         delta_4;

      const Eigen::Vector2d mu_post2 = mu_cont_vec[0].head<2>();
      const Eigen::Matrix2d Sig_post2 = Sig_cont_vec[0].topLeftCorner<2,2>();
      auto upd = update_prior_4x4(Sig0_cont4, mu0_cont4, mu_post2, Sig_post2);

      const auto [mu0_cont_full, Sig0_cont_full] =
        rearrange_prior_seasonal(upd, S, /*mu_nochange=*/mu_cont_vec[0].head(2 + K), /*Sigma_nochange=*/Sig_cont_vec[0].topLeftCorner(2 + K, 2 + K));

      const double w_prev_newleaf = w_cont_vec[0] * p_conti * ( g_cont.back() / g_cont_prev[0] );
      UpdateOut u = update_post_parameter_known(mu0_cont_full, Sig0_cont_full, sigma2, y, w_prev_newleaf, t, S);
      mu_cont_vec.push_back(u.mu_post);
      Sig_cont_vec.push_back(u.Sigma_post);
      w_cont_vec.push_back(u.w);
    }

    // ---- marginalize
    double sum_dis = 0.0; for (double v : w_dis_vec) sum_dis += v;
    double sum_cont_tail = 0.0; for (std::size_t i=1;i<w_cont_vec.size();++i) sum_cont_tail += w_cont_vec[i];
    const double fac = sum_dis + sum_cont_tail;
    if (fac > 0.0) {
      for (double& v : w_dis_vec)  v /= fac;
      for (double& v : w_cont_vec) v /= fac; // note: dividing all is ok; be_change uses w_dis[0]
    }
    // change probability = 1 - w_dis[1] (R's first entry)
    be_change[t-1] = 1.0 - w_dis_vec[0];

    // ---- threshold stopping
    if (be_change[t-1] > thresh) {
      OBCDResultKV out;
      out.uncertainty_dis  = w_dis_vec;
      out.uncertainty_conti= w_cont_vec;
      out.stoppingtime     = t;
      out.change_trace     = be_change;
      out.prune_dis        = prune_dis;
      out.prune_conti      = prune_cont;
      out.candidate_dis    = [&]{auto v=candidate_dis; v.push_back(t); return v;}();
      out.candidate_conti  = [&]{auto v=candidate_cont; v.push_back(t); return v;}();
      out.error_ratio_dis  = err_ratio_dis;
      out.error_ratio_conti= err_ratio_cont;
      out.mu_dis = mu_dis_vec;
      out.mu_conti = mu_cont_vec;
      out.Sigma_dis = Sig_dis_vec;
      out.Sigma_conti = Sig_cont_vec;
      return out;
    }

    // ---- pruning (keep KL mode)
    if (int(w_dis_vec.size() + w_cont_vec.size()) >= Kmax) {
        if (error == "KL") {
            // ---- discrete model pruning ----
            int ind_dis = KL_min_index(
                std::vector<double>(w_dis_vec.begin() + 1, w_dis_vec.end()),
                std::vector<Eigen::MatrixXd>(Sig_dis_vec.begin() + 1, Sig_dis_vec.end()),
                std::vector<Eigen::VectorXd>(mu_dis_vec.begin() + 1, mu_dis_vec.end()),
                "d");
            ind_dis += 1; // because we ignored the first one
            double ratio_dis = w_dis_vec[ind_dis] /
                            std::max(1e-300, w_dis_vec[ind_dis + 1]);
            err_ratio_dis.push_back(ratio_dis);
            prune_dis.push_back(candidate_dis[ind_dis]);

            w_dis_vec[ind_dis + 1] += w_dis_vec[ind_dis]; // combine two. let the w_j+1 = w_j+w_j+1
            w_dis_vec.erase(w_dis_vec.begin() + ind_dis);
            mu_dis_vec.erase(mu_dis_vec.begin() + ind_dis);
            Sig_dis_vec.erase(Sig_dis_vec.begin() + ind_dis);
            candidate_dis.erase(candidate_dis.begin() + ind_dis);

            // ---- continuous model pruning ----
            int ind_cont = KL_min_index(
                std::vector<double>(w_cont_vec.begin() + 1, w_cont_vec.end()),
                std::vector<Eigen::MatrixXd>(Sig_cont_vec.begin() + 1, Sig_cont_vec.end()),
                std::vector<Eigen::VectorXd>(mu_cont_vec.begin() + 1, mu_cont_vec.end()),
                "c");
            ind_cont += 1; // because we ignored the first one
            double ratio_cont = w_cont_vec[ind_cont] /
                                std::max(1e-300, w_cont_vec[ind_cont + 1]);
            err_ratio_cont.push_back(ratio_cont);
            prune_cont.push_back(candidate_cont[ind_cont]);

            w_cont_vec[ind_cont + 1] += w_cont_vec[ind_cont];  // combine two. let the w_j+1 = w_j+w_j+1
            w_cont_vec.erase(w_cont_vec.begin() + ind_cont);
            mu_cont_vec.erase(mu_cont_vec.begin() + ind_cont);
            Sig_cont_vec.erase(Sig_cont_vec.begin() + ind_cont);
            candidate_cont.erase(candidate_cont.begin() + ind_cont);

        } 

        } else {
        prune_dis.push_back(std::numeric_limits<int>::quiet_NaN());
        prune_cont.push_back(std::numeric_limits<int>::quiet_NaN());
        err_ratio_dis.push_back(std::numeric_limits<double>::quiet_NaN());
        err_ratio_cont.push_back(std::numeric_limits<double>::quiet_NaN());
        }
    // append new candidate time
    candidate_dis.push_back(t);
    candidate_cont.push_back(t);
  }

  // finished without threshold crossing
  OBCDResultKV out;
  out.uncertainty_dis  = w_dis_vec;
  out.uncertainty_conti= w_cont_vec;
  out.stoppingtime     = std::numeric_limits<int>::infinity(); // emulate Inf
  out.change_trace     = be_change;
  out.prune_dis        = prune_dis;
  out.prune_conti      = prune_cont;
  out.candidate_dis    = candidate_dis;
  out.candidate_conti  = candidate_cont;
  out.error_ratio_dis  = err_ratio_dis;
  out.error_ratio_conti= err_ratio_cont;
  out.mu_dis = mu_dis_vec;
  out.mu_conti = mu_cont_vec;
  out.Sigma_dis = Sig_dis_vec;
  out.Sigma_conti = Sig_cont_vec;
  return out;
}




// -----------------------------------------
// After outputing the pruned results, recover the whole posterior
// pruned_candidates: the sequence of class labels that were pruned (1-based labels)
// stored_candidate:  the number of classes kept at the time you stored w
// w:                 current (pruned) weights vector (will be expanded)
// ratios:            the recorded w_pruned / w_next ratios, same order as pruned_candidates
// ---------- split_number ----------
inline std::pair<double,double> split_number(double weight, double ratio) {
  // x1 = ratio/(1+ratio)*weight;  x2 = 1/(1+ratio)*weight
  double denom = 1.0 + ratio;
  if (denom <= 0.0) return {0.0, 0.0};
  return { (ratio/denom) * weight, (1.0/denom) * weight };
}

inline std::vector<double>
recover_posterior(const std::vector<int>& pruned_candidates,
                  int stored_candidate,
                  std::vector<double> w,
                  const std::vector<double>& ratios){
  // result[0] = { {1}, {2}, ..., {stored_candidate} }
  // Then for each pruned_candidates[i], merge that class with its next neighbor.
  using VecI  = std::vector<int>;
  using Proc  = std::vector<VecI>;
  std::vector<Proc> result;
  result.reserve(pruned_candidates.size() + 1);

  Proc base;
  base.reserve(stored_candidate);
  for (int i = 1; i <= stored_candidate; ++i) base.push_back(VecI{ i });
  result.push_back(base);

  for (size_t i = 0; i < pruned_candidates.size(); ++i) {
    Proc cur = result.back();
    int target = pruned_candidates[i]; // 1-based label
    int pruned_index = -1;
    for (int j = 0; j < (int)cur.size(); ++j) {
      // find which subset contains 'target'
      if (std::find(cur[j].begin(), cur[j].end(), target) != cur[j].end()) {
        pruned_index = j;
        break;
      }
    }
    if (pruned_index < 0 || pruned_index+1 >= (int)cur.size()) {
      throw std::runtime_error("recover_posterior: invalid pruned_candidates / state.");
    }
    // combine subset j with subset j+1 (append elements)
    VecI combined = cur[pruned_index];
    combined.insert(combined.end(), cur[pruned_index+1].begin(), cur[pruned_index+1].end());
    cur[pruned_index+1] = std::move(combined);
    cur.erase(cur.begin() + pruned_index);
    result.push_back(std::move(cur));
  }
  // for i in 1:length(ratios):
  //   procedure = result[[length(result)-i+1]]
  //   class = index of subset that contains last pruned label
  //   w <- append(w, split_number(w[class], ratios[last-i+1]), after=class)
  //   w <- w[-class]
  const int L = (int)result.size(); // = pruned_candidates.size()+1
  if ((int)ratios.size() != L-1)
    throw std::runtime_error("recover_posterior: ratios length mismatch.");

  for (int i = 0; i < (int)ratios.size(); ++i) {
    const Proc& procedure = result[L - 1 - i];  // zero-based
    int target = pruned_candidates[(int)pruned_candidates.size() - 1 - i];

    // find class index (0-based) whose subset contains target
    int cls = -1;
    for (int j = 0; j < (int)procedure.size(); ++j) {
      if (std::find(procedure[j].begin(), procedure[j].end(), target) != procedure[j].end()) {
        cls = j; break;
      }
    }
    if (cls < 0 || cls >= (int)w.size())
      throw std::runtime_error("recover_posterior: class index out of range.");

    auto parts = split_number(w[cls], ratios[(int)ratios.size() - 1 - i]);
    w.insert(w.begin() + cls + 1, parts.first);
    w.insert(w.begin() + cls + 2, parts.second);
    w.erase(w.begin() + cls);
  }
  return w;
}