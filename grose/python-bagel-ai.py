"""
Simple implementation of Algorithm 1 (Bagel) for the univariate
change-in-mean model with known variance.

Assumptions / simplifications:
- Observations y_t ~ N(mu, obs_var)
- Prior for any segment mean ~ N(mu0, prior_var)
- Known observation variance obs_var
- Priors for change-location ratio terms are assumed equal -> ratio factor = 1
- Merge uses numeric approximation of TV distance on a grid (1D)
References: Bagel paper (Algorithm 1, merge rule, Theorems 1-4). 4
"""
import numpy as np
from math import sqrt, log, exp
from scipy.stats import norm

# --- Utility functions for normal-normal conjugacy (known obs variance) ---


def posterior_normal_known_variance(prior_mean, prior_var, sum_y, n, obs_var):
    """Posterior mean and variance for normal prior + n obs with known variance.
    Data model: y_i ~ N(theta, obs_var), prior theta ~ N(prior_mean, prior_var).
    sum_y: sum of observations in that segment.
    """
    if n == 0:
        return prior_mean, prior_var
    prior_prec = 1.0 / prior_var
    obs_prec = n / obs_var
    post_var = 1.0 / (prior_prec + obs_prec)
    post_mean = post_var * (prior_mean * prior_prec + (sum_y / obs_var))
    return post_mean, post_var


def predictive_logpdf_normal(post_mean, post_var, obs_var, y):
    """Log predictive density for next observation given posterior of mean:
       Predictive variance = obs_var + post_var
    """
    pred_var = obs_var + post_var
    return norm.logpdf(y, loc=post_mean, scale=sqrt(pred_var))


def tv_distance_gaussians_grid(mu1, var1, mu2, var2, ngrid=200, width=6.0):
    """
    Approximate total variation distance between two univariate Gaussians
    via numeric integration on a grid covering both means +/- width * sqrt(var).
    TV = 0.5 * integral |f1 - f2|.
    """
    s1 = sqrt(var1)
    s2 = sqrt(var2)
    left = min(mu1 - width * s1, mu2 - width * s2)
    right = max(mu1 + width * s1, mu2 + width * s2)
    xs = np.linspace(left, right, ngrid)
    f1 = norm.pdf(xs, loc=mu1, scale=s1)
    f2 = norm.pdf(xs, loc=mu2, scale=s2)
    tv = 0.5 * np.trapz(np.abs(f1 - f2), xs)
    return float(tv)


# --- Bagel simplified class for 1D change-in-mean ---


class Bagel1D:
    def __init__(self, mu0=0.0, prior_var=1.0, obs_var=1.0, M=50, c_detect=0.5):
        """
        mu0, prior_var: prior for any segment mean ~ N(mu0, prior_var)
        obs_var: known observation variance
        M: maximum number of distinct posterior components to keep (merge to M)
        c_detect: detection threshold (if 1 - w0,t >= c_detect => announce change)
        """
        self.mu0 = mu0
        self.prior_var = prior_var
        self.obs_var = obs_var
        self.M = M
        self.c_detect = c_detect

        # T holds boundaries: T = {0, i1, i2, ...} but we track only the group indices.
        # Initially only "no-change" component (index 0).
        self.N = 0  # number of distinct posterior components excluding the "no change" 0
        # weights for components w0, w1, ..., wN (un-normalized or normalized)
        self.w = [1.0]  # start with w0 = 1 (no-change prior mass)
        # ratio vectors r(i) for each run of change locations (not needed for 1D demo beyond bookkeeping)
        # we store r as list of arrays (when a component spans several change locations)
        self.r = []  # empty initially because no change candidates
        # For each component i=0..N we store posterior parameters for theta (mean, var)
        # For i=0 (no change): this represents the pre-change parameter posterior
        # For i>=1: posterior for post-change parameter (theta)
        # Each component will maintain (n, sum_y, post_mean, post_var)
        # Note: n and sum_y are counts/sums of observations assigned to that component's posterior.
        self.components = []
        # initialize component 0 (no-change) with zero observations
        self.components.append({'n': 0, 'sum': 0.0,
                                'post_mean': self.mu0, 'post_var': self.prior_var})

        # keep time index
        self.t = 0

    def _add_new_candidate_from_prev_no_change(self):
        """
        Introduce new changepoint candidate for tau = t-1 (i.e. change just before current obs)
        According to Theorem 2, the new candidate posterior can be derived from the no-change posterior.
        For the simple independent priors here, we initialize with prior (no data yet),
        but after seeing the current obs it will be updated immediately in the update step.
        We'll initialize (n=0, sum=0) but store a theta posterior equal to prior.
        """
        self.N += 1
        # component index is N (component 0 is no-change)
        self.w.append(0.0)  # placeholder weight (will be computed)
        # ratio vector for this new run: initially it's a single location, r = [1.0]
        self.r.append(np.array([1.0], dtype=float))
        # create component with zero data (posterior = prior)
        self.components.append({'n': 0, 'sum': 0.0, 'post_mean': self.mu0, 'post_var': self.prior_var})

    def _predictive_logpdf_for_component(self, comp_idx, y):
        """Compute log predictive p(y | component comp_idx) using current component posterior."""
        comp = self.components[comp_idx]
        mu_post = comp['post_mean']
        var_post = comp['post_var']
        return predictive_logpdf_normal(mu_post, var_post, self.obs_var, y)

    def update(self, y):
        """
        Process a single observation y (one iteration of Algorithm 1).
        Returns:
            detected (bool): whether we detected a change this update
            info (dict): current state info (weights, components summary)
        """
        self.t += 1
        # 1) Introducing new changepoint candidate tau = t-1
        # Create new candidate component derived from "no change" posterior (Theorem 2)
        # For simplicity initialize new candidate with prior/posterior equal to prior (will update below)
        self._add_new_candidate_from_prev_no_change()

        # 2) Updating weights using Theorem 1 and predictive densities
        # Calculate raw new weights wi,t = wi,t-1 * (prior-ratio assumed 1) * p(y | component i)
        new_w = []
        log_preds = []
        for i in range(len(self.w)):
            # compute predictive logpdf using current component posterior (which is based on data up to t-1)
            lp = self._predictive_logpdf_for_component(i, y)
            log_preds.append(lp)
            # weights multiplied by predictive (we work in normal space)
            new_w.append(self.w[i] * exp(lp))

        new_w = np.array(new_w, dtype=float)

        # 3) Update parameters of posterior distributions (Theorem 3)
        # For each component, update its sufficient stats deterministically with observed y.
        # For no-change component (i=0): it accumulates pre-change data
        # For post-change components (i>=1): they are for runs that start at different taus;
        # when an observation y arrives, every existing component's posterior must be updated
        # with y if that component corresponds to a model that includes y in its segment.
        # In this simplified representation: we update component i by incrementing its (n, sum)
        # except for the no-change component which also increments.
        # Note: in the full Bagel representation components correspond to *distinct* posteriors that represent multiple taus;
        # here each component corresponds to a single posterior and we treat them directly.
        for i, comp in enumerate(self.components):
            # deterministic update: increment n and sum
            comp['n'] += 1
            comp['sum'] += y
            # recompute posterior mean, var
            mu_post, var_post = posterior_normal_known_variance(
                prior_mean=self.mu0,
                prior_var=self.prior_var,
                sum_y=comp['sum'],
                n=comp['n'],
                obs_var=self.obs_var
            )
            comp['post_mean'] = mu_post
            comp['post_var'] = var_post

        # 4) Normalize new weights
        total = new_w.sum()
        if total <= 0:
            # numerical guard
            new_w = np.ones_like(new_w) / len(new_w)
        else:
            new_w = new_w / total

        # assign new weights
        self.w = list(new_w)

        # 5) Decide if there is a changepoint: if 1 - w0 >= c_detect
        detected = (1.0 - self.w[0]) >= self.c_detect

        # 6) Merging step if number of distinct posteriors N equals M
        # Note: N is number of post-change components; we allow at most self.M components (excluding w0).
        # Current number of distinct components equals len(self.w)-1 (excluding w0)
        while (len(self.w) - 1) > self.M:
            # Find i in {1..N} minimizing wi * TV(fi, fi+1)
            # Here fi is component i, fi+1 is component i+1
            best_i = None
            best_score = float('inf')
            for i in range(1, len(self.w) - 1):
                # components i and i+1 -> take their posterior normals
                mu_i = self.components[i]['post_mean']
                var_i = self.components[i]['post_var'] * 1.0  # post_var
                mu_j = self.components[i + 1]['post_mean']
                var_j = self.components[i + 1]['post_var'] * 1.0
                tv = tv_distance_gaussians_grid(mu_i, var_i, mu_j, var_j, ngrid=300)
                score = self.w[i] * tv
                if score < best_score:
                    best_score = score
                    best_i = i
            if best_i is None:
                break

            i = best_i
            # Combine component i and i+1: new weight, new ratio vector, and merge posterior deterministically
            w_i = self.w[i]
            w_j = self.w[i + 1]
            # new weight at position i+1 after removal will be w_i + w_j and stored at position i+1 (paper combines into j)
            new_weight = w_i + w_j

            # Merge posterior parameters: we form weighted mixture and then "recover" a single Gaussian posterior
            # Simple approach: merge by moment-matching (weight-average means and combine variances)
            mu_i = self.components[i]['post_mean']
            var_i = self.components[i]['post_var']
            mu_j = self.components[i + 1]['post_mean']
            var_j = self.components[i + 1]['post_var']
            # mixture mean
            mix_mean = (w_i * mu_i + w_j * mu_j) / (w_i + w_j)
            # mixture variance = weighted (var + mean^2) - mix_mean^2
            mix_second = (w_i * (var_i + mu_i * mu_i) + w_j * (var_j + mu_j * mu_j)) / (w_i + w_j)
            mix_var = max(1e-12, mix_second - mix_mean * mix_mean)

            # Update: we will remove component i and store the merged posterior in slot (i+1) per paper's prescription
            # Here for code simplicity replace component i+1 with merged and delete i
            self.components[i + 1]['post_mean'] = mix_mean
            self.components[i + 1]['post_var'] = mix_var
            # We cannot reconstruct exact n and sum (they differ per run), so set counts to 1 to keep update working.
            # (This is a simplification: the full algorithm keeps ratio vectors to reconstruct posterior per tau.)
            self.components[i + 1]['n'] = 1
            self.components[i + 1]['sum'] = mix_mean  # approximate

            # update weight and remove i
            self.w[i + 1] = new_weight
            del self.w[i]
            # remove r(i) and components[i]
            if i - 1 >= 0 and (i - 1) < len(self.r):
                del self.r[i - 1]  # ratio vectors list indexes differ by 1 relative to components
            del self.components[i]

        # ensure weights normalized after merging
        w_arr = np.array(self.w)
        if w_arr.sum() <= 0:
            w_arr = np.ones_like(w_arr) / len(w_arr)
        else:
            w_arr = w_arr / w_arr.sum()
        self.w = list(w_arr)

        info = {
            't': self.t,
            'w': np.array(self.w),
            'components': [{'post_mean': c['post_mean'], 'post_var': c['post_var'], 'n': c['n']} for c in self.components],
            'N_distinct_posteriors': len(self.w) - 1
        }
        return detected, info


# --- Quick demo on synthetic data ---
if __name__ == "__main__":
    np.random.seed(2)
    # generate data: mean 0 for first 80 pts then mean 3 afterwards
    obs_var = 1.0
    data = np.concatenate([np.random.normal(0, sqrt(obs_var), size=80),
                           np.random.normal(3.0, sqrt(obs_var), size=120)])

    bagel = Bagel1D(mu0=0.0, prior_var=1.0, obs_var=obs_var, M=20, c_detect=0.5)
    detections = []
    for y in data:
        detected, info = bagel.update(y)
        detections.append(detected)

    # print the first time a detection occurred
    if any(detections):
        t_detect = np.where(detections)[0][0] + 1
        print("First detection at t =", t_detect)
    else:
        print("No detection in sequence.")

    # print final components summary
    print("Final weights:", info['w'])
    for i, comp in enumerate(info['components']):
        print(f"Comp {i}: mean={comp['post_mean']:.3f}, var={comp['post_var']:.3f}, n={comp['n']}")
