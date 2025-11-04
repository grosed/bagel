

update_prior <- function(Sig0,mu0,mu_post,Sigma_post)
{
  Sig_beta_beta <- Sig0[1:1, 1:1]
  Sig_beta_gamma <- Sig0[1:1, 2:2]
  Sig_gamma_beta <- Sig0[2:2, 1:1]
  Sig_gamma_gamma <- Sig0[2:2, 2:2]
  
  mu_beta <- mu0[1:1]
  mu_gamma <- mu0[2:2]
  
  # Update mean
  Sig_beta_beta_inv <- solve(Sig_beta_beta)
  mu_gamma_updated <- mu_gamma + Sig_gamma_beta %*% Sig_beta_beta_inv %*% (mu_post - mu_beta)  # (mu_beta - mu_beta) = 0
  mu_updated <- c(mu_post, mu_gamma_updated)  # so mu_updated = mu0_dis
  
  # Update covariance
  Sig_gamma_gamma_updated <- Sig_gamma_gamma + Sig_gamma_beta %*% Sig_beta_beta_inv %*% 
    (Sigma_post - Sig_beta_beta) %*% Sig_beta_beta_inv %*% Sig_beta_gamma
  Sig_beta_gamma_updated <- Sig_gamma_beta%*%Sig_beta_beta_inv %*%Sigma_post
  Sig_gamma_beta_updated <- Sigma_post %*% Sig_beta_beta_inv %*% Sig_beta_gamma
  Sig_updated <- rbind(
    cbind(Sigma_post, Sig_beta_gamma_updated),
    cbind(Sig_gamma_beta_updated, Sig_gamma_gamma_updated)
  )
  return(list(mu=mu_updated,Sig=Sig_updated))
}