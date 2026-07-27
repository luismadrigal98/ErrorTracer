## Referee check 2: is env_var comparable between the iid and ar() paths?
## Ground truth: with beta ~ 2 and delta = 0.2, V_env should be ~ beta^2*delta^2 = 0.16
## at EVERY lead time, in both models.
suppressMessages({library(ErrorTracer); library(brms)})
set.seed(42); n <- 60; H <- 15
x <- rnorm(n); phi <- 0.9; sig <- 0.5
e <- numeric(n); e[1] <- rnorm(1,0,sig/sqrt(1-phi^2))
for (t in 2:n) e[t] <- phi*e[t-1] + rnorm(1,0,sig)
y <- 1 + 2*x + e
train <- data.frame(y=y, x=x, t=1:n)
newd  <- data.frame(y=NA_real_, x=rnorm(H), t=(n+1):(n+H))

run <- function(form, tag){
  f <- et_fit(form, data=train, chains=2, iter=1500, warmup=750, cores=2,
              refresh=0, seed=1, priors=NULL)
  p <- et_predict(f, newdata=newd, env_noise=list(x=0.2), n_draws=1500, n_perturb=500)
  d <- decompose_uncertainty(p)
  b <- mean(as.matrix(f$fit, variable="b_x")[,1])
  cat("\n########", tag, " posterior beta =", round(b,3),
      " => analytic V_env = beta^2*delta^2 =", round(b^2*0.04,4), "\n")
  # the two arrays that env_var subtracts
  v_pert <- apply(p$lp_perturbed, 2, var)
  v_lp   <- apply(p$posterior_linpred[1:500,,drop=FALSE], 2, var)
  print(round(data.frame(lead=1:H,
      var_lp_perturbed = v_pert,        # manual LP, perturbed x, NO ar error
      var_posterior_linpred = v_lp,     # brms linpred, unperturbed, WITH ar error
      raw_difference = v_pert - v_lp,
      reported_env_var = d$env_var,
      env_mcse = d$v_env_mcse), 4))
  invisible(d)
}
d_iid <- run(y ~ x,                        "IID  (y ~ x)")
d_ar  <- run(y ~ x + ar(time=t, p=1),      "AR(1) (y ~ x + ar())")
cat("\n=== summary: mean reported env_var ===\n")
cat(" iid  :", round(mean(d_iid$env_var),4), "\n")
cat(" ar(1):", round(mean(d_ar$env_var),4),
    "   (leads with env_var floored to exactly 0: ",
    sum(d_ar$env_var == 0), "/", nrow(d_ar), ")\n")
