## Referee check 1: does the AR path really accumulate, and what does
## winsorization do to a near-unit-root predictive variance?
suppressMessages({library(ErrorTracer); library(brms)})
set.seed(42)
n <- 60; H <- 15
x <- rnorm(n)
phi <- 0.9; sig <- 0.5
e <- numeric(n); e[1] <- rnorm(1, 0, sig/sqrt(1-phi^2))
for (t in 2:n) e[t] <- phi*e[t-1] + rnorm(1,0,sig)
y <- 1 + 2*x + e
train <- data.frame(y=y, x=x, t=1:n)
newd  <- data.frame(y=NA_real_, x=rnorm(H), t=(n+1):(n+H))

fit <- et_fit(y ~ x + ar(time = t, p = 1), data = train,
              chains = 2, iter = 1500, warmup = 750, cores = 2,
              refresh = 0, seed = 1, priors = NULL)
cat("=== AR posterior ===\n"); print(summary(fit$fit)$cor_pars)

pr <- et_predict(fit, newdata = newd, env_noise = list(x = 0.2),
                 n_draws = 1500, n_perturb = 500)
d <- decompose_uncertainty(pr)
cat("\n=== decomposition by lead time ===\n")
print(round(d[, c("obs_id","param_var","env_var","residual_var","temporal_var","total_var")], 4))

pp <- pr$posterior_predict
raw  <- apply(pp, 2, var)
wins <- apply(pp, 2, function(c){q<-quantile(c,c(.01,.99),names=FALSE);var(pmin(pmax(c,q[1]),q[2]))})
cat("\n=== raw vs winsorized pp variance ===\n")
print(round(data.frame(lead=1:H, raw=raw, winsor=wins, ratio=raw/wins), 3))

## CI width used for shelf life vs sqrt(total_var) reported in the budget
ci <- pr$credible_intervals; ci90 <- ci[ci$ci_level==0.9,]
cat("\n=== CI width vs budget-implied width ===\n")
print(round(data.frame(lead=1:H, ci_width=ci90$width,
                       implied=2*qnorm(.95)*sqrt(d$total_var),
                       ratio=ci90$width/(2*qnorm(.95)*sqrt(d$total_var))),3))
