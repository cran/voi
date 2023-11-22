## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
set.seed(1) 
nsam <- 10000
inputs <- data.frame(p1 = rnorm(nsam, 1, 1), 
                     p2 = rnorm(nsam, 0, 2))

## -----------------------------------------------------------------------------
outputs_nb <- data.frame(t1 = 0, 
                         t2 = inputs$p1 - inputs$p2)

## -----------------------------------------------------------------------------
outputs_cea <- list( 
  e = data.frame(t1 = 0, t2 = inputs$p1), 
  c = data.frame(t1 = 0, t2 = inputs$p2), 
  k = c(1, 2, 3)
)

## -----------------------------------------------------------------------------
decision_current <- 2
nb_current <- 1
decision_perfect <- ifelse(outputs_nb$t2 < 0, 1, 2)
nb_perfect <- ifelse(decision_perfect == 1, 0, outputs_nb$t2)
(evpi1 <- mean(nb_perfect) - nb_current)

## -----------------------------------------------------------------------------
opp_loss <- nb_perfect - nb_current
mean(opp_loss)

## -----------------------------------------------------------------------------
library(voi)
evpi(outputs_nb)
evpi(outputs_cea)

## -----------------------------------------------------------------------------
prob_correct <- 1 - pnorm(0, 1, sqrt(5))

## -----------------------------------------------------------------------------
mean_truncnorm <- function(mu, sig, lower=-Inf, upper=Inf){ 
  a <- (lower-mu)/sig
  b <- (upper-mu)/sig
  mu + sig * (dnorm(a) - dnorm(b)) / (pnorm(b) - pnorm(a))
}
enb_correct <- mean_truncnorm(1, sqrt(5), lower=0) 
mean_nb_perfect <- enb_correct * prob_correct
(evpi_exact <- mean_nb_perfect - nb_current)

## -----------------------------------------------------------------------------
evppi(outputs_nb, inputs, pars="p1")
evppi(outputs_nb, inputs, pars=c("p1","p2"))

## -----------------------------------------------------------------------------
evppi(outputs_nb, inputs, pars=list("p1","p2"))
evppi(outputs_nb, inputs, pars=list("p1",c("p1","p2")))

## -----------------------------------------------------------------------------
evppi(outputs_cea, inputs, pars=list("p1",c("p1","p2")))

## -----------------------------------------------------------------------------
evppi(outputs_nb, inputs, pars="p1", method="gp", nsim=1000)

## -----------------------------------------------------------------------------
evppi(outputs_nb, inputs, pars="p1", method="earth")

## ----eval=FALSE---------------------------------------------------------------
#  install.packages('INLA', repos='https://inla.r-inla-download.org/R/stable')
#  install.packages('splancs')

## ----eval=FALSE---------------------------------------------------------------
#  evppi(outputs_nb, inputs, pars=c("p1","p2"), method="inla", pfc_struc="iso")

## ----cache=TRUE,message=FALSE-------------------------------------------------
evppi(chemo_nb, chemo_pars, pars=colnames(chemo_pars), method="bart")
evpi(chemo_nb)

## -----------------------------------------------------------------------------
evppi(outputs_nb, inputs, pars=c("p1","p2"), method="gam", gam_formula="s(p1) + s(p2)")

## -----------------------------------------------------------------------------
evppi(outputs_nb, inputs, pars="p1", se=TRUE, B=100)

## -----------------------------------------------------------------------------
evppi(outputs_nb, inputs, pars="p1", n.blocks=20, method="so")

## -----------------------------------------------------------------------------
evppi(outputs_nb, inputs, pars="p1", method="sal")

## -----------------------------------------------------------------------------
model_fn_nb <- function(p1, p2){ 
  c(0, p1 - p2) 
}

## -----------------------------------------------------------------------------
model_fn_cea <- function(p1, p2){ 
  rbind(e = c(0, p1), 
        c = c(0, p2)) 
}

## -----------------------------------------------------------------------------
par_fn <- function(n){
  data.frame(p1 = rnorm(n, 1, 1),
             p2 = rnorm(n, 0, 2))
}

## ----eval=FALSE---------------------------------------------------------------
#  evppi_mc(model_fn_nb, par_fn, pars="p1", ninner=1000, nouter=100)

## ----eval=FALSE---------------------------------------------------------------
#  par_fn_corr <- function(n, p1=NULL){
#    p1_new <- if (is.null(p1)) rnorm(n, 1, 1) else p1
#    data.frame(p1 = p1_new,
#               p2 = rnorm(n, p1_new, 2))
#  }
#  evppi_mc(model_fn_nb, par_fn_corr, pars="p1", ninner=1000, nouter=100)

## -----------------------------------------------------------------------------
datagen_normal <- function(inputs, n=100, sd=1){
  data.frame(xbar = rnorm(nrow(inputs),
                          mean = inputs[,"p1"],
                          sd = sd / sqrt(n)))
}

set.seed(1)
evsi(outputs_nb, inputs, datagen_fn = datagen_normal, n=c(10,100,1000))

## -----------------------------------------------------------------------------
evsi(outputs_nb, inputs, study = "normal_known", n=c(100,1000), pars = "p1")
evsi(outputs_cea, inputs, study = "normal_known", n=c(100,1000), pars = "p1")

## -----------------------------------------------------------------------------
likelihood_normal <- function(Y, inputs, n=100, sig=1){
  mu <- inputs[,"p1"]
  dnorm(Y[,"xbar"], mu, sig/sqrt(n))
}

evsi(outputs_nb, inputs, datagen_fn = datagen_normal, likelihood = likelihood_normal, 
     n=100, pars = "p1", method="is", nsim=1000)

## -----------------------------------------------------------------------------
evsi(outputs_nb, inputs, study = "normal_known", n=100, pars = "p1", method="is", nsim=1000)

## ----message=FALSE------------------------------------------------------------
evsi(outputs_nb, inputs, study = "normal_known", n=10000, pars = "p1", method="mm", Q=30,
     model_fn = model_fn_nb, par_fn = par_fn,
     analysis_args = list(prior_mean=1, prior_sd=1, sampling_sd=1, niter=1000))

## -----------------------------------------------------------------------------
analysis_fn <- function(data, args, pars){
  dat <- list(y=c(data[,"y1"], data[,"y2"]))
  design <- list(n = rep(args$n, 2))
  priors <- list(a1=53, b1=60, mu=log(0.54), sigma=0.3)
  jagsdat <- c(dat, design, priors)
  or_jagsmod <- "
  model {
    y[1] ~ dbinom(p[1], n[1])
    y[2] ~ dbinom(p[2], n[2])
    p[1] <- p1
    p[2] <- odds[2] / (1 + odds[2])
    p1 ~ dbeta(a1, b1)
    odds[1] <- p[1] / (1 - p[1])
    odds[2] <- odds[1] * exp(logor)
    logor ~ dnorm(mu, 1/sigma^2)
  }
  "
  or.jag <- rjags::jags.model(textConnection(or_jagsmod), 
                              data=jagsdat, inits=list(logor=0, p1=0.5), quiet = TRUE)
  update(or.jag, 100, progress.bar="none")
  sam <- rjags::coda.samples(or.jag, c("logor"), 500, progress.bar="none")
  data.frame(logor_side_effects = as.numeric(sam[[1]][,"logor"]))
}

## -----------------------------------------------------------------------------
datagen_fn <- function(inputs, n=100){
  p1 <- inputs[,"p_side_effects_t1"]
  logor <- inputs[,"logor_side_effects"]
  odds1 <- p1 / (1 - p1)
  odds2 <- odds1 * exp(logor)
  p2 <- odds2 / (1 + odds2)
  nsim <- nrow(inputs)
  data.frame(y1 = rbinom(nsim, n, p1),
             y2 = rbinom(nsim, n, p2))
}

## ----cache=TRUE,message=FALSE,eval=requireNamespace("rjags")------------------
ev <- evsi(outputs=chemo_nb, inputs=chemo_pars, 
           method="mm",
           pars="logor_side_effects", 
           pars_datagen = c("p_side_effects_t1", "logor_side_effects"), 
           datagen_fn = datagen_fn, analysis_fn = analysis_fn, 
           n = 100, Q = 10, 
           model_fn = chemo_model_lor_nb, par_fn = chemo_pars_fn)

## -----------------------------------------------------------------------------
p1 <- rbeta(10000, 5, 95)

## ----fig.width=7,fig.height=5-------------------------------------------------
beta <- rnorm(10000, 0.8, 0.4)
p2 <- plogis(qlogis(p1) + beta)
plot(density(p1), lwd=2, xlim=c(0,1), main="")
lines(density(p2), col="red", lwd=2)
legend("topright", col=c("black","red"), lwd=c(2,2), 
       legend=c("Surveyed infection probability", "True infection probability"))

## -----------------------------------------------------------------------------
var(p2)

## -----------------------------------------------------------------------------
inputs <- data.frame(p1, beta)
(evppi_beta <- evppivar(p2, inputs, par="beta"))
(evppi_p1 <- evppivar(p2, inputs, par="p1"))

## -----------------------------------------------------------------------------
sqrt(var(p2)) # or sd(p2)
sqrt(var(p2) - evppi_beta$evppi)
sqrt(var(p2) - evppi_p1$evppi)

## ----fig.width=7,fig.height=5-------------------------------------------------
plot(x=p1, y=p2, pch=".")
mod <- mgcv::gam(p2 ~ te(p1, bs="cr"))
p1fit <- fitted(mod)
lines(sort(p1), p1fit[order(p1)], col="blue")

## -----------------------------------------------------------------------------
p1res <- p2 - p1fit
var(p2) - var(p1res)

## -----------------------------------------------------------------------------
var(p1fit)

## -----------------------------------------------------------------------------
evsivar(p2, inputs, study = "binary", pars="p1", n=c(100,1000,10000))

## -----------------------------------------------------------------------------
inputs_p2 = data.frame(p2 = p2)
evsivar(p2, inputs=inputs_p2, study = "binary", pars="p2", n=c(100, 1000, 10000))

