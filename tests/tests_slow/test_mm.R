
## single parameter

## multi parameter (2)

## multi parameter (convoi1 chemo example)
## todo get the names consistent 


load_all("../../../EVSI")
load_all("../..")

## gen.quantiles works straightforwardly - works by randomly permuting each 1D quantile set. appends sample sizes, to hopefully cover the space of all parameters and sample sizes.  returns a data frame.

## match names with those used in test
## amb is state 1 
## hosp is state 2
## rec is state 3 
## death is state 4 

library(dplyr)
inputs_mm  <- chemo_pars
poi <- c("gamma.hosp", "gamma.dead", "lambda.hosp.rec.TH", "lambda.home.rec.TH","pi1","pi2")
Q <- 50

## TODO break down the evsi_moment_matching function
## using argument formats drawn from chemo example in voi package
## attach other packages as needed 

decision_model_cea <- list(effect=chemo_cea$e, cost=chemo_cea$c, wtp=chemo_cea$k)

evsi_specification <- list(evppi_parameter = c("pi1","pi2",
                                             "gamma.hosp","gamma.dead", 
                                             "lambda.home.rec.TH", "lambda.hosp.rec.TH"), 
                           nested_sample_number = 50,                                     
                           quantile_parameter = c("pi1","pi2",
                                                "gamma.hosp","gamma.dead", 
                                                "lambda.home.rec.TH", "lambda.hosp.rec.TH"))

## TODO explain quantile_parameters. why different from evppi_parameter? 

decision_model_psa <- inputs_mm
decision_model_mcmc <- list(model = model.data.full,
                            prior_data = chemo_auxdata,
                            iter_number = 10000,
                            burnin_number = 3000)
evsi_specification

## (1) Estimate EVPPI

evppi.full <- evppi(chemo_cea, inputs_mm, poi)
evppi.full <- c(0.614339708525307, 8.26902009850676, 32.6788861715972, 46.3974562401871, 
27.3720236731283) # five wtps from 1k to 5k 

## (2) Extract Q (done)
## (3) Generate quantiles 

quants <- EVSI::gen.quantiles(poi, decision_model_psa, Q=Q, N=c(5,1000))

## (4a) Simulate future data 

data.list.full <- list()
for(i in 1:Q){
    sim_data <- rfn(quants[i,], quants[i,"N"])
    data.list.full[[i]] <- sim_data
}

## (4b) Fit Bayesian model to it
## this is all in mm.post.var function. generate args for this

filein <- file.path(tempdir(), fileext="datmodel.txt")
R2OpenBUGS::write.model(decision_model_mcmc$model, filein)
model.stats <- "mm_jags_model.txt"
# data <- data.list.full
N.name <- "n"
N.size <- quants[,"N"]
effects <- chemo_effects_fn
costs <- chemo_costs_fn
data.stats <- decision_model_mcmc$prior_data

argument.names <- unique(c(names(formals(effects)), names(formals(costs))))
## Model parameters that are updated in Bayesian model
monitor <- names(which(sapply(argument.names, grep.fun, model.file.text = readLines(model.stats)) > 0))

## what is "when prior predictive sampling is required???".  debug() to see if this is done.  no.  assume we don't need it for the moment. 

    ## Calculate costs and effects for all simulations
    ## FAFF WITH ARGS HERE. looks better done with do.call 
    ## so what does find.args do.  It takes a df of BUGS outputs and forms list of args to supply to effects(), 
    ## aha also faff with args named eg SE[1] SE[2]
    ## TH is missing in this case -- should just rewrite to use list 

    ## So TODO decide on package formats for model functions and updating functions 
    ## currently accepts one PSA iteration, args could be vectors, scalars, matrices
    ## names should match with poi

    ## Then need function to draw mcmc sample of those pars given simulated data
    ## could convert to list of lists, then loop over outer list below 
    ## is there a utility in rstan or suchlike for this?

## (4c) Rerun PSA using new posterior to update NB dist 

var.prepost<-list()
for(q in 1:Q){
    ## Both data sets must be given to jags
    data.full<-append(data.list.full[[q]], data.stats)
    ss <- list(N.size[q])
    names(ss) <- N.name
    data.full <- append(data.full, ss)
    ## Sample from Bayesian updating
    sample.dat <- running.jags(model.stats, data.full, monitor, n.burnin=3000, n.iter=10000, thin=1)   ## eh why so quick? 
    
    TH <- chemo_auxdata$TH
    n.comparisons <- 2
    incremental.benefit <- matrix(NA, nrow=dim(sample.dat)[1], ncol=2*n.comparisons)
    for(i in 1:dim(sample.dat)[1]){
        formals(effects) <- EVSI::find.args(effects, sample.dat, monitor, i)
        formals(costs) <- EVSI::find.args(costs, sample.dat, monitor, i)
        model.effects <- eval(effects(), envir=parent.frame())
        model.costs <- eval(costs(), envir=parent.frame())
        ref <- 1
        incremental.benefit[i,] <- cbind(model.effects[-ref] - model.effects[ref], model.costs[- ref] - model.costs[ref])
    }
    print(paste("Update", q, "completed"))
    var.prepost[[q]] <- var(incremental.benefit)
}
EVSI.var.full<-list(variance.Q = var.prepost, # list of 50 4x4 matrices 
                    N.size = N.size,
                    evppi = evi,
                    he = he,
                    update = update.func)

### (5): Calculate the sigma^2 = Var[INB_t(theta_s)] - (1/Q)*sum(sigma_q^2)
### Estimate the posterior variance across multiple potential datasets

### (6): Rescale the simulations of eta_t(s) such that their variance is equal to sigma^2;
### (7): Estimates EVSI by using rescaled simulations 

  evsi.heath <- EVSI::evsi.calc(EVSI.var.full)

rfn <- function(quants, n){ 
    X.SE1 <- rbinom(1, n, quants[,"pi1"])
    X.SE2 <- rbinom(1, n, quants[,"pi2"])
    X.N.hosp <- rbinom(1, sum(X.SE1,X.SE2), quants[,"gamma.hosp"])
    X.N.dead <- rbinom(1, X.N.hosp, quants[,"gamma.dead"])
    
    N.home <- sum(X.SE1,X.SE2) - X.N.hosp
    recover.home <- -log(1 - quants[,"lambda.home.rec.TH"])
    T.rec.home <- rexp(N.home, recover.home)
    
    N.hosp <- X.N.hosp - X.N.dead
    recover.hosp <- -log(1 - quants[,"lambda.hosp.rec.TH"])
    T.rec.hosp <- rexp(N.hosp, recover.hosp)
    
    data.n <- list(X.SE1 = X.SE1, X.SE2 = X.SE2,
                   X.N.hosp = X.N.hosp, X.N.dead = X.N.dead,
                   N.home = N.home, T.rec.home = T.rec.home, 
                   N.hosp = N.hosp, T.rec.hosp = T.rec.hosp)
}  
