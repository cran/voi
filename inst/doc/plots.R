## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE,warning=FALSE----------------------------------------------
library(voi)
library(ggplot2)
library(dplyr)

## -----------------------------------------------------------------------------
head(chemo_evsi_or)

## -----------------------------------------------------------------------------
evpi_df <- evpi(outputs = chemo_cea_501)

## -----------------------------------------------------------------------------
evppi_df <- attr(chemo_evsi_or, "evppi")

## -----------------------------------------------------------------------------
ggplot(chemo_evsi_or, aes(x=k, y=evsi, group=n, color=n)) +
  geom_line() + 
  scale_colour_gradient(low="skyblue", high="darkblue") + 
  xlab("Willingness to pay (£)") + ylab("EVSI per person") + 
  geom_line(data=evpi_df, aes(x=k, y=evpi), 
            color="black", lwd=1.5, inherit.aes = FALSE) + 
  geom_line(data=evppi_df, aes(x=k, y=evppi), 
            color="darkblue", lwd=1.5, inherit.aes = FALSE) + 
  labs(color="Sample size") + 
  xlim(0,52000) +
  annotate("text",x=50000,y=125,label="EVPI",hjust=0) +
  annotate("text",x=50000,y=107,label="EVPPI",color="darkblue",hjust=0) +
  annotate("text",x=50000,y=50,label="EVSI",color="darkblue",hjust=0)

## -----------------------------------------------------------------------------
n_use <- c(250, 500, 1000)
chemo_evsi_or %>%
  filter(n %in% n_use) %>%
ggplot(aes(x=k, y=evsi, group=n, color=n)) +
  geom_line(lwd=1.5) + 
  scale_colour_binned(breaks=n_use, 
                      low = "gray80", high = "gray20",
                      guide=guide_legend(title="Sample size",
                                         reverse=TRUE)) + 
  xlab("Willingness to pay (£)") + ylab("EVSI (per person)") + 
  geom_line(data=evppi_df, aes(x=k, y=evppi), 
            color="darkblue", lwd=1.5, lty=3, inherit.aes = FALSE) + 
  annotate("text", x=50000, y=130, col="darkblue", label="EVPPI", hjust=1)

## -----------------------------------------------------------------------------
chemo_evsi_or %>% 
  filter(k %in% seq(10000,50000, by=1000)) %>%
ggplot(aes(x=n, y=evsi, group=k, color=k)) +
  geom_line() + 
  scale_colour_gradient(low="skyblue", high="darkblue") + 
  xlab("Sample size") + ylab("EVSI per person") + 
  labs(color="Willingness-to-pay (£)")

## -----------------------------------------------------------------------------
nbs <- enbs(chemo_evsi_or, 
            costs_setup = c(5000000, 10000000), 
            costs_pp = c(28000, 42000), 
            pop = 46000, 
            time = 10)
nbs %>% 
  filter(k==20000) %>%
  head()

## -----------------------------------------------------------------------------
nbs %>% 
  filter(k==20000) %>%
  ggplot(aes(x=n, y=enbs)) + 
  geom_line() + 
  xlab("Sample size") + ylab("Expected net benefit of sampling") + 
  scale_y_continuous(labels = scales::dollar_format(prefix="£")) 

## -----------------------------------------------------------------------------
nbs %>% 
  filter(k==20000) %>%
  mutate(q975 = qnorm(0.975, enbs, sd),
         q75 = qnorm(0.75, enbs, sd),
         q25 = qnorm(0.25, enbs, sd),
         q025 = qnorm(0.025, enbs, sd)) %>%
  ggplot(aes(y=enbs, x=n)) +
  geom_ribbon(aes(ymin=q025, ymax=q975), fill="gray80") +
  geom_ribbon(aes(ymin=q25, ymax=q75), fill="gray50") +
  geom_line() +
  xlab("Sample size") + ylab("Expected net benefit of sampling") + 
  scale_y_continuous(labels = scales::dollar_format(prefix="£")) + 
  annotate("text",x=1100,y=85000000,label="95% credible interval") +
  annotate("text",x=1250,y=70000000,label="50% credible interval")

## -----------------------------------------------------------------------------
nbs_subset <- nbs %>% 
 filter(k==20000, n %in% seq(100,800,by=100))
eopt <- enbs_opt(nbs_subset, smooth=TRUE, keep_preds=TRUE)

## -----------------------------------------------------------------------------
eopt

## -----------------------------------------------------------------------------
preds <- attr(eopt, "preds")
head(preds, 2)

## -----------------------------------------------------------------------------
ggplot(nbs_subset, aes(x=n, y=enbs)) + 
  geom_point() + 
  geom_line(data=preds) + 
  xlab("Sample size") + ylab("Expected net benefit of sampling") + 
  scale_y_continuous(labels = scales::dollar_format(prefix="£")) + 
  geom_vline(data=eopt, aes(xintercept=nlower), col="blue", lty=2) +
  geom_vline(data=eopt, aes(xintercept=nmax), col="blue", lty=2) +
  geom_vline(data=eopt, aes(xintercept=nupper), col="blue", lty=2) +
  geom_hline(data=eopt, aes(yintercept=enbsmax), col="blue", lty=2)

## -----------------------------------------------------------------------------
attr(nbs, "enbsmax") %>%
  filter(k %in% seq(20000,50000,by=10000))

## -----------------------------------------------------------------------------
p1 <- attr(nbs, "enbsmax") %>%
  ggplot(aes(y=nmax, x=k)) + 
  geom_ribbon(aes(ymin=nlower, ymax=nupper), fill="gray") + 
  geom_line() + 
  ylim(0,800) + 
  xlab("Willingness to pay (£)") + ylab("Optimal sample size") +
  ggtitle("Unsmoothed")
nbs_smoothed <- enbs(chemo_evsi_or, 
                     costs_setup = c(5000000, 10000000), 
                     costs_pp = c(28000, 42000), 
                     pop = 46000, time = 10, smooth=TRUE)
p2 <- attr(nbs_smoothed, "enbsmax") %>%
  ggplot(aes(y=nmax, x=k)) +
  geom_ribbon(aes(ymin=nlower, ymax=nupper), fill="gray") + 
  geom_line() + 
  ylim(0,800) + 
  xlab("Willingness to pay (£)") + ylab("Optimal sample size") +
  ggtitle("Smoothed")
gridExtra::grid.arrange(p1, p2, nrow=1)

## -----------------------------------------------------------------------------
nbs <- enbs(chemo_evsi_or %>% filter(k==20000, n==400), 
            costs_setup = c(5000000, 10000000), 
            costs_pp = c(28000, 42000), 
            pop = seq(0,60000,length.out=50), 
            time = seq(0,10,length.out=50))
ggplot(nbs, aes(y=pop, x=time, fill=pce)) + 
  geom_raster() + 
  labs(fill="Probability") +
  scale_fill_gradient(low="darkblue", high="white") + 
  xlab("Decision horizon (years)") + ylab("Incident population") +
  theme_minimal()

