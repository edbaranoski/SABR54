# SABR 54 scripts

require(stats)
require(tidyverse)
require(plotly)

# 2025 values
lgERA <- 4.15 
# runs_per_IP contains the likelihood for runs scored per inning and how many
# runs
runs_per_IP <- tibble(
  pcts = c(
    0.7468,
    0.1363,
    0.0645,
    0.0297,
    0.0130,
    0.00541,
    0.00242,
    0.000945,
    0.000373,
    0.000169,
    0.0000780,
    0.0000260,
    0.00000433,
    0.00000867,
    0.00000433
  ),
  values = c(0:14)
)
BB_pct <- 0.0843
HBP_pct <- 0.0106
K_pct <- 0.2229
single_pct <- 0.1432
double_pct <- 0.0425
triple_pct <- 0.0034
HR_pct <- 0.031
groudball_pct <- 0.284
flyball_pct <- 0.261

# sample parameters
innings <- 140 # sample innings
WHIPfactor <- 3+1.5  # PAs/inning, assuming WHIP of 1.5
TBF <-innings*WHIPfactor # total batters faced giving innings and WHIPfactor
sample_size <- 10000 # number of sample "seasons"

# ERA variance calculations:
# runs per inning has to be multiplied by 9 to get average contribution to nine-inning ERA
ERA <- binnedProbability(runs_per_IP$pcts,innings,trials=sample_size,
                         weights=9*runs_per_IP$values)

# FIP calculation:
# rate stats per PA need to be converted to rate per inning by multiplying by
# WHIPfactor. FIP is already scaled to ERA, so no furhter weighting is needed
FIP <- binnedProbability(c(BB_pct,HBP_pct,HR_pct,K_pct),TBF,trials=sample_size,
                         weights=WHIPfactor*c(3,3,13,-2))
FIP_c <- lgERA-mean(FIP$values) # FIP requires a constant value to scale to ERA
# note that any constant does not change the variances
FIP$values <- FIP$values+FIP_c

# wOBA calculation:
# weights have to be multiplied by the slope of the ERA vs wOBA curve to put it
# in the proper context. From 2021-2025, this is slope is about 13
ERA_wOBA_slope <-13
wOBA <- binnedProbability(c(BB_pct,HBP_pct,single_pct,double_pct,triple_pct,HR_pct),
                          TBF,trials=sample_size,
                          weights=ERA_wOBA_slope*c(0.701,0.733,0.895,1.27,1.607,2.066))

# K%-BB% calculation: 
# K%-BB% is inversely related to ERA.  The best fit has a slope and an intercept
# that have to be factored into the ERA estimate using K%-BB%.  These values
# were the least squares fit to pitcher data from 2021-2025
ERA_KBB_slope <- -9.67
ERA_KBB_intercept <- 5.47
KBB <- binnedProbability(c(BB_pct,K_pct),TBF,trials=sample_size,
                         weights=ERA_KBB_slope*c(-1,1),
                         added_constant=ERA_KBB_intercept)

# SIERA calculation:
# SIERA has non-linear terms which prevent the analytical calculation of
# variances but binnedProbability can still compute the relative components.  It
# also uses a new parameter netGB of an appropriately signed version of
# groundball percentage minus flyball percentage.
# SIERAprob can still contain the relative frequency for each component piece
# element of SIERA.  
# Note: the covariance matrices in SIERAprob and incomplete because since it
# does not include terms based on netGB
SIERAprob <- binnedProbability(c(HR_pct, BB_pct+HBP_pct, groudball_pct, flyball_pct),
                               TBF,trials=sample_size)
probs <- SIERAprob$count/TBF
netGB <- probs[3,]-probs[4,]
netGBsquared <- sign(netGB)*netGB^2
SIERA_values <- 6.145-16.986*probs[1,]+7.653*probs[1,]^2+11.434*probs[2,]+
  -1.858*netGB-6.664*netGBsquared+10.13*probs[1,]*netGB-5.195*probs[2,]*netGB
SIERA_values <- SIERA_values+(lgERA-mean(SIERA_values))

# wOBA-K% (in various ratios):
# This is trickier, because one has to balance the slope of the ERA to (wOBA-x
# K%)  for the desired value of x
K_wOBA_ratio <- 2
ERA_wOBA_K_slope <- 9
wOBA_K <- binnedProbability(c(BB_pct,HBP_pct,single_pct,double_pct,triple_pct,HR_pct,K_pct),
                          TBF,trials=sample_size,
                          weights=ERA_wOBA_slope*c(0.701,0.733,0.895,1.27,1.607,2.066,-K_wOBA_ratio))

