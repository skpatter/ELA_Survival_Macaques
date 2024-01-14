

library(brms)
library(dplyr)

ddAdult <- read.csv("ddAdult_June_github.csv", header = TRUE, na.strings = "") 
ddEarly <- read.csv("ddEarly_June_github.csv", header = TRUE, na.strings = "") 

## survival models ##
UpdatedPrior_AdultYr = c(set_prior("normal(12,0.4)", class = "Intercept"),
                         set_prior("normal(0,1)", class = "b"),
                         set_prior("cauchy(0,2)", class = "sd"))

# Adult survival #

# each form of early life adversity
cayo_adults1_fullrank <- brm(timeYr | cens(censored) ~ (1 | Birth.Season) + (1 | MaternalID.text) +
                               Sex * (Maternal_loss +  High_temperatureLog + 
                                        RainfallLog + Group_sizeLog + Maternal_kin_countLog + 
                                        Matriline_rank + First_born + CompetingSib + Hurricane)   ,
                             data = ddAdult, prior=UpdatedPrior_AdultYr,
                             iter = 2000, warmup = 1000, chains = 3,
                             family = weibull() )
# rerun with ordinal/continuous rank for supplemental
cayo_adults1_fullrank_ordinal <- brm(timeYr | cens(censored) ~ (1 | Birth.Season) + (1 | MaternalID.text) +
                                       Sex * (Maternal_loss +  High_temperatureLog + 
                                                RainfallLog + Group_sizeLog + Maternal_kin_countLog + 
                                                MatrilineRank + First_born + CompetingSib + Hurricane)   ,
                                     data = ddAdult, prior=UpdatedPrior_AdultYr,
                                     iter = 2000, warmup = 1000, chains = 2,
                                     family = weibull() )

# cumulative early life adversity index
cayo_adults_index <- brm(timeYr | cens(censored) ~ (1 | Birth.Season) + (1 | MaternalID.text) +
                                    Sex * (Cumulative_early_life_adversity) ,
                                  data = ddAdult, prior=UpdatedPrior_AdultYr,
                                  iter = 2000, warmup = 1000, chains = 3,
                                  family = weibull() )

# early life survival #
UpdatedPrior_EarlyIndex2 = c(set_prior("normal(1,0.1)", class = "Intercept"), 
                             set_prior("normal(0,.5)", class = "b"))


# each form of early life adversity
cayo_early <- brm(ExitFourYr | cens(CensorFour) ~ (1 | Birth.Season) + (1 | MaternalID.text) +
                             Sex * (Maternal_loss +  High_temperatureLog + 
                                      RainfallLog + Group_sizeLog + Maternal_kin_countLog + 
                                      Matriline_rank + First_born + Hurricane) ,
                           data = ddEarly, prior=UpdatedPrior_EarlyIndex2,
                           iter = 2000, warmup = 1000, chains = 3,
                           family = weibull() )

# rerun with ordinal/continuous rank for supplemental
cayo_early_ordinal <- brm(ExitFourYr | cens(CensorFour) ~ (1 | Birth.Season) + (1 | MaternalID.text) +
                                     Sex * (Maternal_loss +  High_temperatureLog + 
                                              RainfallLog + Group_sizeLog + Maternal_kin_countLog + 
                                              MatrilineRank + First_born + Hurricane) ,
                                   data = ddEarly, prior=UpdatedPrior_EarlyIndex2,
                                   iter = 2000, warmup = 1000, chains = 2,
                                   family = weibull() )

# cumulative early life adversity index
cayo_early_index <- brm(ExitFourYr | cens(CensorFour) ~ (1 | Birth.Season) + (1 | MaternalID.text) +
                                  Sex * (Cumulative_early_life_adversity) ,
                                data = ddEarly,prior=UpdatedPrior_EarlyIndex2,
                                iter = 2000, warmup = 1000, chains = 3,
                                family = weibull() )

# model comparisons #

cayo_adults1_fullrank_waic <- add_criterion(cayo_adults1_fullrank, "loo")
cayo_adults_index_fullrank_waic <- add_criterion(cayo_adults_index_fullrank, "loo")

adultWAIC <- loo_compare(cayo_adults1_fullrank_waic,
                         cayo_adults_index_fullrank_waic,
                         criterion = "loo")

cayo_early_fullrank_waic <- add_criterion(cayo_early_fullrank, "loo")
cayo_early_indexfullrank_waic <- add_criterion(cayo_early_indexfullrank, "loo")

earlyWAIC <- loo_compare(cayo_early_fullrank_waic,
                         cayo_early_indexfullrank_waic,
                         criterion = "loo")


# maternal social connectedness #

dsocAdult <- read.csv("dsocAdult_june_github.csv", header = TRUE, na.strings = "") 
dsocEarly <- read.csv("dsocEarly_june_github.csv", header = TRUE, na.strings = "")

UpdatedPrior_EarlyCSI = c(set_prior("normal(1,0.1)", class = "Intercept"),
                          set_prior("normal(0,.5)", class = "b"),
                          set_prior("cauchy(0,2)", class = "sd"))

cayo_early_CSI <- brm(ExitFourYr | cens(CensorFour) ~ (1 | Birth.Season) + (1|MaternalID.text) +
                             Sex * s_csiff ,
                           data = dsocEarly, prior=UpdatedPrior_EarlyCSI,
                           iter = 3000, warmup = 1000, chains = 3,
                           family = weibull() )

cayo_adult_CSI <- brm(ExitAgeYr | cens(censored) ~ (1 | Birth.Season) + (1|MaternalID.text) +
                             Sex * s_csiff ,
                           data = dsocAdult, prior=UpdatedPrior_AdultCSIYr,
                           iter = 3000, warmup = 1000, chains = 3,
                           family = weibull() )



#####################
## Pedigree Models ##
#####################

# dd = data
# A = relatedness matrix

# with pedigree:
mped2early <- brm(ExitFourYr | cens(CensorFour) ~  Cumulative_early_life_adversity + 
                    (1|gr(AnimalID.text, cov=A)) ,
                  data = dd,
                  data2 = list(A = A),
                  control=list(adapt_delta=0.99),
                  family = weibull() )

# without pedigree
m2early <- brm(ExitFourYr | cens(CensorFour) ~  Cumulative_early_life_adversity  ,
               data = dd,
               control=list(adapt_delta=0.99),
               iter = 1000, warmup = 100, chains = 2,
               family = weibull() )


