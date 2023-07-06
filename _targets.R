library(targets)
library(tarchetypes)
library(tidyverse)
source("R/functions.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("brms", "tidybayes", "tidyverse"))

# full workflow
list(
  
  #### Power analysis ####
  
  # run power analysis
  tar_target(power, runPowerAnalysis(nSims = 100, nPerCell = 140)),
  
  #### Model compilation ####
  
  # compile Stan models
  tar_target(compiledModel1, compileModel1(dSim)),
  tar_target(compiledModel2, compileModel2(dSim)),
  tar_target(compiledModel3, compileModel3(dSim)),
  
  #### Simulation ####
  
  # simulate data
  tar_target(dSim, simulateData()),
  # fit models
  tar_target(simModel1, update(compiledModel1, newdata = dSim,chains = 4, cores = 4, seed = 1)),
  tar_target(simModel2, update(compiledModel2, newdata = dSim,chains = 4, cores = 4, seed = 1)),
  tar_target(simModel3, update(compiledModel3, newdata = dSim,chains = 4, cores = 4, seed = 1)),
  # plot results
  tar_target(plotSim1, plotModel1(dSim, simModel1, filename = "figures/simulation/simModel1.pdf")),
  tar_target(plotSim2, plotModel2(dSim, simModel2, filename = "figures/simulation/simModel2.pdf")),
  tar_target(plotSim3, plotModel3(dSim, simModel3, filename = "figures/simulation/simModel3.pdf")),
  # posterior samples
  tar_target(postSim1, posterior_samples(simModel1)),
  tar_target(postSim2, posterior_samples(simModel2)),
  tar_target(postSim3, posterior_samples(simModel3)),
  # posterior predictions and estimands
  tar_target(cond, expand_grid(timeCond = c("Delay","Pressure"), frameCond = c("Control","Need","Debt"))),
  tar_target(fittedSim2, fitted(simModel2, newdata = cond, summary = FALSE)),
  tar_target(estimandSim2.1, fittedSim2[,4] - fittedSim2[,1]),
  tar_target(estimandSim2.2, fittedSim2[,5] - fittedSim2[,2]),
  tar_target(estimandSim2.3, fittedSim2[,6] - fittedSim2[,3]),
  tar_target(estimandSim3.1, postSim3$`b_dg_frameCondControl:logtimeSpent`),
  tar_target(estimandSim3.2, postSim3$`b_dg_frameCondNeed:logtimeSpent`),
  tar_target(estimandSim3.3, postSim3$`b_dg_frameCondDebt:logtimeSpent`),
  
  #### Print session info ####
  
  # print session info for reproducibility
  tar_target(sessionInfo, writeLines(capture.output(sessionInfo()), "sessionInfo.txt"))
  
)
