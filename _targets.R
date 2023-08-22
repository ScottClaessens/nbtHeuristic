library(targets)
library(tarchetypes)
library(tidyverse)
source("R/functions.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("brms", "cowplot", "readxl", "tidybayes", "tidyverse"))

# full workflow
list(
  
  #### Power analysis ####
  
  # run power analysis
  tar_target(power, runPowerAnalysis(nSims = 2, nPerCell = 140)), # !! to change back to nSims = 100
  
  #### Model compilation ####
  
  # compile Stan models
  tar_target(compiledModel1, compileModel1(dSim)),
  tar_target(compiledModel2, compileModel2(dSim)),
  tar_target(compiledModel3, compileModel3(dSim)),
  
  #### Simulation ####
  
  # simulate data
  tar_target(dSim, simulateData()),
  # fit models
  tar_target(simModel1, update(compiledModel1, newdata = dSim, chains = 4, cores = 4, seed = 1)),
  tar_target(simModel2, update(compiledModel2, newdata = dSim, chains = 4, cores = 4, seed = 1)),
  tar_target(simModel3, update(compiledModel3, newdata = dSim, chains = 4, cores = 4, seed = 1)),
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
  
  #### Real data ####
  
  # load data
  tar_target(fileData, "data/nbtHeuristic_cleanData.csv", format = "file"),
  tar_target(d, read_csv(fileData)),
  # plot sample
  tar_target(plotSample, plotSampleProlific(d)),
  # fit models
  tar_target(m1, update(compiledModel1, newdata = d, chains = 4, cores = 4, seed = 1)),
  tar_target(m2, update(compiledModel2, newdata = d, chains = 4, cores = 4, seed = 1)),
  tar_target(m3, update(compiledModel3, newdata = d, chains = 4, cores = 4, seed = 1)),
  # plot results
  tar_target(plot1, plotModel1(d, m1, filename = "figures/analysis/model1.pdf")),
  tar_target(plot2, plotModel2(d, m2, filename = "figures/analysis/model2.pdf")),
  tar_target(plot3, plotModel3(d, m3, filename = "figures/analysis/model3.pdf")),
  # posterior samples
  tar_target(post1, posterior_samples(m1)),
  tar_target(post2, posterior_samples(m2)),
  tar_target(post3, posterior_samples(m3)),
  # posterior predictions and estimands
  tar_target(fitted2, fitted(m2, newdata = cond, summary = FALSE)),
  tar_target(estimand2.1, fitted2[,4] - fitted2[,1]),
  tar_target(estimand2.2, fitted2[,5] - fitted2[,2]),
  tar_target(estimand2.3, fitted2[,6] - fitted2[,3]),
  tar_target(estimand3.1, post3$`b_dgAmountGiven_frameCondControl:logdgPageSubmitSecs`),
  tar_target(estimand3.2, post3$`b_dgAmountGiven_frameCondNeed:logdgPageSubmitSecs`),
  tar_target(estimand3.3, post3$`b_dgAmountGiven_frameCondDebt:logdgPageSubmitSecs`),
  
  #### Print session info ####
  
  # print session info for reproducibility
  tar_target(sessionInfo, writeLines(capture.output(sessionInfo()), "sessionInfo.txt"))
  
)
