# custom functions

# function for running power analysis with zoib model
runPowerAnalysis <- function(nSims = 100, nPerCell = 100, seed = 1) {
  # set seed
  set.seed(seed)
  # initialise zoib model
  f <- bf(
    y ~ 0 + condition,
    phi ~ 0 + condition,
    zoi ~ 0 + condition,
    coi ~ 0 + condition,
    family = zero_one_inflated_beta()
  )
  m <- brm(
    formula = f,
    data = data.frame(y = c(0.5, 0.5), condition = c("control", "treatment")),
    prior = c(
      prior(normal(0, 2), class = b),
      prior(normal(2, 1), class = b, dpar = phi),
      prior(normal(0, 2), class = b, dpar = zoi),
      prior(normal(0, 2), class = b, dpar = coi)
    ),
    chains = 0
  )
  # vector for sim results
  power <- rep(NA, nSims)
  # simulation
  for (i in 1:nSims) {
    # 1. simulate data - descriptive stats from Rand et al. 2012 (Study 6)
    d <- data.frame(
      y = c(
        rnorm(nPerCell, mean = 0.2088, sd = 0.1442),
        rnorm(nPerCell, mean = 0.2698, sd = 0.1406)
      ),
      condition = rep(c("TimeDelay","TimePressure"), each = nPerCell)
    )
    # 2. outcome bounded between 0 and 1
    d$y <- ifelse(d$y < 0, 0, d$y)
    d$y <- ifelse(d$y > 1, 1, d$y)
    # 3. fit zoib model
    model <- update(m, newdata = d, cores = 4, seed = seed)
    # 4. difference?
    pred <- 
      fitted(
        model,
        newdata = data.frame(condition = unique(d$condition)),
        summary = FALSE
      )
    diff <- pred[,1] - pred[,2]
    power[i] <- quantile(diff, 0.025) > 0 | quantile(diff, 0.975) < 0
  }
  # calculate overall power
  power <- mean(power)
  return(power)
}

# function to simulate data
simulateData <- function(N = 840, seed = 1) {
  # simulating from experiment DAG:
  # timeCond -> timeSpent -> dg
  # timeSpent <- frameCond -> dg
  # timeSpent <- U -> dg
  # sim seed
  set.seed(seed)
  # randomise participants to experimental conditions
  timeCond  <- sample(c("Pressure","Delay"), size = N, replace = TRUE)
  frameCond <- sample(c("Control","Need","Debt"), size = N, replace = TRUE)
  # simulate unobserved confound
  u <- rnorm(N)
  # set strength of confound
  k <- 0.2
  # time spent (lognormal)
  timeSpent <- 
    ifelse(
      # participants respond faster in the "pressure" condition and
      # participants also respond a little faster in the "need" condition
      timeCond == "Pressure",
      rlnorm(N, meanlog = ifelse(frameCond == "Need", 1.5, 1.8) + k*u, sdlog = 0.4),
      rlnorm(N, meanlog = ifelse(frameCond == "Need", 2.7, 3.0) + k*u, sdlog = 0.5)
    )
  # dictator game (zoib)
  zoi <- rbinom(N, prob = 0.30, size = 1) # probability zero or one
  coi <- rbinom(N, prob = 0.05, size = 1) # probability one, given 0/1
  mu <- ifelse(
    frameCond == "Control",
    # time spent is unrelated to giving in the control condition
    brms::inv_logit_scaled(rnorm(N, - 1 - k*u, 0.1)),
    ifelse(
      frameCond == "Need",
      # in the need condition, those who are faster give more
      brms::inv_logit_scaled(rnorm(N, - as.numeric(scale(log(timeSpent))) - k*u, 0.1)),
      # in the debt condition, those who are slower give more
      brms::inv_logit_scaled(rnorm(N, - 1 + as.numeric(scale(log(timeSpent))) - k*u, 0.1))
    )
  )
  phi <- 10
  dg <-
    ifelse(
      zoi == 1,
      # if zero-one inflation, choose zero or one
      ifelse(coi == 0, 0, 1),
      # else, simulate (0, 1)
      rbeta(N, shape1 = mu*phi, shape2 = (1-mu)*phi)
    )
  # full data frame
  out <- data.frame(id = 1:N, timeCond, frameCond, timeSpent, dg)
  return(out)
}

# compile model 1
compileModel1 <- function(d) {
  # model formula
  lnorm_model <- bf(
    timeSpent ~ 0 + timeCond:frameCond,
    sigma ~ 0 + timeCond:frameCond,
    family = lognormal()
  )
  # compile model
  out <- brm(
    formula = lnorm_model,
    data = d,
    prior = c(
      prior(normal(1.5, 0.5), class = b),
      prior(normal(-3, 1), class = b, dpar = sigma)
    ),
    chains = 0
  )
  return(out)
}

# compile model 2
compileModel2 <- function(d) {
  # model formula
  zoib_model <- bf(
    dg ~ 0 + timeCond:frameCond,
    phi ~ 0 + timeCond:frameCond,
    zoi ~ 0 + timeCond:frameCond,
    coi ~ 0 + timeCond:frameCond,
    family = zero_one_inflated_beta()
  )
  # compile model
  out <- brm(
    formula = zoib_model,
    data = d,
    prior = c(
      prior(normal(0, 2), class = b),
      prior(normal(2, 1), class = b, dpar = phi),
      prior(normal(0, 2), class = b, dpar = zoi),
      prior(normal(0, 2), class = b, dpar = coi)
    ),
    chains = 0
  )
  return(out)
}

# compile model 3
compileModel3 <- function(d) {
  # model formulae
  iv_model1 <- bf(dg ~ 0 + frameCond + frameCond:log(timeSpent))
  iv_model2 <- bf(log(timeSpent) ~ 0 + timeCond:frameCond)
  # compile model
  out <- brm(
    formula = iv_model1 + iv_model2 + set_rescor(TRUE),
    data = d,
    prior = c(
      prior(normal(0, 1), class = b, resp = dg),
      prior(normal(0, 1), class = b, resp = logtimeSpent)
    ),
    chains = 0
  )
  return(out)
}

# plot results of model 1
plotModel1 <- function(d, model1, filename) {
  # experimental conditions
  conditions <- expand_grid(
    timeCond = c("Delay","Pressure"),
    frameCond = c("Control","Need","Debt")
  )
  # get posterior predictions on outcome scale
  post <-
    fitted(
      model1,
      newdata = conditions,
      summary = FALSE
    )
  # plot results
  out <-
    conditions %>%
    mutate(post = apply(post, 2, function(x) as.data.frame(x))) %>%
    unnest(post) %>%
    ggplot() +
    geom_jitter(
      data = d,
      aes(x = timeCond, y = log(timeSpent)),
      size = 0.8,
      width = 0.2,
      colour = "lightgrey"
      ) +
    stat_interval(
      aes(x = timeCond, y = log(x)),
      size = 10,
      alpha = 0.5
      ) +
    facet_wrap(. ~ frameCond) +
    scale_colour_brewer() +
    scale_y_continuous(
      name = "Time spent making decision in seconds (log scale)",
      breaks = log(c(2.5, 5, 10, 20, 40, 80)),
      labels = function(x) exp(x)
      ) +
    xlab("Timing condition") +
    theme_classic() +
    theme(legend.position = "none")
  # save
  ggsave(filename = filename, plot = out, height = 4, width = 5)
  return(out)
}

# plot results of model 2
plotModel2 <- function(d, model2, filename) {
  # experimental conditions
  conditions <- expand_grid(
    timeCond = c("Delay","Pressure"),
    frameCond = c("Control","Need","Debt")
  )
  # get posterior predictions on outcome scale
  post <-
    fitted(
      model2,
      newdata = conditions,
      summary = FALSE
    )
  # plot results
  out <-
    conditions %>%
    mutate(post = apply(post, 2, function(x) as.data.frame(x))) %>%
    unnest(post) %>%
    ggplot() +
    geom_jitter(
      data = d,
      aes(x = timeCond, y = dg),
      size = 0.8,
      width = 0.2,
      colour = "lightgrey"
    ) +
    stat_interval(
      aes(x = timeCond, y = x),
      size = 10,
      alpha = 0.5
    ) +
    facet_wrap(. ~ frameCond) +
    scale_colour_brewer() +
    scale_y_continuous(
      name = "Proportion given in Dictator Game",
      limits = c(0, 1)
    ) +
    xlab("Timing condition") +
    theme_classic() +
    theme(legend.position = "none")
  # save
  ggsave(filename = filename, plot = out, height = 4, width = 5)
  return(out)
}

# plot results of model 3
plotModel3 <- function(d, model3, filename) {
  # plot
  out <-
    plot(
      conditional_effects(model3),
      points = TRUE,
      plot = FALSE
    )[[3]] +
    scale_y_continuous(
      name = "Proportion given in Dictator Game",
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      limits = c(-0.15, 1.15)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_continuous(name = "Time spent making decision in seconds") +
    guides(
      colour = guide_legend(title = "Framing"),
      fill = guide_legend(title = "Framing")
      ) +
    theme_classic()
  # move layers
  out$layers <-
    list(
      out$layers[[3]],
      out$layers[[4]],
      out$layers[[1]],
      out$layers[[2]]
    )
  # save
  ggsave(filename = filename, plot = out, height = 4, width = 5)
  return(out)
}
