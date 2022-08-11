# load libraries
library(rethinking)
library(dagitty)
library(tidyverse)
library(xtable)

# load data
test_dataset <- read_csv('output/test_dataset.csv')

dat <- list(
  g = test_dataset$expression,
  G = test_dataset$mean_expression,
  r = test_dataset$score,
  D = as.numeric(as.factor(test_dataset$drug)),
  C = as.numeric(as.factor(test_dataset$cell_id))
  )

# eda

## run
m5 <- ulam(
  flist = alist(
    # r as normal
    r ~ dnorm(mu, sigma),
    # regress r on g and G
    mu <- beta1[D,C] * g + beta2[D,C] * G,
    # declare priors
    matrix[D, C]:beta1 ~ dnorm(0, 1),
    matrix[D, C]:beta2 ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = dat,
  chains = 4,
  cores = 4,
  seed = 1234
)

write_rds(m5, 'output/m5.rds')

## visualize
png('presentation/figures/trank_5.png')
trankplot(m5, ask = FALSE, max_rows = 1000)
trankplot(m5)
dev.off()

precis(m5, depth = 4)

# precis table
print(xtable(precis(m5, depth = 3)),
      booktabs = TRUE,
      file = 'presentation/tables/model_5.tex')

# plot the samples
plot(test_dataset$score ~ test_dataset$expression, type = 'n')
abline(lm(test_dataset$score ~ test_dataset$expression))

boxplot(test_dataset$score ~ as.numeric(as.factor(test_dataset$hallmark)))
boxplot(test_dataset$score ~ as.numeric(as.factor(test_dataset$cell_id)))
boxplot(test_dataset$score ~ as.numeric(as.factor(test_dataset$drug)))

# plot the effects
# the intercept
set.seed(1234)
alphas <- rnorm(1000, post5$alpha, post5$sigma)
dens(alphas,
     xlab = 'The intercept (alpha)')

## the beta1
set.seed(1234)
plot(NULL, xlim = c(-2, 2), ylim = c(0, 2),
     xlab = 'Effect of g on r', ylab = 'Density')
for (i in 1:5) {
  for (j in 1:4) {
    d <- rnorm(10000, post5$beta1[ , i, j], post5$sigma)
    dens(d, add = TRUE, col = j)
    abline(v = mean(d), lty = 2)
  }
}

## the beta2
set.seed(1234)
plot(NULL, xlim = c(-2, 2), ylim = c(0, 2),
     xlab = 'Effect of g on r', ylab = 'Density')
for (i in 1:5) {
  for (j in 1:4) {
    d <- rnorm(10000, post5$beta2[ , i, j], post5$sigma)
    dens(d, add = TRUE, col = j)
    abline(v = mean(d), lty = 2)
  }
}

## plot the predictions
gseq <- seq(from = min(dat$g), to = max(dat$g), len = 50)
Gseq <- seq(from = min(dat$G), to = max(dat$G), len = 50)
mu <- link(m5, data = list(g = gseq, gbar = mean(gseq), D = 2, C = 3, G = Gseq, Gbar = mean(Gseq)))
plot(dat$r ~ dat$g, type = 'n')
plot(gseq, apply(mu, 2, mean, rm.na = TRUE))

# marginalize