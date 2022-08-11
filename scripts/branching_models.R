# load libraries
library(rethinking)
library(dagitty)
library(tidyverse)
library(xtable)

# simulate
N <- 1000
set.seed(123)

g <- rgamma(N, 2, 1)
dens(g)
r <- c(rnorm(N/2, 1, .5), rnorm(N/2, -1, .5))
dens(r)
r <- r - .2 * g
dens(r)

plot(r ~ g)
abline(lm(r ~ g), col = 'red')

m1 <- ulam(
  flist = alist(
    r|r > 0 ~ dnorm(mu, sigma),
    r|r < 0 ~ dnorm(mu, sigma),
    mu <- beta * g,
    beta ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = list(
    r = r,
    g = g
  ),
  chains = 4,
  cores = 4
)

precis(m1)
post1 <- extract.samples(m1)

set.seed(1234)
beta <- rnorm(10000, post1$beta, post1$sigma)
dens(beta)

# branching model
## simulate
g <- rpois(N, 20)
dens(g)

set.seed(1234)
GI <- sample(c(1,2,3), N, replace = TRUE)
g <- rpois(N, 5) + 2 * GI
G <- aggregate(g, list(G = GI), mean)$x[GI]
r <- rnorm(N) - .5 * g - G

# replicated for 3 Cs
g <- rep(g, 3*3)
GI <- rep(GI, 3*3)
G <- rep(G, 3*3)
r <- rep(r, 3*3)
D <- rep(rep(c(1,2,3), each = N), 3)
C <- rep(c(1,2,3), each = N*3)

# introduce the effect of D, module C
r <- ifelse(D == GI, r * - .5 * D * C, r)

boxplot(r ~ paste(D, GI, C))

## run
m5 <- ulam(
  flist = alist(
    # r as normal
    r ~ dnorm(mu, sigma),
    # regress r on g and G
    mu <- alpha + beta1[D,C] * (g - gbar) + beta2[D,C] * (G - Gbar),
    # declare priors
    alpha ~ dnorm(0, 1),
    matrix[D, C]:beta1 ~ dnorm(0, 1),
    matrix[D, C]:beta2 ~ dnorm(0, 1),
    sigma ~ dunif(0, 10)
  ),
  data = list(
    N = N,
    g = g,
    G = G,
    r = r,
    D = D,
    C = C,
    gbar = mean(g),
    Gbar = mean(G)
  ),
  chains = 4,
  cores = 4,
  seed = 1234
)

## visualize
png('presentation2/figures/trank_5.png')
trankplot(m5, ask = FALSE, max_rows = 1000)
dev.off()

precis(m5, depth = 3)

print(xtable(precis(m5, depth = 3)),
      booktabs = TRUE,
      file = 'presentation2/tables/model_5.tex')

post5 <- extract.samples(m5)


png('presentation2/figures/model_5.png',
    height = 4, width = 7, units = 'in', res = 300)

par(mfrow = c(1, 2))

# plot(graphLayout(dag5))

boxplot(r ~ paste(D, GI, C, sep = '/'), lwd = 2,
        xlab = 'Drug effect (D), by Group (G) & Cell (C)',
        ylab = 'response (r)',
        main = 'ifelse(D == G,r*-.5*D*C,r)')

set.seed(1234)
D1G1C1 <- link(m5, data = list(D = 1, G = 1, C = 1, gbar = mean(g), Gbar = mean(g)))
D1G1C2 <- link(m5, data = list(D = 1, G = 1, C = 2, gbar = mean(g), Gbar = mean(g)))
D1G1C3 <- link(m5, data = list(D = 1, G = 1, C = 3, gbar = mean(g), Gbar = mean(g)))

dens(D1G1C3, lwd = 2,
     ylim = c(0, .1),
     xlim = c(0, 60),
     main = 'Post. Pred. Dist.',
     xlab = 'Effect D1 on G1, by C')
dens(D1G1C2, lwd = 2, col = 2, add = TRUE)
dens(D1G1C1, lwd = 2, col = 4, add = TRUE)

text(mean(D1G1C3), .09, 'C3')
text(mean(D1G1C2), .09, 'C2', col = 2)
text(mean(D1G1C1), .09, 'C1', col = 4)

dev.off()
