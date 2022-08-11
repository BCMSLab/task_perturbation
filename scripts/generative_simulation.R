# load libraries
library(rethinking)
library(dagitty)
library(tidyverse)
library(xtable)

# model 1: r on g ----
# dag
dag1 <- dagitty("dag{ g -> r }")

## simulate
N <- 100

set.seed(1234)
g <- rpois(N, 5)
r <- rnorm(N) - .5 * g

## run
m1 <- ulam(
  flist = alist(
    # r as normal
    r ~ dnorm(mu, sigma),
    # regress r on g
    mu <- alpha + beta * (g - gbar),
    # declare priors
    alpha ~ dnorm(0, 1),
    beta ~ dnorm(0, 1),
    sigma ~ dunif(0, 5)
  ),
  data = list(
    N = N,
    g = g,
    r = r,
    gbar = mean(g)
  ),
  chains = 4,
  cores = 4,
  seed = 1234
)

## visualize
png('presentation2/figures/trank_1.png')
trankplot(m1)
dev.off()

precis(m1)

print(xtable((precis(m1))),
      booktabs = TRUE,
      file = 'presentation2/tables/model_1.tex')

post1 <- extract.samples(m1)


png('presentation2/figures/model_1.png',
    height = 4, width = 7, units = 'in', res = 300)

par(mfrow = c(1, 2))

# plot(graphLayout(dag1))

plot(r ~ g,
     xlab = 'gene (g)',
     ylab = 'response (r)',
     main = 'r = rnorm(N) - .5*g',
     lwd = 2)
abline(lm(r ~ g))

plot(r ~ g, col = 'white',
     xlab = 'gene (g)',
     ylab = 'response (r)',
     main = 'Post. Pred. Dist.')
gseq <- seq(from = min(g), to = max(g), len = 100) 
set.seed(123)
mu <- link(m1, data = list(g = gseq, gbar = mean(g)))
lines(gseq, apply(mu, 2, mean))
shade(apply(mu, 2, PI, prob = .99), gseq)

set.seed(123)
rsim <- sim(m1, data = list(g = gseq, gbar = mean(g)))
shade(apply(rsim, 2, PI, prob = .95), gseq)

dev.off()


# model 2: indexing on G ----
# dag
dag2 <- dagitty("dag{ g -> r g -> G G -> r }")

## simulate
N <- 100

set.seed(1234)
G <- sample(c(1,2,3), N, replace = TRUE)
table(G)
g <- rpois(N, 5) + 2 * G
boxplot(g ~ G)
r <- rnorm(N) - .5 * g - G
boxplot(r ~ G)
## run
m2 <- ulam(
  flist = alist(
    # r as normal
    r ~ dnorm(mu, sigma),
    # regress r on g
    mu <- alpha[G] + beta * (g - gbar),
    # declare priors
    # indexing alpha on G
    alpha[G] ~ dnorm(0, 1),
    beta ~ dnorm(0, 1),
    sigma ~ dunif(0, 5)
  ),
  data = list(
    N = N,
    g = g,
    G = G,
    r = r,
    gbar = mean(g)
  ),
  chains = 4,
  cores = 4,
  seed = 1234
)

## visualize
png('presentation2/figures/trank_2.png')
trankplot(m2)
dev.off()

precis(m2, depth = 3)

print(xtable(precis(m2, depth = 3)),
      booktabs = TRUE,
      file = 'presentation2/tables/model_2.tex')

post2 <- extract.samples(m2)

png('presentation2/figures/model_2.png',
    height = 4, width = 7, units = 'in', res = 300)

par(mfrow = c(1, 2))

# plot(graphLayout(dag2))

# plot(r ~ g,
#      xlab = 'gene (g)',
#      ylab = 'response (r)',
#      main = 'r = rnorm(N) - .5*g')
# 
# plot(r ~ g, col = 'white',
#      xlab = 'gene (g)',
#      ylab = 'response (r)',
#      main = 'Post. Pred. Dist.')
# 
# gseq <- seq(from = min(g), to = max(g), len = 10)
# set.seed(123)
# mu <- link(m2, data = list(g = gseq, gbar = mean(g)))
# lines(gseq, apply(mu, 2, mean))
# shade(apply(mu, 2, PI, prob = .99), gseq)
# 
# set.seed(123)
# rsim <- sim(m2, data = list(g = gseq, gbar = mean(g)))
# shade(apply(rsim, 2, PI, prob = .99), gseq)

boxplot(r ~ G,
        xlab = 'Group (G)',
        ylab = 'response (r)',
        main = 'g = g+2*G; r = r-2*G',
        lwd = 2)

G1 <- rnorm(10000, post2$alpha[, 1], post2$sigma)
G2 <- rnorm(10000, post2$alpha[, 2], post2$sigma)
G3 <- rnorm(10000, post2$alpha[, 3], post2$sigma)

dens(G1,
     xlim = c(-11, -1),
     ylim = c(0, .5),
     lwd = 2,
     xlab = 'Effect of G on r',
     main = 'Post. Pred. Dist.')
dens(G2, lwd = 2, col = 2, add = TRUE)
dens(G3, lwd = 2, col = 4, add = TRUE)

text(mean(G1)+1, .45, 'G1')
text(mean(G2), col = 2, .45, 'G2')
text(mean(G3)-1, col = 4, .45, 'G3')

dev.off()

# model 3: include r on G ----
# dag
dag3 <- dagitty("dag{ g -> r g -> G G -> r }")

## simulate
N <- 100

set.seed(1234)
G <- sample(c(1,2,3), N, replace = TRUE)
table(G)
g <- rpois(N, 5) + 2 * G
boxplot(g ~ G)
G <- aggregate(g, list(G = G), mean)$x[G]

r <- rnorm(N) - .5 * g - G

## run
m3 <- ulam(
  flist = alist(
    # r as normal
    r ~ dnorm(mu, sigma),
    # regress r on g and G
    mu <- alpha + beta1 * (g - gbar) + beta2 * (G - Gbar),
    # declare priors
    alpha ~ dnorm(0, 1),
    beta1 ~ dnorm(0, 1),
    beta2 ~ dnorm(0, 1),
    sigma ~ dunif(0, 5)
  ),
  data = list(
    N = N,
    g = g,
    G = G,
    r = r,
    gbar = mean(g),
    Gbar = mean(G)
  ),
  chains = 4,
  cores = 4,
  seed = 1234
)

## visualize
png('presentation2/figures/trank_3.png')
trankplot(m3)
dev.off()

precis(m3)

print(xtable(precis(m3)),
      booktabs = TRUE,
      file = 'presentation2/tables/model_3.tex')

post3 <- extract.samples(m3)


png('presentation2/figures/model_3.png',
    height = 4, width = 7, units = 'in', res = 300)

par(mfrow = c(1, 2))

# plot(graphLayout(dag3))

plot(r ~ G, lwd = 2,
     xlab = 'Group (G)',
     ylab = 'response (r)',
     main = 'r = r - G/mean(g)')

plot(r ~ G, col = 'white', lwd = 2,
     xlab = 'Group (G)',
     ylab = 'response (r)',
     main = 'Post. Pred. Dist.')

Gseq <- seq(from = min(G), to = max(G), len = 5)
set.seed(123)
mu <- link(m3, data = list(G = Gseq, Gbar = mean(G), gbar = mean(g)))

lines(Gseq, apply(mu, 2, mean))
shade(apply(mu, 2, PI, prob = .99), Gseq)

set.seed(123)
rsim <- sim(m3, data = list(G = Gseq, Gbar = mean(G), gbar = mean(g)))
shade(apply(rsim, 2, PI, prob = .99), Gseq)

dev.off()

# model 4: index on D (multi-level) ----
# dag
dag4 <- dagitty("dag{ g -> r g -> G G -> r D -> r}")

## simulate
N <- 100

set.seed(1234)
GI <- sample(c(1,2,3), N, replace = TRUE)
g <- rpois(N, 5) + 2 * GI
G <- aggregate(g, list(G = GI), mean)$x[GI]
r <- rnorm(N) - .5 * g - G

# replicated for 3 Ds
g <- rep(g, 3)
GI <- rep(GI, 3)
G <- rep(G, 3)
r <- rep(r, 3)
D <- rep(c(1,2,3), each = N)

# introduce the effect of D
r <- ifelse(D == GI, r * - .5 * D, r)

boxplot(r ~ paste(D, GI))

## run
m4 <- ulam(
  flist = alist(
    # r as normal
    r ~ dnorm(mu, sigma),
    # regress r on g and G
    mu <- alpha + beta1[D] * (g - gbar) + beta2[D] * (G - Gbar),
    # declare priors
    alpha ~ dnorm(0, 1),
    beta1[D] ~ dnorm(0, 1),
    beta2[D] ~ dnorm(0, 1),
    sigma ~ dunif(0, 10)
  ),
  data = list(
    N = N,
    g = g,
    G = G,
    r = r,
    D = D,
    gbar = mean(g),
    Gbar = mean(G)
  ),
  chains = 4,
  cores = 4,
  seed = 1234
)

## visualize
png('presentation2/figures/trank_4.png')
trankplot(m4)
dev.off()

precis(m4, depth = 2)

print(xtable(precis(m4, depth = 2)),
      booktabs = TRUE,
      file = 'presentation2/tables/model_4.tex')

post4 <- extract.samples(m4)


png('presentation2/figures/model_4.png',
    height = 4, width = 7, units = 'in', res = 300)

par(mfrow = c(1, 2))

# plot(graphLayout(dag4))

boxplot(r ~ paste(D, GI, sep = '/'), lwd = 2,
        xlab = 'Drug effect (D), by Group (G)',
        ylab = 'response (r)',
        main = 'ifelse(D == G,r*-.5*D,r)')

set.seed(1234)
D1G1 <- link(m4, data = list(D = 1, G = 1, gbar = mean(g), Gbar = mean(g)))
D1G2 <- link(m4, data = list(D = 1, G = 2, gbar = mean(g), Gbar = mean(g)))
D1G3 <- link(m4, data = list(D = 1, G = 3, gbar = mean(g), Gbar = mean(g)))

dens(D1G1, lwd = 2,
     ylim = c(0, .1),
     main = 'Post. Pred. Dist.',
     xlab = 'Effect D1 on G')
dens(D1G2, lwd = 2, col = 2, add = TRUE)
dens(D1G3, lwd = 2, col = 4, add = TRUE)

text(mean(D1G1), .08, 'G1')
text(mean(D1G2), .09, 'G2', col = 2)
text(mean(D1G3), .1, 'G3', col = 4)

dev.off()

# model 5: index on C (multi-multi-level) ----
# dag
dag5 <- dagitty("dag{ g -> r g -> G G -> r D -> r C -> r}")

## simulate
N <- 100

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
