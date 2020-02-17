rm(list=ls())

library(coda)

logit <- function(theta, a, b){log((theta-a)/(b-theta))}
logit.inv <- function(z, a, b){b-(b-a)/(1+exp(z))}


# Read in data ------------------------------------------------------------
n <- 1566
alpha <- read.table("../../data/alpha")[1:n,1]
beta <- read.table("../../data/beta")[1:n,1]
beta.fitted <- matrix(scan("beta-fitted-1"), nrow=n, byrow=T)
beta.w <- matrix(scan("beta-w-1"), nrow=n, byrow=T)
n.spline.bf <- 8
n.samples <- ncol(beta.fitted)
n.phi <- 1
n.params <- n.spline.bf + n.phi + 1 + 1
beta.s <- t(matrix(scan("beta-theta-1"), nrow=n.params, byrow=T))
dim(beta.s)
colnames(beta.s) <- c(paste0("spline.",1:n.spline.bf), paste0("phi.", 1:n.phi),
                      "sigma.sq", "rho")



n.samples <- ncol(beta.fitted)
burn.in <- floor(0.5*n.samples)
sub <- (burn.in+1):n.samples

plot(window(mcmc(beta.s), start = burn.in), density = FALSE)
plot(mcmc(beta.s), density = FALSE)
summary(window(mcmc(beta.s), start=burn.in))

y.hat <- apply(beta.fitted[,sub], 1, mean)
w.hat <- apply(beta.w[, sub], 1, mean)

# Assess model fitS
plot(beta, y.hat, pch = 19, cex = 0.5, ylab = "Fitted Values", xlab = "True values",
     cex.lab = 1.5, ylim = c(0, 1), xlim = c(0, 1))
lines(beta, beta, lty = 2, col = 'red')


# Obtain credible intervals for fitted values
fitted.quants <- summary(window(mcmc(t(beta.fitted)), start = burn.in))$quantiles
# Output model fit with credible intervals --------------------------------
png("../../../figures/fit-m1-stage2.png", width = 480, height = 480)
par(mar=c(6.5,4.1,4.1,2.1))
plot(beta, fitted.quants[, 3], pch = 19, cex = 0.5, ylab = "Posterior Medians",
     xlab = "True Values", cex.lab = 1.5, bty = 'n', las = 1)
segments(beta, fitted.quants[, 1], beta, fitted.quants[, 5], col = 'grey')
points(beta, fitted.quants[, 3], pch = 19, cex = 0.5)
text(0.18, -0.1, "Model 1", xpd = TRUE, cex = 1.5)
dev.off()


# Spline Curve ------------------------------------------------------------
alpha.vals <- read.table("../../data/alpha-spline-curve")[, 1]
graph.n.iter <- 2500
graph.n <- length(alpha.vals)
mk <- 29
alpha.graph.fitted <- matrix(scan("curve/alpha-figure-vals"), nrow = graph.n.iter, byrow = F)

quant <- function(x){quantile(x, prob=c(0.5, 0.025, 0.975))}
orig.graph.fitted <- apply(alpha.graph.fitted, 2, quant)

graph.fitted <- matrix(0, nrow = 3, ncol = graph.n)
for (i in 1:graph.n) {
  for (j in 1:nrow(graph.fitted)) {
    graph.fitted[j, i] <- mean(orig.graph.fitted[j, ((i-1)*mk + 1):(i*mk)])
  }
}

png("../../../figures/spline-fit-m1-stage2.png", width = 480, height = 480)
par(mar=c(6.5,4.1,4.1,2.1))
plot(alpha.vals, graph.fitted[1,], axes=F, type="l", ylim=range(beta),
     xlab="Anthropogenic Noise", ylab="Biological Sounds", cex.lab = 1.5, bty = 'n')
polygon(c(alpha.vals,rev(alpha.vals)), c(graph.fitted[2,],rev(graph.fitted[3,])), col =
          "lightsteelblue2", border = FALSE)
lines(alpha.vals, graph.fitted[1,])
points(alpha, beta, pch = 19, cex = 0.5)
axis(1, las = 1)
axis(2, las = 2)
text(0.5, -0.1, "Model 1", xpd = TRUE, cex = 1.5)
dev.off()

counter <- 0
for (i in 1:n) {
  col <- which.min(abs(alpha[i] - alpha.vals))
  if (beta[i] > graph.fitted[2, col] & beta[i] < graph.fitted[3, col]) {
    counter <- counter + 1
  }
}
coverage.95 <- counter / n * 100
coverage.95

# Assess Convergence ------------------------------------------------------
beta.s.1 <- t(matrix(scan("beta-theta-1"), nrow=n.params, byrow=T))
beta.s.2 <- t(matrix(scan("beta-theta-2"), nrow=n.params, byrow=T))
beta.s.3 <- t(matrix(scan("beta-theta-3"), nrow=n.params, byrow=T))

beta.s.1.post.burn <- window(beta.s.1, start = burn.in)
beta.s.2.post.burn <- window(beta.s.2, start = burn.in)
beta.s.3.post.burn <- window(beta.s.3, start = burn.in)

# Join together the three chains for all parameters
beta.s.all <- rbind(beta.s.1.post.burn, beta.s.2.post.burn, beta.s.3.post.burn)

sum.chain.all <- summary(mcmc(beta.s.all))

# Calculate empirical 95% CI intervals
sum.chain.1 <- summary(window(mcmc(beta.s.1), start = burn.in))
sum.chain.2 <- summary(window(mcmc(beta.s.2), start = burn.in))
sum.chain.3 <- summary(window(mcmc(beta.s.3), start = burn.in))
within.1 <- sum.chain.1$quantiles[, 5] - sum.chain.1$quantiles[, 1]
within.2 <- sum.chain.2$quantiles[, 5] - sum.chain.2$quantiles[, 1]
within.3 <- sum.chain.3$quantiles[, 5] - sum.chain.3$quantiles[, 1]

within.vals <- data.frame(within.1, within.2, within.3)

total.length <- sum.chain.all$quantiles[, 5] - sum.chain.all$quantiles[, 1]
mean.within.vals <- apply(within.vals, 1, mean)
r.hat.interval <- total.length / mean.within.vals

r.hat.interval

chains <- mcmc.list(mcmc(beta.s.1), mcmc(beta.s.2),
                    mcmc(beta.s.3))
plot(chains, density = FALSE)

# For predictions ---------------------------------------------------------

write.table(t(beta.s)[, sub], "beta-theta-post-burn", row.names = FALSE,
                        col.names = FALSE, sep = "\t")


