rm(list = ls())

library(coda)

logit <- function(theta, a, b) {
  log((theta-a)/(b-theta))
}
logit.inv <- function(z, a, b) {
  b-(b-a)/(1+exp(z))
}



# Read in Data ------------------------------------------------------------

n <- 1566
n.spline.bf <- 8
n.phi <- 1
n.lambda <- 6
n.theta <- n.spline.bf + n.phi + n.lambda + 1
spline.seq <- read.table("../../data/splineSequence")[1:n, 1]
alpha <- read.table("../../data/alpha")[1:n, 1]
beta <- read.table("../../data/beta")[1:n, 1]

beta.s <- t(matrix(scan("beta-theta-1"), nrow = n.theta, byrow = T))
dim(beta.s)
beta.fitted <- matrix(scan("beta-fitted-1"), nrow = n, byrow = T)
beta.w <- matrix(scan("beta-w-1"), nrow = n, byrow = T)


colnames(beta.s) <- c(paste0("spline.", 1:n.spline.bf), "phi", "lambda.11",
                       "lambda.21", "lambda.31", "lambda.22", "lambda.32", 
                       "lambda.33", "rho")

n.samples <- ncol(beta.fitted)
burn.in <- floor(.25*n.samples)
sub <- burn.in:(n.samples - 1)

plot(mcmc(beta.s), density=F)

vals <- summary(window(mcmc(beta.s), start=burn.in))
vals$quantiles

y.hat <- apply(beta.fitted[,sub], 1, mean)

par(mfrow = c(1,1))
plot(beta, y.hat, pch = 19, cex = 0.5, xlab = "True Values", 
     ylab = "Fitted Values", cex.lab = 1.5, xlim = c(0, 1), ylim = c(0, 1))
lines(beta, beta)


# Obtain credible intervals for fitted values
fitted.quants <- summary(window(mcmc(t(beta.fitted)), start = burn.in))$quantiles
# Output model fit with credible intervals --------------------------------
png("../../../figures/fit-m2-stage2.png", width = 480, height = 480)
par(mar=c(6.5,4.1,4.1,2.1))
plot(beta, fitted.quants[, 3], pch = 19, cex = 0.5, ylab = "Posterior Medians",
     xlab = "True Values", cex.lab = 1.5, bty = 'n', las = 1)
segments(beta, fitted.quants[, 1], beta, fitted.quants[, 5], col = 'grey')
points(beta, fitted.quants[, 3], pch = 19, cex = 0.5)
text(0.18, -0.1, "Model 2", xpd = TRUE, cex = 1.5)
dev.off()


# New curve ---------------------------------------------------------------
alpha.vals <- read.table("../../data/alpha-spline-curve")[, 1]
graph.n <- length(alpha.vals)
graph.n.iter <- 1875
mk <- 87
alpha.graph.fitted <- matrix(scan("curve/alpha-figure-vals"), nrow = graph.n.iter, byrow = F)

quant <- function(x){quantile(x, prob=c(0.5, 0.025, 0.975))}
orig.graph.fitted <- apply(alpha.graph.fitted, 2, quant)

graph.fitted <- matrix(0, nrow = 3, ncol = graph.n)
for (i in 1:graph.n) {
  for (j in 1:nrow(graph.fitted)) {
    graph.fitted[j, i] <- mean(orig.graph.fitted[j, ((i-1)*mk + 1):(i*mk)])
  }
}

png("../../../figures/spline-fit-m2-stage2.png", width = 480, height = 480)
par(mar=c(6.5,4.1,4.1,2.1))
plot(alpha.vals, graph.fitted[1,], axes=F, type="l", ylim = range(beta),
     xlab="Anthropogenic Noise", ylab="Biological Sounds", bty = 'n', cex.lab = 1.5)
polygon(c(alpha.vals,rev(alpha.vals)), c(graph.fitted[2,],rev(graph.fitted[3,])), col =
          "lightsteelblue2", border = FALSE)
lines(alpha.vals, graph.fitted[1,])
points(alpha, beta, pch = 19, cex = 0.5)
axis(1, las = 1)
axis(2, las = 1)
text(0.5, -0.1, "Model 2", xpd = TRUE, cex = 1.5)
dev.off()

# Credible Interval Model Comparison Computation --------------------------

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
beta.s.1 <- t(matrix(scan("beta-theta-1"), nrow=n.theta, byrow=T))
beta.s.2 <- t(matrix(scan("beta-theta-2"), nrow=n.theta, byrow=T))
beta.s.3 <- t(matrix(scan("beta-theta-3"), nrow=n.theta, byrow=T))

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

# Output for Beta Model ---------------------------------------------------
# write.table(y.hat, "beta-fitted-mean", row.names = FALSE,
#             col.names = FALSE, sep = "\t")

# Output for Prediction ---------------------------------------------------
write.table(t(beta.s)[, sub], "beta-theta-post-burn", 
            row.names = FALSE, col.names = FALSE, sep = "\t")
