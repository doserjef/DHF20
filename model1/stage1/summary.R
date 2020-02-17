rm(list = ls())

library(coda)

logit <- function(theta, a, b) {
  log((theta-a)/(b-theta))
}
logit.inv <- function(z, a, b) {
  b-(b-a)/(1+exp(z))
}



# Read in Data ------------------------------------------------------------
options("scipen"=100, "digits"=4)
n <- 1566
n.spline.bf <- 5
n.phi <- 1
n.theta <- n.spline.bf + n.phi + 1 + 1
x <- read.table("../../data/roadCovariate")[1:n, 1]
alpha <- read.table("../../data/alpha")[1:n, 1]

alpha.s <- t(matrix(scan("alpha-theta-1"), nrow = n.theta, byrow = T))
dim(alpha.s)
alpha.fitted <- matrix(scan("alpha-fitted-1"), nrow = n, byrow = T)
alpha.w <- matrix(scan("alpha-w-1"), nrow = n, byrow = T)
# Z <- t(matrix(scan("alpha-spline-1"), nrow = n.spline.bf, byrow = T))
Z <- as.matrix(read.table("alpha-spline"))

colnames(alpha.s) <- c(paste0("spline.", 1:n.spline.bf), paste0("phi.", 1:n.phi), 
                       "sigma.sq", "rho")

n.samples <- nrow(alpha.s)
burn.in <- floor(.1*n.samples)
sub <- burn.in:(n.samples - 1)

# plot(mcmc(alpha.s), density=F)

summary(window(mcmc(alpha.s), start=burn.in))

y.hat <- apply(alpha.fitted[,sub], 1, mean)

par(mfrow = c(1,1))
# plot(alpha, y.hat, pch = 19, cex = 0.5, xlab = "True Values", 
#      ylab = "Fitted Values", cex.lab = 1.5, xlim = c(0, 1), ylim = c(0, 1))
# lines(alpha, alpha, col = 'red')

# Obtain credible intervals for fitted values
fitted.quants <- summary(window(mcmc(t(alpha.fitted)), start = burn.in))$quantiles



# Obtain values for plotting ----------------------------------------------
rc.vals <- read.table("../../data/rc-curve")[, 1]
graph.n.iter <- 5000
graph.n <- length(rc.vals)
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

png("../../../figures/spline-fit-m1-stage1.png", width = 480, height = 480)
par(mar=c(6.5,4.1,4.1,2.1))
plot(rc.vals, graph.fitted[1,], axes=F, type="l", ylim=range(alpha),
     xlab="Road Covariate", ylab="Anthropogenic Noise", cex.lab = 1.5, bty = 'n')
polygon(c(rc.vals,rev(rc.vals)), c(graph.fitted[2,],rev(graph.fitted[3,])), col =
          "lightsteelblue2", border = FALSE)
lines(rc.vals, graph.fitted[1,])
points(x, alpha, pch = 19, cex = 0.5)
axis(1, las = 1)
axis(2, las = 1)
text(3.0, -0.27, "Model 1", xpd = TRUE, cex = 1.5)
dev.off()


# Output model fit with credible intervals --------------------------------
curr <- par()$mar
png("../../../figures/fit-m1-stage1.png", width = 480, height = 480)
par(mar=c(6.5,4.1,4.1,2.1))
plot(alpha, fitted.quants[, 3], pch = 19, cex = 0.5, ylab = "Posterior Medians", 
     xlab = "True Values", cex.lab = 1.5, bty = 'n', las = 1, 
     col = rgb(red = 0, green = 0, blue = 0, alpha = 0.3))
segments(alpha, fitted.quants[, 1], alpha, fitted.quants[, 5], col = 'grey')
points(alpha, fitted.quants[, 3], pch = 19, cex = 0.5, 
       col = rgb(red = 0, green = 0, blue = 0, alpha = 0.3))
text(0.5, -0.27, "Model 1", xpd = TRUE, cex = 1.5)
dev.off()


# Assess Convergence ------------------------------------------------------
alpha.s.1 <- t(matrix(scan("alpha-theta-1"), nrow=n.theta, byrow=T))
alpha.s.2 <- t(matrix(scan("alpha-theta-2"), nrow=n.theta, byrow=T))
alpha.s.3 <- t(matrix(scan("alpha-theta-3"), nrow=n.theta, byrow=T))

alpha.s.1.post.burn <- window(alpha.s.1, start = burn.in)
alpha.s.2.post.burn <- window(alpha.s.2, start = burn.in)
alpha.s.3.post.burn <- window(alpha.s.3, start = burn.in)

# Join together the three chains for all parameters
alpha.s.all <- rbind(alpha.s.1.post.burn, alpha.s.2.post.burn, alpha.s.3.post.burn)

sum.chain.all <- summary(mcmc(alpha.s.all))

# Calculate empirical 95% CI intervals
sum.chain.1 <- summary(window(mcmc(alpha.s.1), start = burn.in))
sum.chain.2 <- summary(window(mcmc(alpha.s.2), start = burn.in))
sum.chain.3 <- summary(window(mcmc(alpha.s.3), start = burn.in))
within.1 <- sum.chain.1$quantiles[, 5] - sum.chain.1$quantiles[, 1]
within.2 <- sum.chain.2$quantiles[, 5] - sum.chain.2$quantiles[, 1]
within.3 <- sum.chain.3$quantiles[, 5] - sum.chain.3$quantiles[, 1]

within.vals <- data.frame(within.1, within.2, within.3)

total.length <- sum.chain.all$quantiles[, 5] - sum.chain.all$quantiles[, 1]
mean.within.vals <- apply(within.vals, 1, mean)
r.hat.interval <- total.length1 / mean.within.vals

r.hat.interval

chains <- mcmc.list(mcmc(alpha.s.1), mcmc(alpha.s.2), 
                    mcmc(alpha.s.3))
plot(window(chains, start = burn.in), density = FALSE)
plot(chains, density = FALSE)

plot(mcmc(alpha.s.3), density = FALSE)


# Ensure fitted values are correct ----------------------------------------
mu.samples.pb <- alpha.mu[, sub]
plot(alpha, mu.samples.pb[, 500])

# Output for Beta Model ---------------------------------------------------
# write.table(y.hat, "alpha-fitted-mean", row.names = FALSE,
#             col.names = FALSE, sep = "\t")

write.table(alpha.fitted[, sub], "alpha-samples-post-burn",
            row.names = FALSE, col.names = FALSE, sep = "\t")

write.table(t(alpha.s)[, sub], "alpha-theta-post-burn", row.names = FALSE,
            col.names = FALSE, sep = "\t")

# For graph
write.table(t(alpha.s)[, 20001:25000], "alpha-theta-post-burn-graph", row.names = FALSE, 
            col.names = FALSE, sep = "\t")
