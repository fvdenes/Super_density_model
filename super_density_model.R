# -----------------------------------------------------------------------------------
# Code to simulate data and test superdensity model for count data using QPAD offsets
# 
# 
# -----------------------------------------------------------------------------------

library(R2jags)

# Simulate overall density (superdensity) at each station
R <- 100    # The number of sites
T <- 3      # The number of surveys at each site

# Simulate the vegetation and the 
veg <- sort(runif(R, -1.5, 1.5))    # The standardized vegetation at each site
veg2 <- veg^2
alpha <- 2
beta1 <- 2
beta2 <- -2
lam <- exp(alpha + beta1*veg + beta2*(veg^2))
N <- rpois(R, lam)
sup.dens <- N/3.14
hist(sup.dens, main = "", xlab = "Superdensity")
plot(veg, sup.dens, xlab = "Standardized vegetation index", ylab = "Superdensity", frame.plot = FALSE)
points(veg, (lam/3.14), type="l")

# Simulate the offsets
exp_off <- matrix(NA, R, T)
for(i in 1:R){
  for(j in 1:T){
    exp_off[i, j] <- rgamma(1, 2, 2)
  }
}

off <- log(exp_off)

# Simulate the precentage of superdensity available at each survey
psisurv <- matrix(NA, R, T)
for(i in 1:R){
  psisurv[i, ] <- runif(T, 0.6, 1)
}

# Simulate the actual density during each survey
delta <- sup.dens*psisurv

# Simulate the lambda for count data at each site
lambda <- delta*exp_off

# Simulate the counts at each site
y <- matrix(NA, R, T)
for(i in 1:R){
  for(j in 1:T){
    y[i, j] <- rpois(1, lambda[i, j])
  }
}
hist(y)
plot(veg, y[,1])


mod <- glm(y[, 1] ~ veg + I(veg^2), family = poisson("log"), offset = off[, 1])
mod$coefficients

plot(veg, sup.dens)
points(veg, exp(predict(mod.naive, data = veg)), type="l")
pred.mod <- exp(mod$coefficients[1] + mod$coefficients[2]*veg + mod$coefficients[3]*(veg^2))
points(veg, pred.mod, type="l", lwd = 2, col="red")
points(veg, lam/3.14, type="l", lty = 2, lwd = 2)



# Build the hierarchical model incorporating changes in density throughout the season
# -----------------------------------------------------------------------------------

# First, run a basic model to compare to the mle model with the offset

sink("qpad0.txt")
cat("

model{
  
  alpha ~ dnorm(0, 0.001)I(-5, 5)
  beta1.d ~ dnorm(0, 0.001)I(-5, 5)
  beta2.d ~ dnorm(0, 0.001)
  
  for(i in 1:R){
    sup.dens[i] <- exp(alpha + beta1.d*veg[i] + beta2.d*veg2[i])
    log(lambda[i]) <- log(sup.dens[i]) + off[i]
    y[i] ~ dpois(lambda[i])
  }
}


", fill=TRUE)
sink()


data <- list(y=y[, 1], off=off[,1], R = R, veg=veg, veg2=veg2)

params <- c("alpha", "beta1.d", "beta2.d", "sup.dens")

dinit <- y[, 1]/off[, 1]

inits <- function(){list(alpha=rnorm(1, 0, 1), beta1.d=rnorm(1, 0, 1), beta2.d=rnorm(1, 0, 1))}

ni<-100
nc<-3
nb<-10
nt<-1

mod0 <- jags(data, inits, params, "qpad0.txt", n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt)

mod0$BUGSoutput$mean[1:3]
mod$coefficients


plot(sup.dens, mod0$BUGSoutput$mean$sup.dens)
plot(sup.dens, pred.mod)



sink("qpad1.txt")
cat("
    
    model{
    
    alpha ~ dnorm(0, 0.001)I(-5, 5)
    beta1.d ~ dnorm(0, 0.001)I(-5, 5)
    beta2.d ~ dnorm(0, 0.001)
    
    
    for(i in 1:R){
    sup.dens[i] <- exp(alpha + beta1.d*veg[i] + beta2.d*veg2[i])
    for(j in 1:T){
      psi[i, j] ~ dunif(0, 1)
      log(delta[i, j]) <- log(sup.dens[i]) + log(psi[i, j])
    
      log(lambda[i, j]) <- log(delta[i, j]) + off[i, j]
      y[i, j] ~ dpois(lambda[i, j])
      }
    }
    
  }
    
    
", fill=TRUE)
sink()


data <- list(y=y, off=off, R = R, T=T, veg=veg, veg2=veg2)

mod1 <- jags(data, inits, params, "qpad1.txt", n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt)
print(mod1, dig=2)

par(mfrow=c(1, 2))
plot(sup.dens, mod0$BUGSoutput$mean$sup.dens)
plot(sup.dens, mod1$BUGSoutput$mean$sup.dens)


par(mfrow=c(1, 1))
plot(veg, sup.dens)
points(veg, mod0$BUGSoutput$mean$sup.dens, type="l")
points(veg, mod1$BUGSoutput$mean$sup.dens, type="l", col="red", lwd=2)
points(veg, lam/3.14, type="l", lty = 2, lwd = 2)


K <- 7    # The number of years of surveys
yr.effect <- rnorm(K, 0, 0.5)

lam.yr<-matrix(NA, R, K)
for(i in 1:R){
  for(k in 1:K){
    lam.yr[i, k] <- exp(alpha+yr.effect[k] + beta1*veg[i] + beta2*veg2[i])
  }
}
lam.yr

N.yr <- matrix(NA, R, K)
for(i in 1:R){
  N.yr[i, ] <- rpois(K, lam.yr[i, ])
}

sup.dens.yr <- N.yr/3.14


# Simulate the offsets
exp_off.yr <- array(NA, dim = c(R, K, T))
for(i in 1:R){
  for(k in 1:K){
    for(j in 1:T){
      exp_off.yr[i, k, j] <- rgamma(1, 2, 2)
    }
  }
}

off.yr <- log(exp_off.yr)


# Simulate the precentage of superdensity available at each survey in each year
psisurv.yr <- array(NA, dim = c(R, K, T))
for(i in 1:R){
  for(k in 1:K){
    psisurv.yr[i, k, ] <- runif(T, 0.6, 1)
  }
}

# Simulate the actual density during each survey
delta.yr <- array(NA, dim = c(R, K, T))
for(k in 1:K){
  delta.yr[, k, ] <- sup.dens.yr[, k]*psisurv.yr[, k, ]
}

# Simulate the lambda for count data at each site
lambda.yr <- array(NA, dim = c(R, K, T))
for(k in 1:K){
  lambda.yr[, k, ] <- delta.yr[, k, ]*exp_off.yr[, k, ]
}


# Simulate the counts at each site
y.yr <- array(NA, dim = c(R, K, T))
for(i in 1:R){
  for(k in 1:K){
    for(j in 1:T){
      y.yr[i, k, j] <- rpois(1, lambda.yr[i, k, j])
    }
  }
}

hist(y.yr[,1,])
plot(veg, y.yr[,1,1])


sink("qpad2.txt")
cat("
    
    model{
    for(k in 1:K){
      alpha[k] ~ dnorm(0, 0.001)I(-5, 5)
      }
    
    beta1.d ~ dnorm(0, 0.001)I(-5, 5)
    beta2.d ~ dnorm(0, 0.001)
    
    
    for(k in 1:K){
      for(i in 1:R){
        sup.dens[i, k] <- exp(alpha[k] + beta1.d*veg[i] + beta2.d*veg2[i])
        for(j in 1:T){
          psi[i, k, j] ~ dunif(0, 1)
          log(delta[i, k, j]) <- log(sup.dens[i, k]) + log(psi[i, k, j])
    
          log(lambda[i, k, j]) <- log(delta[i, k, j]) + off[i, k, j]
          y[i, k, j] ~ dpois(lambda[i, k, j])
        }
      }
    # Derived quantities
    meanD[k] <- mean(sup.dens[, k])
    }
    
  }
    
    
", fill=TRUE)
sink()


data <- list(y=y.yr, off=off.yr, R = R, T=T, K = K, veg=veg, veg2=veg2)

params <- c("alpha", "beta1.d", "beta2.d", "meanD")

inits <- function(){list(alpha=rnorm(K, 0, 1), beta1.d=rnorm(1, 0, 1), beta2.d=rnorm(1, 0, 1))}

mod2 <- jags(data, inits, params, "qpad2.txt", n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt)
print(mod2, dig=2)

plot(apply(sup.dens.yr, 2, mean), mod2$BUGSoutput$mean$meanD)

par(mfrow=c(1, 2))

plot(apply(sup.dens.yr, 2, mean), mod2$BUGSoutput$mean$meanD, xlab = "True mean density", ylab = "Estimated mean density", pch = 16, bty = "l")
