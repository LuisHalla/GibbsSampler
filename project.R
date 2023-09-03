#### LIBRARIES ####
library(LaplacesDemon)

#### FUNCTIONS ####

# Mixture
#' @param n number of samples
#' @param k number of mixtures
#' @param p mixture proportions --> len=k
#' @param mu mixture means --> len=k
#' @param sig mixture sds --> len=k
rmix <- function(n, k, p, mu, sig){
  
  j <- sample(1:k, size=n, prob=p, replace=T) # j = 1:k
  
  return(rnorm(n, mu[j], sig[j]))
}

# Normalize fc_c
normalize <- function(x) x/sum(x)

# Full conditional c
#' @param n number of samples
#' @param k number of mixtures
#' @param p mixture proportions --> len=k
#' @param mu mixture means --> len=k
#' @param sig mixture sds --> len=k
#' @param Y mixture --> len=n
fc_c <- function(n, k, p, mu, sig, Y){

  dev <- outer(mu, Y, "-") # kxn matrix whose rows are deviations of Y_j (len=n) from mu_j
  fc <- matrix(nrow=k, ncol=n)

  for(j in 1:k){

    normal <- sqrt(sig[j]/(2*pi)) * exp( (-1/2) * dev[j,]^2 * sig[j] )
    fc[j,] <- p[j] * sqrt(sig[j]) * normal
  }

  fc <- apply(fc, 2, normalize) # normalize columns so they are probabilities

  c <- numeric(n)
  for(i in 1:n) c[i] <- sample(1:k, size=1, prob=fc[,i], replace=T)

  return(c)
}


# Full conditional p
#' @param alpha prior shape Dirichlet
#' @param nj counts of class allocations (1,...,k) --> len=k
fc_p <- function(alpha, nj) rdirichlet(1, alpha+nj)


# Full conditional mu
#' @param k number of mixtures
#' @param m prior mean Normal
#' @param v prior sd Normal
#' @param c class allocations --> len=n
#' @param nj counts of class allocations (1,...,k) --> len=k
#' @param sig mixture sds --> len=k
#' @param Y mixture --> len=n
fc_mu <- function(k, m, v, c, nj, sig, Y){
  
  mu <- numeric(k)

  for(j in 1:k){
    
    zj <- nj[j] * sig[j] * v^2
    zj <- zj / (zj + 1)
    yj <- ifelse(nj[j]==0, 0, mean(Y[c==j]))
    
    mu_post <- zj*yj + (1 - zj)*m
    sig_post <- v^2 / (nj[j]*sig[j] + 1)
    
    mu[j] <- rnorm(1, mean=mu_post, sd=sqrt(sig_post))
  }
  
  return(mu)
}

# Full conditional sig
#' @param k number of mixtures
#' @param g prior shape Gamma
#' @param l prior rate Gamma
#' @param c class allocations --> len=n
#' @param nj counts of class allocations (1,...,k) --> len=k
#' @param mu mixture means --> len=k
#' @param Y mixture --> len=n
fc_sig <- function(k, g, l, c, nj, mu, Y){
  
  sig <- numeric(k)
  
  for(j in 1:k){
    
    shape_post <- g + nj[j]/2
    rate_post <- l + sum((Y[c==j] - mu[j])^2) / 2
      
    sig[j] <- rgamma(1, shape=shape_post, rate=rate_post)
  }
  
  return(sig)
}


# Gibbs Sampler
GS <- function(nsim, Y, n, k, alpha, m, v, g, l){
  
  # List of results of GS
  GS <- list(c=matrix(nrow=nsim, ncol=n),
             p=matrix(nrow=nsim, ncol=k),
             mu=matrix(nrow=nsim, ncol=k),
             sig=matrix(nrow=nsim, ncol=k))
  
  # Gen p1, mu1, sig1 using conjugate priors
  GS$p[1,] <- rdirichlet(1, alpha)
  GS$mu[1,] <- rnorm(k, m, v)
  GS$sig[1,] <- rgamma(k, g, rate=l)
  
  # Gen c1 using p1, mu1, sig1
  GS$c[1,] <- fc_c(n, k, GS$p[1,], GS$mu[1,], GS$sig[1,], Y)

  for(i in 2:nsim){ # Step
    
    # Latent variable c allows the computation of class allocations nj
    nj <- colSums(outer(GS$c[i-1,], 1:k, FUN="=="))
    
    # Gen p with nj using full conditionals
    GS$p[i,] <- fc_p(alpha, nj)
    # Gen mu with nj using full conditionals
    GS$mu[i,] <- fc_mu(k, m, v, GS$c[i-1,], nj, GS$sig[i-1,], Y)
    # Gen sig with nj, mu using full conditionals
    GS$sig[i,] <- fc_sig(k, g, l, GS$c[i-1,], nj, GS$mu[i,], Y)
    
    # Gen c with p, mu, sigma using full conditional
    GS$c[i,] <- fc_c(n, k, GS$p[i,], GS$mu[i,], GS$sig[i,], Y)
    # Latent variable --> independence
  }
  
  return(GS)
}




#### MAIN ####

# Globals
nsim <- 1e3 # n steps
n <- 1e3 # n points in mixture
k <- 3 # n mixtures

# Hyperparameters
alpha <- rep(1, k) # shape Dirichlet
m <- 0 # mean Normal
v <- 1 # sd Normal
g <- 1 # shape Gamma
l <- 1 # rate Gamma

# Mixture: fix p0, mu0, sig0 by user and gen mixture Y
p0 <- c(.1, .2, .7)
mu0 <- c(-4, 0, 4)
sig0 <- rep(1, k)

# p0 <- rdirichlet(1, alpha)
# mu0 <- rnorm(k, m, v)
# sig0 <- rgamma(k, g, rate=l)

Y <- rmix(n, k, p0, mu0, sig0)
hist(Y)

# Gibbs Sampler
GS_mix <- GS(nsim, Y, n, k, alpha, m, v, g, l)

# Plots
plot(GS_mix$mu[,1], ylim=c(min(mu0)-2,max(mu0)+2), type="l",
     xlab='nsim', ylab='mu',
     main=paste('Mixture means\n', n, 'observations'))
lines(GS_mix$mu[,2],col=2)
lines(GS_mix$mu[,3],col=3)

plot(GS_mix$p[,1], ylim=c(0,1), type="l",
     xlab='nsim', ylab='p',
     main=paste('Mixture probabilities\n', n, 'observations'))
lines(GS_mix$p[,2],col=2)
lines(GS_mix$p[,3],col=3)




#### APPLICATION TO FAITHFUL DATA ####
plot(faithful$eruptions, faithful$waiting)


# ERUPTIONS

# Globals
nsim <- 1e3 # n steps
n <- length(faithful$eruptions) # n points in mixture
k <- 2 # n mixtures

# Hyperparameters
alpha <- rep(1, k) # shape Dirichlet
m <- 0 # mean Normal
v <- 1 # sd Normal
g <- 1 # shape Gamma
l <- 1 # rate Gamma

# Gibbs Sampler eruptions
hist(faithful$eruptions)
GS_e <- GS(nsim, Y=faithful$eruptions, n, k, alpha, m, v, g, l)

# Plots
ymin <- min(GS_e$mu) - 1
ymax <- max(GS_e$mu) + 1

plot(GS_e$mu[,1], ylim=c(ymin, ymax), type="l",
     xlab='nsim', ylab='mu',
     main=paste('eruptions means\n', n, 'observations'))
lines(GS_e$mu[,2], col=2)

plot(GS_e$p[,1], ylim=c(0,1), type="l",
     xlab='nsim', ylab='p',
     main=paste('eruptions probabilities\n', n, 'observations'))
lines(GS_e$p[,2], col=2)


# WAITING

# Globals
nsim <- 1e3 # n steps
n <- length(faithful$waiting) # n points in mixture
k <- 2 # n mixtures

# Hyperparameters
alpha <- rep(1, k) # shape Dirichlet
m <- 70 # mean Normal
v <- 20 # sd Normal
g <- 1 # shape Gamma
l <- 1 # rate Gamma

# Gibbs Sampler waiting
hist(faithful$waiting)
GS_w <- GS(nsim, Y=faithful$waiting, n, k, alpha, m, v, g, l)

# Plots
ymin <- min(GS_w$mu) - 1
ymax <- max(GS_w$mu) + 1

plot(GS_w$mu[,1], ylim=c(ymin, ymax), type="l",
     xlab='nsim', ylab='mu',
     main=paste('waiting means\n', n, 'observations'))
lines(GS_w$mu[,2], col=2)

plot(GS_w$p[,1], ylim=c(0,1), type="l",
     xlab='nsim', ylab='p',
     main=paste('waiting probabilities\n', n, 'observations'))
lines(GS_w$p[,2], col=2)

