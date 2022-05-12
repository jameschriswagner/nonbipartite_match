library(Matrix)
library(MASS)
library(boot)


# make Y continuous, and dependent on the covariates
# have one normal confounding variable, with X1-3 centered around that
# f(X) = X^2 where X is the confounding variable
# Compare: G-computation generated with each method
#          

# find dataset for average temperature and GDP

# set parameters
set.seed(98)
N = 1000
simdata = data.frame(Y = rep(NA, N), Z = rep(NA, N), V = rep(NA, N),
                     X1 = rep(NA, N), X2 = rep(NA, N), X3 = rep(NA, N))

# simulate confounder according to a standard normal
simdata[3] = rnorm(N, 0, 1)

# simulate covariates according to MVN
Sig = matrix(data = c(1, 0.5, 0,
                      0.5, 1, -0.5,
                      0, -0.5, 1),
             nrow = 3, ncol = 3, byrow = T)
for(i in 1:N){
  v = simdata[i,3]
  simdata[i,4:6] = mvrnorm(n = N, mu = c(v, v, 0), Sigma = Sig)
}

# simulate continuous treatment
for(i in 1:N){
  mu = simdata[i,3]
  var = 1 + simdata[i,3] ^ 2
  simdata[i,2] = rnorm(1, mu, var)
}

# simulate response as a linear additive effect of treatment
errors = rnorm(N, 0, 1)
alpha = 0
beta1 = c(0.5, 0.5, 0)
beta2 = 0.2
simdata[1] = alpha + t(beta1 %*% t(simdata[4:6])) + beta2 * simdata$Z + errors




# estimating "propensity" score
# define generalized propensity score to be P(Z > mean(Z))
# estimate with a logistical regression model
d = mean(simdata$Z)
simdata$PZ <- 1 * (simdata$Z > d)
propensity <- glm(PZ ~ X1 + X2 + X3,
                  data = simdata,
                  family = "binomial")
simdata$eX <- inv.logit(predict(propensity, newdata = simdata))

