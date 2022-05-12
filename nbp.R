library(methods)
library(Matrix)
library(MASS)
library(Hmisc)
library(stats)
library(utils)
library(nbpMatching)

# getting some "observed" simulation data
obs = simdata[c(1:2, 4:8)]
covar = simdata[4:7]
n = nrow(covar)

# matching via Mahalanobis distance
mat1 = gendistance(covar, ndiscard = floor(n / 2))$dist
mahalanobis.match = subsetMatches(nonbimatch(distancematrix(mat1, nphantoms = floor(n / 2))),
                                  phantom = T,
                                  infinite = T)

# match using squared propensity distance
mat2 = matrix(data = rep(NA, n ^ 2),
             nrow = n, ncol = n,
             byrow = T)
for(k in 1:n){
  for(l in 1:n){
    if(obs$Z[k] == obs$Z[l]){
      mat2[k,l] = Inf
    } else {
      mat2[k,l] = (obs$eX[k] - obs$eX[l])^2
    }
  }
}
mat2 = make.phantoms(mat2, nphantoms = floor(n / 2))
propensity.match = subsetMatches(nonbimatch(distancematrix(mat2, nphantoms = floor(n / 2))),
                                 phantom = T,
                                 infinite = T)

# match using mahalanobis distance with propensity caliper of 0.5 SD
# infinite distance if same treatment
caliper = 0.5 * sd(obs$eX)
mat3 = mat1
for(i in 1:n){
  for(j in 1:n){
    if(obs$Z[i] == obs$Z[j]){
      mat3[i,j] = Inf
    } else if(abs(obs$eX[i] - obs$eX[j]) > caliper) {
      mat3[i,j] = Inf
    }
  }
}
mat3 = make.phantoms(mat3, nphantoms = floor(n / 4))
caliper.match = subsetMatches(nonbimatch(distancematrix(mat3, nphantoms = floor(n / 4))),
                              phantom = T,
                              infinite = T)
length(caliper.match$Group2.ID)

# near-far matching from Bo Lu, et. al using mahalanobis distance
# with propensity caliper
mat4 = mat1
nearfar <- function(a, b, delta){
  denom = (obs[a,2] - obs[b,2])^2
  return(delta / denom)
}
for(i in 1:n){
  for(j in 1:n){
    if(obs$Z[i] == obs$Z[j]){
      mat4[i,j] = Inf
    } else if(abs(obs$eX[i] - obs$eX[j]) > caliper) {
      mat4[i,j] = Inf
    } else {
      mat4[i,j] = nearfar(i, j, mat1[i,j])
    }
  }
}
mat4 = make.phantoms(mat4, nphantoms = floor(n / 4))
nearfar1.match = subsetMatches(nonbimatch(distancematrix(mat4, nphantoms = floor(n / 4))),
                               phantom = T,
                               infinite = T)
length(nearfar1.match$Group2.ID)

# near-far matching using different method
# with scalar penalty of 1.0 if Z-Z' < SD(Z)
# inspired by Baiocchi et al.
mat5 = mat1
for(i in 1:n){
  for(j in 1:n){
    if(obs$Z[i] == obs$Z[j]){
      mat5[i,j] = Inf
    } else if(abs(obs$Z[i] - obs$Z[j]) < sd(obs$Z)) {
      mat5[i,j] = mat5[i,j] + 1.0
    }
  }
}
mat5 = make.phantoms(mat5, nphantoms = floor(n / 4))
nearfar2.match = subsetMatches(nonbimatch(distancematrix(mat5, nphantoms = floor(n / 4))),
                               phantom = T,
                               infinite = T)

length(nearfar2.match$Group2.ID)
length(propensity.match$Group2.ID)





