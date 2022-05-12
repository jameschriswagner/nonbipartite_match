# create a vector stating whether an observation is higher or lower than
# its paired counterpart
# 0 for lower, 1 for higher, 2 for unpaired
highlow <- function(og_df, pairs){
  hl = rep(2, nrow(og_df))
  for(i in 1:nrow(pairs)){
    z1 = og_df$Z[pairs$Group1.Row[i]]
    z2 = og_df$Z[pairs$Group2.Row[i]]
    
    if(z1 >= z2){
      hl[pairs$Group1.Row[i]] = 1
      hl[pairs$Group2.Row[i]] = 0
    } else {
      hl[pairs$Group1.Row[i]] = 0
      hl[pairs$Group2.Row[i]] = 1
    }
  }
  return(hl)
}

stupid_highlow = rep(0, nrow(obs))
med = median(obs$Z)
for(i in 1:nrow(obs)){
  if(obs$Z[i] >= med){
    stupid_highlow[i] = 1
  } else {
    stupid_highlow[i] = 0
  }
}
obs$stupid_hl = stupid_highlow
obs$V = simdata$V


# summary statistics
obs$highlow.mahalanobis = highlow(obs, mahalanobis.match)
obs$highlow.propensity = highlow(obs, propensity.match)
obs$highlow.caliper = highlow(obs, caliper.match)
obs$highlow.nf1 = highlow(obs, nearfar1.match)
obs$highlow.nf2 = highlow(obs, nearfar2.match)




round(sapply(high, mean), 2);round(sapply(high, sd), 2)
round(sapply(low, mean), 2);round(sapply(low, sd), 2)

sum(1 * (obs$highlow.caliper == 2))\



# hypothesis testing
Y1 = obs$Y[obs$highlow.nf2 == 1]
Y2 = obs$Y[obs$highlow.nf2 == 0]
wilcox.test(Y1, Y2, paired=F)








