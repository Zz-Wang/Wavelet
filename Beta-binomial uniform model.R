####################
# Beta-Binomial-Uniform model
# from Don Berry's HM (FDA "white paper", Sec 3.1)
####################

model{
  for( i in 1:I) {
    x[i] ~ dbin(p[i] , n[i])
    p[i] ~ dbeta(a,b)
  }
  
  a ~ dunif(0,10)
  b ~ dunif(0,10)
  #   a ~ dgamma(2,2)
  #   b ~ dgamma(2,2)
  
}  #  end of BUGS code


#Data:
list(x = c(20, 4, 11, 10, 5, 36, 9, 7, 4, NA),
     n = c(20, 10, 16, 19, 14, 46, 10, 9, 6, 1), I=10)


# Inits:
list(a=4, b=2, p = c(.5, .5, .5, .5, .5, .5, .5, .5, .5, .5))
