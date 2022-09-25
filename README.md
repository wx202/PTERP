# PTERP
PTE (proportion of treatment explained) and RP (relative power) estimates for optimally-transformed surrogate 

install.packages("PTERP")

library(PTERP)

data=PTERP::data

output=PTERP(data,ncut=c(50,100,150,200,500,1000),n.resam=5)
