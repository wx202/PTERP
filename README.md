# PTERP
PTE (proportion of treatment explained) and RP (relative power) estimates to evaluate the surrogacy of a surrogate for the primary outcome 

install.packages("PTERP")

library(PTERP)

data=PTERP::data

output=PTERP(data,ncut=c(50,100,150,200,500,1000),n.resam=5)

Wang, X., Parast, L., Tian, L. and Cai, T., 2022. Towards Optimal Use of Surrogate Markers to Improve Power. arXiv preprint arXiv:2209.08414.
