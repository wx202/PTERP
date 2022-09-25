########
# PTERE calculate PTE and RE to evaluate the surrogacy of a surrogate S for the primary outcome Y.
# Input: data (matrix): first column observed outcome yob, second column observed surrogate sob,
#                       third column treatment indicator aob;
#        ncut: sample sizes RE calculated at;
#        n.resam: number of perturbation resampling.
# Output: PTE estimate and RE estimates for different sample sizes,
#         and corresponding standard error estimates.
########
# install.packages("survival")
# library(survival)
PTERP=function(data,ncut=c(50,100,150,200,500,1000),n.resam=500){


  Kern.FUN <- function(zz,zi,bw)
  {
    out = (VTM(zz,length(zi))- zi)/bw
    dnorm(out)/bw

  }


  VTM<-function(vc, dm){
    matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
  }


  gopt=function(yob, sob, aob, n){
    if (min(range(sob[aob==1])==range(sob))>0){ supp1=range(sob[aob==1])
    }else if (min(sob[aob==1])-min(sob)==0){supp1=c(quantile(sob[aob==1],0),quantile(sob[aob==1],0.99))
    }else if (max(sob[aob==1])-max(sob)==0){ supp1=c(quantile(sob[aob==1],0.01),quantile(sob[aob==1],1))
    } else {print("You need to switch the groups.")}
    ind1=which.min(abs(supp1[1]-s)):which.min(abs(supp1[2]-s))
    supp0=range(sob[aob==0])
    ind0=which.min(abs(supp0[1]-s)):which.min(abs(supp0[2]-s))
    indc=intersect(ind1,ind0)
    ind1=setdiff(ind1,indc)
    ind0=setdiff(ind0,indc)

    bw = 1.06*sd(sob)*n^(-1/5)/(n^0.06)
    kern = Kern.FUN(zz=s,zi=sob,bw)
    f0.hat=apply(kern*(1-aob),2,mean)/mean((1-aob))
    f1.hat=apply(kern*aob,2,mean)/mean(aob)
    m1.s.hat=apply(yob*kern*aob,2,sum)/apply(kern*aob,2,sum)
    m0.s.hat=apply(yob*kern*(1-aob),2,sum)/apply(kern*(1-aob),2,sum)
    integrand<-(f0.hat^2/f1.hat)[indc]
    nc=length(indc)-1
    K2=(integrand[1] + integrand[nc+1] + 2*sum(integrand[seq(2,nc,by=2)]) + 4 *sum(integrand[seq(3,nc-1, by=2)]) )*step/3
    integrand<-(f0.hat*(m0.s.hat-m1.s.hat))[indc]
    c.hat=(integrand[1] + integrand[nc+1] + 2*sum(integrand[seq(2,nc,by=2)]) + 4 *sum(integrand[seq(3,nc-1, by=2)]) )*step/3

    if (length(ind0)>=5){
      indstar=indc[length(indc)]*( mean(abs(indc[1]-ind0))>=mean(abs(indc[length(indc)]-ind0)) )+
        indc[1]*( mean(abs(indc[1]-ind0))<mean(abs(indc[length(indc)]-ind0)) )
      integrand<-(f0.hat)[ind0]
      n0=length(ind0)-1
      K1=(integrand[1] + integrand[n0+1] + 2*sum(integrand[seq(2,n0,by=2)]) + 4 *sum(integrand[seq(3,n0-1, by=2)]) )*step/3
      lambda=(K2+K1*(f0.hat/f1.hat)[indstar])^{-1}*c.hat+K1*(K2+K1*(f0.hat/f1.hat)[indstar])^{-1}*((m0.s.hat-m1.s.hat)[indstar])
      c=(1+(K1/K2)*(f0.hat/f1.hat)[indstar])^{-1}*( (m1.s.hat-m0.s.hat)[indstar]+((f0.hat/f1.hat)[indstar]/K2) *c.hat)
      g.s.hat=m1.s.hat+(f0.hat/f1.hat)*lambda
      sgs=rbind(cbind(s[indc],g.s.hat[indc]),cbind(s[ind1],g.s.hat[ind1]),cbind(s[ind0],(m0.s.hat+c)[ind0]))
      g.s.hat=sgs[order(sgs[,1]),2]
    } else{
      lambda=(K2)^{-1}*c.hat
      g.s.hat=m1.s.hat+(f0.hat/f1.hat)*lambda
      sgs=rbind(cbind(s[indc],g.s.hat[indc]),cbind(s[ind1],g.s.hat[ind1]),cbind(s[ind0],(m0.s.hat)[ind0]))
      g.s.hat=sgs[order(sgs[,1]),2]
    }

    out=cbind(sort(sgs[,1]),g.s.hat)
  }


  gen.perturb.weights=function(data.num, n, num.perturb=500){
    set.seed(data.num)
    #matrix(rexp(n*num.perturb, rate=1), nrow=n, ncol=num.perturb)
    index = sapply(1:num.perturb,function(x) sample(1:n,n,replace=T))
    apply(index,2,function(x) tabulate(x,nbins=n))
  }


  #intcox package cannot handle weights, so we used bootstrap
  gen.bootstrap.weights=function(data.num, n, num.perturb=500){
    set.seed(data.num)
    sapply(1:num.perturb,function(x) sample(1:n,n,replace=T))
    #weights = apply(index,2,function(x) tabulate(x,nbins=n))
    #list(index=index,weights=weights)
  }


  gopt.perturb=function(yob, sob, aob, n,v){
    if (min(range(sob[aob==1])==range(sob))>0){ supp1=range(sob[aob==1])
    }else if (min(sob[aob==1])-min(sob)==0){supp1=c(quantile(sob[aob==1],0),quantile(sob[aob==1],0.99))
    }else if (max(sob[aob==1])-max(sob)==0){ supp1=c(quantile(sob[aob==1],0.01),quantile(sob[aob==1],1))
    } else {print("You need to switch the groups.")}
    ind1=which.min(abs(supp1[1]-s)):which.min(abs(supp1[2]-s))
    supp0=range(sob[aob==0])
    ind0=which.min(abs(supp0[1]-s)):which.min(abs(supp0[2]-s))
    indc=intersect(ind1,ind0)
    ind1=setdiff(ind1,indc)
    ind0=setdiff(ind0,indc)

    bw = 1.06*sd(sob)*n^(-1/5)/(n^0.06)
    kern = Kern.FUN(zz=s,zi=sob,bw)
    f0.hat=apply(v*kern*(1-aob),2,mean)/mean(v*(1-aob))
    f1.hat=apply(v*kern*aob,2,mean)/mean(v*aob)
    m1.s.hat=apply(v*yob*kern*aob,2,sum)/apply(v*kern*aob,2,sum)
    m0.s.hat=apply(v*yob*kern*(1-aob),2,sum)/apply(v*kern*(1-aob),2,sum)
    integrand<-(f0.hat^2/f1.hat)[indc]
    nc=length(indc)-1
    K2=(integrand[1] + integrand[nc+1] + 2*sum(integrand[seq(2,nc,by=2)]) + 4 *sum(integrand[seq(3,nc-1, by=2)]) )*step/3
    integrand<-(f0.hat*(m0.s.hat-m1.s.hat))[indc]
    c.hat=(integrand[1] + integrand[nc+1] + 2*sum(integrand[seq(2,nc,by=2)]) + 4 *sum(integrand[seq(3,nc-1, by=2)]) )*step/3

    if (length(ind0)>=5){
      indstar=indc[length(indc)]*( mean(abs(indc[1]-ind0))>=mean(abs(indc[length(indc)]-ind0)) )+
        indc[1]*( mean(abs(indc[1]-ind0))<mean(abs(indc[length(indc)]-ind0)) )
      integrand<-(f0.hat)[ind0]
      n0=length(ind0)-1
      K1=(integrand[1] + integrand[n0+1] + 2*sum(integrand[seq(2,n0,by=2)]) + 4 *sum(integrand[seq(3,n0-1, by=2)]) )*step/3
      lambda=(K2+K1*(f0.hat/f1.hat)[indstar])^{-1}*c.hat+K1*(K2+K1*(f0.hat/f1.hat)[indstar])^{-1}*((m0.s.hat-m1.s.hat)[indstar])
      c=(1+(K1/K2)*(f0.hat/f1.hat)[indstar])^{-1}*( (m1.s.hat-m0.s.hat)[indstar]+((f0.hat/f1.hat)[indstar]/K2) *c.hat)
      g.s.hat=m1.s.hat+(f0.hat/f1.hat)*lambda
      sgs=rbind(cbind(s[indc],g.s.hat[indc]),cbind(s[ind1],g.s.hat[ind1]),cbind(s[ind0],(m0.s.hat+c)[ind0]))
      g.s.hat=sgs[order(sgs[,1]),2]
    } else{
      lambda=(K2)^{-1}*c.hat
      g.s.hat=m1.s.hat+(f0.hat/f1.hat)*lambda
      sgs=rbind(cbind(s[indc],g.s.hat[indc]),cbind(s[ind1],g.s.hat[ind1]),cbind(s[ind0],(m0.s.hat)[ind0]))
      g.s.hat=sgs[order(sgs[,1]),2]
    }

    out=cbind(sort(sgs[,1]),g.s.hat)
  }


  resam6<- function(v,data){
    n=nrow(data)
    indexindex=sample(n, n/2, replace = FALSE)
    data1=data[indexindex,]
    data2=data[-indexindex,]
    n1=nrow(data1);n2=nrow(data2)

    v1=v[1:n1]; v2=v[(n1+1):(n1+n2)]
    ######## re2
    ####
    yob=data1[,1]; sob=data1[,2]; aob=data1[,3]; n=n1; v=v1
    from = min(sob); to = max(sob); step=((to - from)/nn); s=seq(from, to, by = step)
    out=tryCatch(gopt.perturb(yob=data1[,1], sob=data1[,2], aob=data1[,3], n=n1,v=v1),
                 error = function(e)(NA) )
    if (any(is.na(out))) {return(rep(NA,1+length(ncut)))}
    s=out[,1]
    g.s.hat=out[,2]

    ####
    yob=data2[,1]; sob=data2[,2]; aob=data2[,3]; n=n2; v=v2
    mu0=mean(as.numeric(v)*yob*(1-aob))/mean(as.numeric(v)*(1-aob)); mu1=mean(as.numeric(v)*yob*(aob))/mean(as.numeric(v)*aob)
    sig=sqrt( 4*(mean( as.numeric(v)*aob*(yob-mu1)^2 ) + mean( as.numeric(v)*(1-aob)*(yob-mu0)^2 )) )
    tempind=c(sapply(1:n, function(kk){which.min(abs(sob[kk]-s))}))
    mu0g=mean(as.numeric(v)*g.s.hat[tempind]*(1-aob))/mean(as.numeric(v)*(1-aob)); mu1g=mean(as.numeric(v)*g.s.hat[tempind]*aob)/mean(as.numeric(v)*aob)
    sig.g=sqrt( 4*(mean( as.numeric(v)*aob*(g.s.hat[tempind]-mu1g)^2 ) + mean( as.numeric(v)*(1-aob)*(g.s.hat[tempind]-mu0g)^2 )) )

    ## PTE RE
    causal=mean(as.numeric(v)*yob*aob)/mean(as.numeric(v)*aob)-mean(as.numeric(v)*yob*(1-aob))/mean(as.numeric(v)*(1-aob))
    causals=mean(as.numeric(v)*g.s.hat[tempind]*aob)/mean(as.numeric(v)*aob)-mean(as.numeric(v)*g.s.hat[tempind]*(1-aob))/mean(as.numeric(v)*(1-aob))

    pte2new=causals/causal
    re2new.es=rep(NA,length(ncut))
    for (j in 1:length(ncut)){
      re2new.es[j]=(1-( pnorm(1.96-sqrt(ncut[j])*causals/sig.g)-pnorm(-1.96-sqrt(ncut[j])*causals/sig.g) ) )/
        (1-( pnorm(1.96-sqrt(ncut[j])*causal/sig)-pnorm(-1.96-sqrt(ncut[j])*causal/sig) ) )
    }

    ######## re1
    ####
    yob=data2[,1]; sob=data2[,2]; aob=data2[,3]; n=n2; v=v2
    from = min(sob); to = max(sob); step=((to - from)/nn); s=seq(from, to, by = step)
    out=tryCatch(gopt.perturb(yob=data2[,1], sob=data2[,2], aob=data2[,3], n=n2,v=v2),
                 error = function(e)(NA) )
    if (any(is.na(out))) {return(rep(NA,1+length(ncut)))}
    s=out[,1]
    g.s.hat=out[,2]

    ####
    yob=data1[,1]; sob=data1[,2]; aob=data1[,3]; n=n1; v=v1
    mu0=mean(as.numeric(v)*yob*(1-aob))/mean(as.numeric(v)*(1-aob)); mu1=mean(as.numeric(v)*yob*(aob))/mean(as.numeric(v)*aob)
    sig=sqrt( 4*(mean( as.numeric(v)*aob*(yob-mu1)^2 ) + mean( as.numeric(v)*(1-aob)*(yob-mu0)^2 )) )
    tempind=c(sapply(1:n, function(kk){which.min(abs(sob[kk]-s))}))
    mu0g=mean(as.numeric(v)*g.s.hat[tempind]*(1-aob))/mean(as.numeric(v)*(1-aob)); mu1g=mean(as.numeric(v)*g.s.hat[tempind]*aob)/mean(as.numeric(v)*aob)
    sig.g=sqrt( 4*(mean( as.numeric(v)*aob*(g.s.hat[tempind]-mu1g)^2 ) + mean( as.numeric(v)*(1-aob)*(g.s.hat[tempind]-mu0g)^2 )) )

    ## PTE RE
    causal=mean(as.numeric(v)*yob*aob)/mean(as.numeric(v)*aob)-mean(as.numeric(v)*yob*(1-aob))/mean(as.numeric(v)*(1-aob))
    causals=mean(as.numeric(v)*g.s.hat[tempind]*aob)/mean(as.numeric(v)*aob)-mean(as.numeric(v)*g.s.hat[tempind]*(1-aob))/mean(as.numeric(v)*(1-aob))

    pte1new=causals/causal
    re1new.es=rep(NA,length(ncut))
    for (j in 1:length(ncut)){
      re1new.es[j]=(1-( pnorm(1.96-sqrt(ncut[j])*causals/sig.g)-pnorm(-1.96-sqrt(ncut[j])*causals/sig.g) ) )/
        (1-( pnorm(1.96-sqrt(ncut[j])*causal/sig)-pnorm(-1.96-sqrt(ncut[j])*causal/sig) ) )
    }

    ptenew=(pte1new+pte2new)/2
    renew.es=(re1new.es+re2new.es)/2

    out=c(ptenew,renew.es)

  }


  # main function
  # PTERE=function(data, ncut=c(50,100,150,200,500,1000) ){
  n=nrow(data)
  set.seed(12345678)
  indexindex=sample(n, n/2, replace = FALSE)
  data1=data[indexindex,]
  data2=data[-indexindex,]
  n1=nrow(data1);n2=nrow(data2)
  nn=200

  ######## re2
  ####
  yob=data1[,1]; sob=data1[,2]; aob=data1[,3]; n=n1
  from = min(sob); to = max(sob); step=((to - from)/nn); s=seq(from, to, by = step)
  temp=tryCatch(gopt(yob=data1[,1], sob=data1[,2], aob=data1[,3], n=n1),
                error = function(e)(NA) )
  if (any(is.na(temp))) {return(NULL)}
  s=temp[,1]
  g.s.hat=temp[,2]

  ####
  yob=data2[,1]; sob=data2[,2]; aob=data2[,3]; n=n2
  mu0=mean(yob*(1-aob))/mean(1-aob); mu1=mean(yob*(aob))/mean(aob)
  sig=sqrt( 4*(mean( aob*(yob-mu1)^2 ) + mean( (1-aob)*(yob-mu0)^2 )) )
  tempind=c(sapply(1:n, function(kk){which.min(abs(sob[kk]-s))}))
  mu0g=mean(g.s.hat[tempind]*(1-aob))/mean(1-aob); mu1g=mean(g.s.hat[tempind]*aob)/mean(aob)
  sig.g=sqrt( 4*(mean( aob*(g.s.hat[tempind]-mu1g)^2 ) + mean( (1-aob)*(g.s.hat[tempind]-mu0g)^2 )) )

  ## PTE RE
  causal=mean(yob*aob)/mean(aob)-mean(yob*(1-aob))/mean(1-aob)
  causals=mean(g.s.hat[tempind]*aob)/mean(aob)-mean(g.s.hat[tempind]*(1-aob))/mean(1-aob)
  pte2new.es=causals/causal
  re2new.es=rep(NA,length(ncut))
  for (j in 1:length(ncut)){
    re2new.es[j]=(1-( pnorm(1.96-sqrt(ncut[j])*causals/sig.g)-pnorm(-1.96-sqrt(ncut[j])*causals/sig.g) ) )/
      (1-( pnorm(1.96-sqrt(ncut[j])*causal/sig)-pnorm(-1.96-sqrt(ncut[j])*causal/sig) ) )
  }

  ######## re1
  ####
  yob=data2[,1]; sob=data2[,2]; aob=data2[,3]; n=n2
  from = min(sob); to = max(sob); step=((to - from)/nn); s=seq(from, to, by = step)
  temp=tryCatch(gopt(yob=data2[,1], sob=data2[,2], aob=data2[,3], n=n2),
                error = function(e)(NA) )
  if (any(is.na(temp))) {return(NULL)}
  s=temp[,1]
  g.s.hat=temp[,2]

  ####
  yob=data1[,1]; sob=data1[,2]; aob=data1[,3]; n=n1
  mu0=mean(yob*(1-aob))/mean(1-aob); mu1=mean(yob*(aob))/mean(aob)
  sig=sqrt( 4*(mean( aob*(yob-mu1)^2 ) + mean( (1-aob)*(yob-mu0)^2 )) )
  tempind=c(sapply(1:n, function(kk){which.min(abs(sob[kk]-s))}))
  mu0g=mean(g.s.hat[tempind]*(1-aob))/mean(1-aob); mu1g=mean(g.s.hat[tempind]*aob)/mean(aob)
  sig.g=sqrt( 4*(mean( aob*(g.s.hat[tempind]-mu1g)^2 ) + mean( (1-aob)*(g.s.hat[tempind]-mu0g)^2 )) )

  ## PTE RE
  causal=mean(yob*aob)/mean(aob)-mean(yob*(1-aob))/mean(1-aob)
  causals=mean(g.s.hat[tempind]*aob)/mean(aob)-mean(g.s.hat[tempind]*(1-aob))/mean(1-aob)
  pte1new.es=causals/causal
  re1new.es=rep(NA,length(ncut))
  for (j in 1:length(ncut)){
    re1new.es[j]=(1-( pnorm(1.96-sqrt(ncut[j])*causals/sig.g)-pnorm(-1.96-sqrt(ncut[j])*causals/sig.g) ) )/
      (1-( pnorm(1.96-sqrt(ncut[j])*causal/sig)-pnorm(-1.96-sqrt(ncut[j])*causal/sig) ) )
  }

  ptenew.es=(pte1new.es+pte2new.es)/2
  renew.es=(re2new.es+re1new.es)/2

  #### variance
  n=n1+n2
  re=n.resam
  v=matrix(rexp(n*re),nrow=n)
  g.s.re=apply(v,2,resam6,data)
  temp=g.s.re[1,];ind=which((temp<1)*(temp>0)>0)
  se=c(apply(g.s.re[,ind],1,sd,na.rm=T) )

  ##
  out=data.frame( (cbind(ptenew.es,renew.es[1],renew.es[2],renew.es[3],renew.es[4],renew.es[5],renew.es[6],
                         se[1],se[2],se[3],se[4],se[5],se[6],se[7] )) )
  colnames(out)[2:7]=paste0('rp',ncut)
  colnames(out)[8]=paste0('pte.se')
  colnames(out)[9:14]=paste0('rp.se',ncut)
  return(out)
}
