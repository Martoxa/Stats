## Calculates an Empirical Orthogonal Function analysis. Inputs are the analysis matrix (Mx), and if the data is scaled, FALSE is default
EOF<-function(Mx,scl=FALSE){
  Mx<-scale(Mx, center = TRUE, scale = scl)
  c_Mx<-cov(Mx)
  svd_Mx<-svd(c_Mx)
  signalsMx<-Mx%*%svd_Mx$v
  explainedMx<-100*(svd_Mx$d^2/sum(svd_Mx$d^2))
  out<-list("SVD"=svd_Mx,"Signals"=signalsMx,"Explained"=explainedMx)
  return(out)
}

## 
RuleN<-function(Mx,iter,alph){
  simulatedEig<-as.data.frame(replicate(dim(Mx)[2],rnorm(iter)))
  simulatedExp<-as.data.frame(replicate(dim(Mx)[2],rnorm(iter)))
  for(i in 1:iter){
    rnM<-replicate(dim(Mx)[2],rnorm(dim(Mx)[1]))
    rnM<-rnM*sqrt(diag(var(Mx)))
    Modesrn<-EOF(rnM)
    simulatedEig[i,]<-Modesrn$SVD$d
    simulatedExp[i,]<-Modesrn$Explained
  }
  simulatedEig<-apply(simulatedEig,2,sort,decreasing=FALSE)
  simulatedExp<-apply(simulatedExp,2,sort,decreasing=FALSE)
  alp<-floor(alph*iter)
  out<-list("THEig"=simulatedEig[alp,],"THExp"=simulatedExp[alp,])
  return(out)
}
