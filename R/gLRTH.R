gLRTH<-function(n0, n1, n2, m0, m1, m2, ...){
  n<-n0+n1+n2
  m<-m0+m1+m2
  p0<-(n2+m2+(n1+m1)/2)/(n+m)
  pH<-(m2+m1/2)/m
  pD<-(n2+n1/2)/n
  
  LRT.stat<-ifelse(4*n0*n2>n1^2, 
                   2*log(((1-pH)/(1-p0))^(2*m0)*(pH*(1-pH)/p0/(1-p0))^m1*(pH/p0)^(2*m2)*
                           (n0/(n*(1-p0)^2))^n0*(n1/(n*2*p0*(1-p0)))^n1*(n2/(n*p0^2))^n2),
                   2*log(((1-pH)/(1-p0))^(2*m0)*(pH*(1-pH)/p0/(1-p0))^m1*(pH/p0)^(2*m2)*
                           ((1-pD)/(1-p0))^(2*n0)*(pD*(1-pD)/p0/(1-p0))^n1*(pD/p0)^(2*n2)))
  
  pval<-(pchisq(LRT.stat, 1, lower.tail = F) + pchisq(LRT.stat, 2, lower.tail = F))/2
  
  return(list(chisq.stat=LRT.stat, pval=pval))
}