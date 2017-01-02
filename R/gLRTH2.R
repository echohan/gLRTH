#' The Function for the Likelihood Ratio Test for Genome-Wide Association under Genetic Heterogeneity with long format of genotype and disease status data
#'
#' @return The test statistic and asymptotic p-value for the likelihood ratio test under genetic heterogeneity
#' @param g an n x 1 vector contains the number of minor allele (i.e., 0, 1 or 2) for a given SNP, where n is the number of observations.
#' @param d an n x 1 vector contains disease status, 1 is case and 0 is control, where n is the number of observations.
#' @author Xiaoxia Han and Yongzhao Shao
#' @references
#' Qian M., Shao Y. (2013) A Likelihood Ratio Test for Genome-Wide Association under Genetic Heterogeneity.
#' Annals of Human Genetics, 77(2): 174-182.
#' @examples
#' disease<-c(rep(1, 200), rep(0, 200))
#' geno1<-c(rbinom(n=50, size=2, prob=0.5), rbinom(n=150, size=2, prob=0.23))
#' geno2<-rbinom(n=200, size=2, prob=0.5)
#' geno<-c(geno1, geno2)
#' gLRTH2(g=geno, d=disease)
#' @export
#' @importFrom stats complete.cases pchisq


gLRTH2<-function(g,d){
  data = data.frame(g = g, d = d)
  data = as.data.frame(data)
  data = data[complete.cases(data), ]
  n0 = nrow(subset(data, g == 0 & d == 1 ))
  n1 = nrow(subset(data, g == 1 & d == 1 ))
  n2 = nrow(subset(data, g == 2 & d == 1 ))
  m0 = nrow(subset(data, g == 0 & d == 0 ))
  m1 = nrow(subset(data, g == 1 & d == 0 ))
  m2 = nrow(subset(data, g == 2 & d == 0 ))

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

  return(list(statistic=LRT.stat, pval=pval))
}
