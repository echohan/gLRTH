\name{gLRTH2}
\alias{gLRTH2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Function for the Likelihood Ratio Test for Genome-Wide Association under Genetic Heterogeneity with long format of genotype and disease status data
}
\description{
Returns test statistic and asymptotic p-value for the likelihood ratio test under genetic heterogeneity proposed by Qian and Shao 2013.
}
\usage{
gLRTH2(g, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{a n x 1 vector contains the number of minor allele (i.e., 0, 1 or 2) for a given SNP, where n is the number of observations.
}
  \item{d}{a n x 1 vector contains disease status, 1 is case and 0 is control, where n is the number of observations.
}
}
\details{
Missing values in g (genotype) or d (disease status)  will be removed.
}
\value{
\item{statistic}{The test statistic for the likelihood ratio test under genetic heterogeneity}
\item{pval}{The asymptotic p-value for the for the likelihood ratio test under genetic heterogeneity}
}
\references{
Qian, M., Shao, Y. (2013) A Likelihood Ratio Test for Genome-Wide Association under Genetic Heterogeneity
}
\author{
Xiaoxia Han and Yongzhao Shao
}


\examples{
disease<-c(rep(1, 200), rep(0, 200))
geno1<-c(rbinom(n=50, size=2, prob=0.5), rbinom(n=150, size=2, prob=0.23))
geno2<-rbinom(n=200, size=2, prob=0.5)
geno<-c(geno1, geno2)
gLRTH2(g=geno, d=disease)
}
