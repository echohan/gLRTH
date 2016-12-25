\name{gLRTH}
\alias{gLRTH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Function for the Likelihood Ratio Test for Genome-Wide Association under Genetic Heterogeneity with genotype frequencies as input values
}
\description{
Returns test statistic and asymptotic p-value for the likelihood ratio test under genetic heterogeneity proposed by Qian and Shao 2013.
}
\usage{
gLRTH(n, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n0}{AA genotype frequency in case
}
 \item{n1}{Aa genotype frequency in case
}
 \item{n2}{aa genotype frequency in case
}
  \item{m0}{AA genotype frequency in control
}
  \item{m1}{Aa genotype frequency in control
}
  \item{m0}{aa genotype frequency in control
}
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
gLRTH(n0=2940, n1=738, n2=53, m0=3601, m1=1173, m2=117)
}


