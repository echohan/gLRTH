% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gLRTH.R
\name{gLRTH}
\alias{gLRTH}
\title{The Function for the Likelihood Ratio Test for Genome-Wide Association under Genetic Heterogeneity with genotype frequencies as input values}
\usage{
gLRTH(n0, n1, n2, m0, m1, m2)
}
\arguments{
\item{n0}{AA genotype frequency in case}

\item{n1}{Aa genotype frequency in case}

\item{n2}{aa genotype frequency in case}

\item{m0}{AA genotype frequency in control}

\item{m1}{Aa genotype frequency in control}

\item{m2}{aa genotype frequency in control}
}
\value{
The test statistic and asymptotic p-value for the likelihood ratio test under genetic heterogeneity
}
\description{
The Function for the Likelihood Ratio Test for Genome-Wide Association under Genetic Heterogeneity with genotype frequencies as input values
}
\examples{
gLRTH(n0=2940, n1=738, n2=53, m0=3601, m1=1173, m2=117)
}
\author{
Xiaoxia Han and Yongzhao Shao
}
\references{
Qian M., Shao Y. (2013) A Likelihood Ratio Test for Genome-Wide Association under Genetic Heterogeneity.
Annals of Human Genetics, 77(2): 174-182.
}

