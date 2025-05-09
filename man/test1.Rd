% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test1.R
\name{test1}
\alias{test1}
\title{Basic draft of Test1}
\usage{
test1(x, setting = "u", c = 0.3, B = 1000)
}
\arguments{
\item{x}{numeric vector. The input data of observed sample, assumed to be sorted.}

\item{setting}{character. Specifies the null distribution. Currently only
`'u'` (Uniform) is supported (default 'u').}

\item{c}{numeric. Tuning parameter for penalty terms (default 0.3).}

\item{B}{integer. Number of Monte Carlo replicates for p-value calculation
(default 1000). Higher values increase precision but take longer.}
}
\value{
Test results, an object of class `"htest"` containing the following components:
  \item{statistic}{The value of the test statistic (Lambda).}
  \item{p.value}{The p-value for the test, computed via simulation.}
  \item{method}{A character string indicating the test performed.}
  \item{data.name}{A character string giving the name of the data.}
  \item{parameters}{A named numeric vector containing `c` (tuning parameter).}
  \item{null.value}{Description of the null hypothesis distribution.}
}
\description{
Performs the P-test in the paper "(1999) Testing uniformity versus a monotone density".
  Calculates a test statistic (Lambda) and a p-value.
}
\examples{
set.seed(123)
# Data close to uniform
uniform_data <- sort(runif(50))
test_uniform <- test1(uniform_data, c = 0.3, B = 1000)
print(test_uniform)

# Data with a monotone trend (e.g., Beta(2,1))
monotone_data <- sort(rbeta(50, 2, 1))
test_monotone <- test1(monotone_data, c = 0.3, B = 1000)
print(test_monotone)

}
\references{
Jiayang Sun. Michael Woodroofe. "Testing uniformity versus a monotone density."
Ann. Statist. 27 (1) 338 - 360, February 1999. https://doi.org/10.1214/aos/1018031114
}
