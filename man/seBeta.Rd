\name{seBeta}
\alias{seBeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Standard Errors and CIs for Standardized Regression Coefficients
}
\description{
  Computes Normal Theory and ADF Standard Errors and CIs for Standardized Regression Coefficients
}
\usage{
  seBeta(X, y, cov.x = NULL, cov.xy = NULL, var.y = NULL, Nobs = NULL,
         alpha = 0.05, estimator = "ADF", digits = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
    Matrix of predictor scores.
  }
  \item{y}{
    Vector of criterion scores.
  }
  \item{cov.x}{
    Covariance or correlation matrix of predictors.
  }
  \item{cov.xy}{
    Vector of covariances or correlations between predictors and criterion.
  }
  \item{var.y}{
    Criterion variance.
  }
  \item{Nobs}{
    Number of observations.
  }
  \item{alpha}{
    Desired Type I error rate; default = .05.
  }
  \item{estimator}{
    'ADF' or 'Normal' confidence intervals - requires raw X and raw y; default = 'ADF'.
  }
  \item{digits}{
    Number of significant digits to print; default = 3.
  }

}

\value{
  \item{cov.Beta}{
    Normal theory or ADF covariance matrix of standardized regression coefficients.
  }
  \item{se.Beta}{
    standard errors for standardized regression coefficients.
  }
  \item{alpha}{
    desired Type-I error rate.
  }
  \item{CI.Beta}{
    Normal theory or ADF (1-alpha)% confidence intervals for standardized regression coefficients.
  }
  \item{estimator}{
    estimator = "ADF" or "Normal".
  }


}

\author{
Jeff Jones and Niels Waller
}

\references{
Jones, J. A, and Waller, N. G. (2015). The Normal-Theory and Asymptotic Distribution-Free (ADF) covariance matrix of standardized regression coefficients: Theoretical extensions and finite sample behavior. Psychometrika, 80, 365-378.
}
%%\note{
  %%  ~~further notes~~
    %%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

  %%\seealso{
    %% ~~objects to See Also as \code{\link{help}}, ~~~
      %%}
\examples{
library(MASS)

set.seed(123)

R <- matrix(.5, 3, 3)
diag(R) <- 1
X <- mvrnorm(n = 200, mu = rep(0, 3), Sigma = R, empirical = TRUE)
Beta <- c(.2, .3, .4)
y <- X \%*\% Beta + .64 * scale(rnorm(200))
seBeta(X, y, Nobs = 200, alpha = .05, estimator = 'ADF')

# 95\% CIs for Standardized Regression Coefficients:
#
#        lbound estimate ubound
# beta_1  0.104    0.223  0.341
# beta_2  0.245    0.359  0.473
# beta_3  0.245    0.360  0.476
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Statistics }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


\alias{seBeta}



