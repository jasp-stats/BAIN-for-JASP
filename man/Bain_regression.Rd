\name{Bain_regression}
\alias{Bain_regression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Bayes factors of informative hypotheses for regression models
}
\description{
  This function computes approximated adjusted fractional Bayes factors for hypothesis tests in regression models.
}
\usage{
Bain_regression(formula, data, ERr = NULL, IRr = NULL, ...,
                covariates_hypo = NULL, standardize = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    A formula that specifies the regression model. See \code{\link{formula}} or \code{\link{lm}} for details.
  }

  \item{data}{
    A data frame that contains the variables in the regression model.
  }

  \item{ERr}{
    A matrix representing equality constraints in a hypothesis. If no equality constraints, set ERr = NULL or ERr = matrix(0,0,0).
    See \code{\link{Bain}} for details.
  }
  \item{IRr}{
    A matrix representing inequality or about equality constraints in a hypothesis. If no inequality and about equality constraints, set IRr = NULL or IRr=matrix(0,0,0). See \code{\link{Bain}} for details.
  }
  \item{\dots}{
    Sets of matrices representing the equality and/or inequality constraints for the extra informative hypotheses. For example, if there are two informative hypotheses, then ... will be two matrices for the second hypothesis: ERr2 and IRr2. See \code{\link{Bain}} for details.
  }

  \item{covariates_hypo}{
    A character vector or a formula that specifies names of predictors of which the coefficients are in the hypotheses. For example, if the hypothesis compares cofficients of predictors "Age" and "Sex" of the dependent variable "Income", then covariates_hypo = c("Age", "Sex") or covariates_hypo = ~ Age + Sex. Please note that the order of the variables specified by covariates_hypo should be in line with the restriction matrices specified by ERr and IRr. For example, if covariates_hypo = ~ Age + Sex, then the first column in ERr or IRr should represent the coefficient of Age and the second should represent the coefficient of Sex.
The defaul value is covariates_hypo = NULL, which assumes all predictors specified in formula are used in the hypotheses.
  }

  \item{standardize}{
    logical. If FALSE (default), then the estimates and covariance matrix of coefficients will not be standardized. If TRUE, then the estimates and covariance matrix of coefficients will be standardsized. When testing whether coefficients are zero, positive, or negative, we suggest standardize = FALSE. When comparing coefficients, we suggest standardize = TRUE. The standardization uses function "seBeta" in package "fungible", see reference for the theory.
}

}

\value{
  \item{fit}{Fit for each hypothesis under comparison}
  \item{complexity}{Complexity for each hypothesis under comparison}
  \item{BF}{Bayes factor for each hypothesis compared to its complement}
  \item{PMPa}{Posterior probability for each hypothesis excluding unconstrained hypothesis}
  \item{PMPb}{Posterior probability for each hypothesis including unconstrained hypothesis}
}


\references{
Jones, J.A. and Waller, N.G. (2015). The normal theory and asymptotic distribution-free (ADF) covariance matrix of standardized regression coefficients: theorical extensions and finite sample behavior. Psychometrika, 80, 365-378.
}

\author{
  Xin Gu, Herbert Hoijtink, Joris Mulder
}
%%\note{
  %%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%  \seealso{
    %% ~~objects to See Also as \code{\link{help}}, ~~~
%%  }
\examples{
  ##Example 1
  #Input
  data(cars)

  #Hypothesis
  #H0: theta = 0
  ERr<-matrix(c(1,0),nrow=1,ncol=2,byrow = TRUE)
  IRr<-NULL

  res<-Bain_regression(dist ~ speed, data = cars, ERr, IRr) #run
  #Results are printed.
  #Results for fit, complexity, Bayes factor, and PMPs are saved in "res":

  plot(res)
  #Results for PMPs are plotted.


  ##Example 2
  #Input
  data(LifeCycleSavings)

  #Hypothesis 1
  #H1: theta1 = theta2 = theta3 = 0
  ERr1<-matrix(c(1,0,0,0,0,1,0,0,0,0,1,0),nrow=3,ncol=4,byrow = TRUE)
  IRr1<-NULL

  #H2: theta1 < 0, theta2 = 0, theta3 > 0
  ERr2<-matrix(c(0,1,0,0),nrow=1,ncol=4,byrow = TRUE)
  IRr2<-matrix(c(-1,0,0,0,0,0,1,0),nrow=2,ncol=4,byrow = TRUE)

  #H3: theta1 < theta2 < theta3
  ERr3<-NULL
  IRr3<-matrix(c(-1,1,0,0,0,-1,1,0),nrow=2,ncol=4,byrow = TRUE)

  #run
  res<-Bain_regression(sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings,
  ERr1, IRr1, ERr2, IRr2, ERr3, IRr3,
  covariates_hypo = ~ pop15 + dpi + ddpi, standardize = TRUE)
  ## assume "pop75" is not in the hypothesis.
  #Results are printed.
  #Results for fit, complexity, Bayes factor, and PMPs are saved in "res":

  plot(res)
  #Results for PMPs are plotted.

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


\alias{Bain_regression}



