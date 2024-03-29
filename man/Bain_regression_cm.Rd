\name{Bain_regression_cm}
\alias{Bain_regression_cm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Bayes factors of informative hypotheses for regression models
}
\description{
  This function computes approximated adjusted fractional Bayes factors for hypothesis tests in regression models.
}
\usage{
  Bain_regression_cm(formula, data, hyp, standardize = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    A formula that specifies the regression model. See \code{\link{formula}} or \code{\link{lm}} for details.
  }

  \item{data}{
    A data frame that contains the variables in the regression model.
  }

  \item{hyp}{
    A character string that specifies hypotheses under evaulation. Variable names should match the names of regression predictors. For example, if formula = Y ~ X1 + X2 + X3, variable names should consist of "X1", "X2" and "X3", but cannot contain e.g., "X4". Variable names in one constraint have to be linked by +, -, =, >, <. Different constraints in one hypothesis can be separated by "&". Competing hypotheses can be separated by ';'. For example, hyp = "X1>0&X2>X3;X1=X2=X3=0". See examples and \code{\link{create_matrices}} for details.
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

  hyp = "speed = 0"
  res<-Bain_regression_cm(dist ~ speed, cars, hyp) #run
  #Results are printed.
  #Results for fit, complexity, Bayes factor, and PMPs are saved in "res":

  plot(res)
  #Results for PMPs are plotted.


  ##Example 2
  #Input
  data(LifeCycleSavings)

  hyp<-"pop15 = dpi = ddpi = 0;
        pop15 < 0 & dpi = 0 & ddpi > 0;
        pop15 < dpi < ddpi"

  #run
  res<-Bain_regression_cm(sr ~ pop15 + pop75 + dpi + ddpi, LifeCycleSavings,
                          hyp, standardize = TRUE)
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


\alias{Bain_regression_cm}



