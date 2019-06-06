\name{Bain_ancova_cm}
\alias{Bain_ancova_cm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Bayes factors of informative hypotheses for ANCOVA models (user-friendly hypotheses input)
}
\description{
  This function computes approximated adjusted fractional Bayes factors for ancova test
  with user-friendly hypotheses input.
}
\usage{
Bain_ancova_cm(X, dep_var, covariates, group, hyp)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{X}{
    A data frame that contains dependent variable, covariates, and group factor.
  }

  \item{dep_var}{
    A character that indicates the name of the dependent variable in X. Only one dependent variable can be specified.
  }

  \item{covariates}{
    A character vector that indicates names of covariates in X.
  }

  \item{group}{
    A character that indicates the name of the group factor in X. Only one group factor can be specified.
  }

  \item{hyp}{
    A character string that specifies hypotheses under evaulation. Variable names in the hypotheses have to be specified as 'group.n' with n = 1,2,3... denoting group numbers (e.g., group.1, group.2 and group.3), or 'group.category' (e.g., Sex.M and Sex.F for group factor 'Sex' with two categories 'M' and 'F'). Group/categpory names can contain: letters, numbers, and _. Variable names in one constraint have to be linked by +, -, =, >, <, e.g., group.1=group.2>group.3, or Sex.M<Sex.F. Different constraints in one hypothesis can be separated by "&", e.g., hyp = "Sex.F>0&Sex.M<0". Competing hypotheses can be separated by ';', e.g., hyp = "Sex.M=Sex.F; Sex.M>Sex.F". See examples and \code{\link{create_matrices}} for details.
  }
}

\value{
  \item{fit}{Fit for each hypothesis under comparison}
  \item{complexity}{Complexity for each hypothesis under comparison}
  \item{BF}{Bayes factor for each hypothesis compared to its complement}
  \item{PMPa}{Posterior probability for each hypothesis excluding unconstrained hypothesis}
  \item{PMPb}{Posterior probability for each hypothesis including unconstrained hypothesis}
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
  data(ChickWeight)
  #Hypotheses
  hyp<-"Diet.1=Diet.2=Diet.3=Diet.4;
        Diet.1<Diet.2<Diet.3<Diet.4;
        Diet.1<Diet.2<Diet.4<Diet.3"

  res<-Bain_ancova_cm(ChickWeight,"weight","Time","Diet",hyp) #run
  #Results are printed.
  #Results for fit, complexity, Bayes factor, and PMPs are saved in "res":

  plot(res)
  #Results for PMPs are plotted.

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


\alias{Bain_ancova_cm}



