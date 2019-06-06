\name{Bain_ancova}
\alias{Bain_ancova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Bayes factors of informative hypotheses for ANCOVA models
}
\description{
  This function computes approximated adjusted fractional Bayes factors for ancova test.
}
\usage{
Bain_ancova(X, dep_var, covariates, group,
            ERr = NULL, IRr = NULL, ...)
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

  \item{ERr}{
    A matrix representing equality constraints in a hypothesis. If no equality constraints, set ERr = NULL or ERr = matrix(0,0,0).
    See \code{\link{Bain}} for details.
    Number of columns should equal to number of groups plus 1. If one or more groups are not included in the hypothesis, please specify corresponding columns with zeros.
  }
  \item{IRr}{
    A matrix representing inequality or about equality constraints in a hypothesis. If no inequality and about equality constraints, set IRr = NULL or IRr = matrix(0,0,0). See \code{\link{Bain}} for details.
    Number of columns should equal to number of groups plus 1. If one or more groups are not included in the hypothesis, please specify corresponding columns with zeros.
  }
  \item{\dots}{
    Sets of matrices representing the equality and/or inequality constraints for the extra informative hypotheses. For example, if there are two informative hypotheses, then ... will be two matrices for the second hypothesis: ERr2 and IRr2. See \code{\link{Bain}} for details.
  }
}

\details{
Be aware that the constraints matrices ERr/IRr have to be constructed in accordance with the order of the groups according to the lm function. For example, if the groups are labelled a1, a2, a3, then the order is clear, but if the groups are labelled yellow, blue and gray, then the order is blue, gray, yellow.
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
  #Input
  #Build in data: Weight versus age (Time) of chicks on different diets (Diet).
  data(ChickWeight)

  #Hypotheses

  #H1: theta1=theta2=theta3=theta4
  ERr1<-matrix(c(1,-1,0,0,0,0,1,-1,0,0,0,0,1,-1,0),nrow=3,ncol=5,byrow = TRUE)
  IRr1<-NULL

  #H2: theta1<theta2<theta3<theta4
  ERr2<-NULL
  IRr2<-matrix(c(-1,1,0,0,0,0,-1,1,0,0,0,0,-1,1,0),nrow=3,ncol=5,byrow = TRUE)

  #H3: theta1<theta2<theta4<theta3
  ERr3<-NULL
  IRr3<-matrix(c(-1,1,0,0,0,0,-1,0,1,0,0,0,1,-1,0),nrow=3,ncol=5,byrow = TRUE)

  res<-Bain_ancova(ChickWeight,"weight","Time","Diet",ERr1,IRr1,ERr2,IRr2,ERr3,IRr3) #run
  #Results are printed.
  #Results for fit, complexity, Bayes factor, and PMPs are saved in "res":

  plot(res)
  #Results for PMPs are plotted.
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


\alias{Bain_ancova}



