\name{Bain_ttestData}
\alias{Bain_ttestData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Bayesian t test using approximate adjusted fractional Bayes factors
}
\description{
  This function computes approximated adjusted fractional Bayes factors one-sample or two-sample t test.
}
\usage{
Bain_ttestData(x, y = NULL, nu = 0, type = 1, paired=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A vector for the sample
  }

  \item{y}{
    A vector for the sample in the second group (optional).
  }

  \item{nu}{
   A numeric number for the null value in t test. Default number is nu = 0.
  }

  \item{type}{
    An integer number (1-4) indicating the type of the t test.
    "1" means theta = 0 vs theta != 0.
    "2" means theta = 0 vs theta > 0.
    "3" means theta = 0 vs theta < 0.
    "4" means theta > 0 vs theta < 0.
    "5" means theta = 0 vs theta > 0 vs theta < 0.
    "theta" is the mean parameter of the sample.
  }
  \item{paired}{
   Paired = TRUE for paired t test. Default is paried = FALSE.
  }
}

\value{
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  Return a list that contains Bayes factor and posterior model probabilities.

  %% ...
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
  #Example 1
  #One sample t test:
  #H0: theta=0
  x<-c(-1,1,2,-1,-0.5)
  res<-Bain_ttestData(x)

  plot(res)
  #Results for PMPs are plotted.

  #Other types
  Bain_ttestData(x,type=2)
  Bain_ttestData(x,type=3)
  Bain_ttestData(x,type=4)
  Bain_ttestData(x,type=5)


  #Example2
  #Two sample t test:
  #H0: theta1=theta2
  x<-c(-1,1,2,-1,-0.5)
  y<-c(0,1,2,0,1)
  res<-Bain_ttestData(x,y)

  plot(res)
  #Results for PMPs are plotted.


  #Example3
  #paired t test:
  #H0: theta1=theta2
  x<-c(-1,1,2,-1,-0.5)
  y<-c(0,1,2,0,1)
  res<-Bain_ttestData(x,y,paired=TRUE)

  plot(res)
  #Results for PMPs are plotted.


}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


\alias{Bain_ttestData}



