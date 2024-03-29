\name{Bain_anova_cm}
\alias{Bain_anova_cm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Bayes factors of informative hypotheses for ANOVA models (user-friendly hypotheses input)
}
\description{
  This function computes approximated adjusted fractional Bayes factors for anova test
  with user-friendly hypotheses input.
}
\usage{
Bain_anova_cm(X, dep_var = NULL, group = NULL, hyp)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{X}{
    A data frame that contains dependent variable and group factor.
  }

  \item{dep_var}{
    A character that indicates the name of the dependent variable in X.
  }

  \item{group}{
    A character that indicates the name of the group factor in X.
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
  #Example 1
  #Hypotheses
  #H1: group.1=group.2=group.3   #group means
  #H2: group.1>group.2>group.3
  #H3: group.1<group.2<group.3

  #Input
  X<-data.frame(c(rnorm(20,0,1),rnorm(40,.5,2),rnorm(60,1,5)),
                c(rep(1,20),rep(2,40),rep(3,60)))
  names(X)<-c("dep_var","group")
  head(X)

  hyp<-"group.1=group.2=group.3;
        group.1>group.2>group.3;
        group.1<group.2<group.3"

  res<-Bain_anova_cm(X,"dep_var","group",hyp) #run
  #Results are printed.
  #Results for fit, complexity, Bayes factor, and PMPs are saved in "res":

  plot(res)
  #Results for PMPs are plotted.

  #Example 2
  #Hypotheses

  #Input
  df<-data.frame(c(rbinom(20,10,0.5),rbinom(40,10,0.6),rbinom(30,10,0.6),rbinom(30,10,0.7)),
                c(rep("A",20),rep("B",40),rep("AB",30),rep("O",30)))
  names(df)<-c("ambition","bloodtype")
  head(df)

  hyp<-"bloodtype.O = bloodtype.AB = bloodtype.A = bloodtype.B;
        bloodtype.O > bloodtype.AB > bloodtype.A = bloodtype.B"

  res<-Bain_anova_cm(df,"ambition","bloodtype",hyp) #run
  #Results are printed.
  #Results for fit, complexity, Bayes factor, and PMPs are saved in "res":

  plot(res)
  #Results for PMPs are plotted.

}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


\alias{Bain_anova_cm}



