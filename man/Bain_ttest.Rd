\name{Bain_ttest}
\alias{Bain_ttest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian t test using approximate adjusted fractional Bayes factors
}
\description{
This function computes approximated adjusted fractional Bayes factors for one-sample or two-sample t test.
}
\usage{
Bain_ttest(estimate, variance, n, nu=0, type=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estimate}{
A number (one-sample t test) or a vector (two-sample t test) giving the estimate(s) of the group mean(s).
}

  \item{variance}{
A number (one-sample t test) or a vector (two-sample t test) giving the variance(s) of the group mean(s).
}

  \item{n}{
A positive integer number (one-sample t test) or vector (two-sample t test) indicating the sample size in the data.
 }
  \item{nu}{
A numeric number for the null value in t test. Default number is nu = 0.

}

  \item{type}{
An integer number (1-4) indicating the type of the t test.
"1" means H0: mu = 0 vs Hu: mu != 0.
"2" means H0: mu = 0 vs H1: mu > 0.
"3" means H0: mu = 0 vs H2: mu < 0.
"4" means H1: mu > 0 vs H2: mu < 0.
"5" means H0: mu = 0 vs H1: mu > 0 vs H2: mu < 0.
}

}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
Return a list that contains Bayes factor and posterior model probability for each hypothesis.

%% ...
}

\author{
Xin Gu, Herbert Hoijtink, Joris Mulder
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
#Example 1
#One sample t test:
#H0: mu=0
#Input
estimate<-0.5    #Estimate of theta
variance<-1    #Variance of theta
n<-20      #samplesize
res<-Bain_ttest(estimate,variance,n,nu=0,type=1) #Run

#Output printed
#Hypotheses H0: mu = 0 vs Hu: mu != 0
#
#t test result
# BF_0u   PMP_0   PMP_u
# 3.947   0.798   0.202

plot(res)
#Results for PMPs are plotted.

#Other types
Bain_ttest(estimate,variance,n,nu=0,type=2) #Run
Bain_ttest(estimate,variance,n,nu=0,type=3) #Run
Bain_ttest(estimate,variance,n,nu=0,type=4) #Run
Bain_ttest(estimate,variance,n,nu=0,type=5) #Run


#Example2
#Two sample t test:
#H0: mu1=mu2
#Input
estimate<-c(0.2,0.5)
variance<-c(0.01,0.02)
n<-c(20,40)
res<-Bain_ttest(estimate,variance,nu=0,n) #Run

#Output printed
#Hypotheses H0: mu1 = mu2 vs Hu: mu1 != mu2

#t test result
# BF_0u   PMP_0   PMP_u
# 1.822   0.646   0.354

plot(res)
#Results for PMPs are plotted.
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


\alias{Bain_ttest}



