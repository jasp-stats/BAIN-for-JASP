\name{Bain}
\alias{Bain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayes factors for informative hypotheses
}
\description{
This function computes approximated adjusted fractional Bayes factors (AAFBFs) for equality, inequality, and about equality constrained hypotheses. The function can also be used for computing the AAFBF in the case of multiple groups.
}
\usage{
Bain(estimate, Sigma, grouppara = 0, jointpara = 0, n,
     ERr = NULL, IRr = NULL, ..., seed=100, print = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estimate}{
A vector giving the estimates of the parameters used in the hypotheses.
}

  \item{Sigma}{
A matrix (or list of matrices) giving the covariance matrix of the parameters used in the hypotheses. If grouppara = 0, it is a covariance matrix (of the joint parameters). If grouppara > 0, it is a list (with length P) of covariance matrices for group specific and joint parameters for group 1 until group P, where P denotes the number of groups.
}

  \item{grouppara}{
A postive integer number indicating the number of group specific parameters. For example, if there is no group specific parameter, then grouppara = 0 (default value). If there is only a group specific mean parameter, then grouppara = 1. If there is a group specific intercept and regression coefficient, then grouppara = 2. If data are sampled from one population, that is, there is only one group, then always set grouppara = 0.
}

  \item{jointpara}{
A positive integer number indicating the number of joint parameters. For example, if there is one joint regression coefficients among groups, then jointpara = 1.

}

  \item{n}{
A positive integer number or vector indicating the sample size in the data. If it is a number, then it assumes that data come from one population. If it is a vector, then it assumes that data come from multiple populations, and each element in the vector denotes the group size for the corresponding parameters specified. See Example 2.
}

  \item{ERr}{
A matrix representing equality constraints in a hypothesis. If no equality constraints, set ERr = NULL or ERr = matrix(0,0,0). See details and examples.
}
  \item{IRr}{
A matrix representing inequality or about equality constraints in a hypothesis. If no inequality and about equality constraints, set IRr = NULL or IRr=matrix(0,0,0). See details and examples.
}
  \item{\dots}{
Sets of matrices representing the equality and/or inequality constraints for the extra informative hypotheses. For example, if there are two informative hypotheses, then ... will be two matrices for the second hypothesis: ERr2 and IRr2. See Example 2.
}

  \item{seed}{
The seed for randomly sampling used in Gibbs sampler. Setting different seeds will render slightly different results. Default seed is 100.
}
  \item{print}{
If print is TRUE, then the result will be printed in a table format.
}


}
\details{
The vector "estimate" and the matrix (or list) "Sigma" should be obtained before evaluating hypotheses in Bain. For examples, see Gu, Mulder, and Hoijtink (2017) and Hoijtink, Gu, and Mulder (unpublished) and under the MGBain button at http://informative-hypotheses.sites.uu.nl/ were also the papers and other information can be found.

The length of "estimate" should be equal to grouppara * P + jointpara, where P denotes the number of groups or the length of list "Sigma".

"ERr" is used to specify equality constraints, and "IRr" is used to specify inequality constraints. The general form of the constraints is ER * t^T = r and IR * t^T = r where t is a vector containing the estimates, ER is a length of t by number of equality constraints matrix, IR is a length of t by number of inequality constraints matrix, and r (used with ER) has a length equal to the number of equality constraints, and r (used with IR) has a length equal to the number of inequality constraints. ER, IR, and r contain real numbers.

%%For example, an equality constraint t1=t2+1 can be specified as ER*(t1,t2)^T = r, where ER = (1,-1), and r = 1. %%Then ERr = (ER,r) = (1,-1,1), which is provided to Bain using the command ERr<-matrix(c(1,-1,1),nrow=1,ncol=3%%,byrow = TRUE).


%%Another example, if the inequality constraints are t1 > t2 and t2 > t3 then IR * (t1,t2,t3)^T = r where IR= \cr
%%1 -1 0 \cr
%%0 \sspace 1 -1 \cr
%%and r = (0 0).
%%Then IRr =\cr
%%1 -1  0 0 \cr
%%0  1 -1 0 \cr
%%which is provided to Bain using the command:\cr
%%IRr<-matrix(c(1,-1,0,0,0,1,-1,0),nrow=2,ncol=4,byrow = TRUE)
}

\value{
Returns a list that contains hypothesis testing results (i.e., Bayes factors (BFs), relative fit (f), relative complexity (c), posterior model probabilities (PMPs)), the approximated posterior and prior covariance matrices, the approximated posteriors means, and the fraction (b).

}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%\item{$fit}{the relative fit of each informative hypothesis compared to the unconstrained hypothesis}
%%\item{$complexity}{the relative complexity of each informative hypothesis compared to the unconstrained hypothesis}
%%\item{$BF}{the Bayes factor for each informative hypothesis against its complement}
%%\item{$PMPa}{the posterior model probability of each hypothesis, when considering all the completing hypotheses excluding the unconstrained hypothesis and assuming equal prior model probabilities for all the hypotheses}
%%\item{$PMPb}{the posterior model probability of each hypothesis, when considering all the completing hypotheses including the unconstrained hypothesis and assuming equal prior model probabilities for all the hypotheses}

%%Note that if there are two or more informative hypotheses under consideration, the fit, complexity, and BF will be represented for each informative hypothesis.

%% ...
\references{
Gu, X., Mulder, J., and Hoijtink, H. (2017). Approximated adjusted fractional Bayes
factors: A general method for testing informative hypotheses. British Journal of Mathematical and Statistical Psychology.

Hoijtink, H., Gu, X., and Mulder, J. (unpublished). Multiple group Bayesian evaluation of informative hypotheses.
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
#One group:
#Example 1:
#Hypothesis
#H1: theta1>theta2>theta3   #coefficients in regression model

#Input
estimate<-c(3,2,1)  #estimates of coefficients
Sigma<-matrix(c(3,0,0,0,2,0,0,0,1),3,3,byrow = TRUE) #covariance matrix of coefficients

n<-50 #samples size

#H1
ERr<-NULL
IRr<-matrix(c(1,-1,0,0,0,1,-1,0),nrow=2,ncol=4,byrow = TRUE)

res<-Bain(estimate=estimate,grouppara=0,jointpara=3,Sigma=Sigma,n=n,ERr,IRr) #run
#Results are printed.
#Results for fit, complexity, Bayes factor, and PMPs are saved in "res":

plot(res)
#Results for PMPs are plotted.


#Multiple groups
#Example 2
#t test:
#Hypotheses:
#H1: theta1=theta2   #group means
#H2: theta1>theta2
#H3: theta1<theta2

#Input
estimate<-c(0,0)         #Estimates of group means theta1 and theta2.

cov1<-matrix(c(.5),1,1)
cov2<-matrix(c(.1),1,1)
Sigma<-list(cov1,cov2) #List of variances of group means

n<-c(22,37)           #samplesize

#H1:
ERr1<-matrix(c(1,-1,0),nrow=1,ncol=3,byrow = TRUE)
IRr1<-NULL

#H2
ERr2<-NULL
IRr2<-matrix(c(1,-1,0),nrow=1,ncol=3,byrow = TRUE)

#H3
ERr3<-NULL
IRr3<-matrix(c(-1,1,0),nrow=1,ncol=3,byrow = TRUE)

res<-Bain(estimate,Sigma,grouppara=1,jointpara=0,n=n,ERr1,IRr1,ERr2,IRr2,ERr3,IRr3) #run
#Results are printed.
#Results for fit, complexity, Bayes factor, and PMPs are also saved in "res":

plot(res)
#Results for PMPs are plotted.

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\alias{Bain}



