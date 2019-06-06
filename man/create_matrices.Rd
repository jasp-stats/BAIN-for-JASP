\name{create_matrices}
\alias{create_matrices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create Bain (in)equality constraint matrices
}
\description{
  Parses a character string describing a set of BAIN informative hypotheses,
  and returns BAIN (in)equality constraint matrices
}
\usage{
create_matrices(varnames, hyp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varnames}{
   A character (vector of characters), containing names of variables used in the hypotheses.
    %Object of class \code{\link{lm}}, from which the model parameters are extracted.
  }

  \item{hyp}{
   A character string, containing a Bain hypothesis (see Details).
  }

}

\details{
Informative hypotheses specified as a character string by "hyp" should adhere
to the following simple syntax:
\itemize{
   \item Competing hypotheses are separated by ";".
         Thus, "a=b;a>b" means that H1: a=b, and H2: a>b.
   \item Each individual hypothesis consists of a (series of) (in)equality
         constraint(s). Every single (in)equality constraint is of the form
         "R1*mu1 + R2*mu2+... = r", where capital Rs refer to numeric scaling
         constants, must refer to the names of parameters in the model, and
         the lower case r refers to a constant. Standard mathematical
         simplification rules apply; thus, "R1*mu1 = R2*mu2" is equivalent to
         "R1*mu1 - R2*mu2 = 0".
   \item Multiple unrelated constraints within one hypothesis can be chained
         by "&". Thus, "a=b&c=d" means that H1: a=b AND c=d.
   \item Multiple related constraints within one hypothesis can be chained by
         repeating the (in)equality operators "=", "<", or ">". Thus, "a<b<c"
         means that H1: a < b AND b < c.
   \item Parameters can be grouped by placing them in a parenthesized, comma
         separated list. Thus, "(a,b)>c" means that H1: a > c AND b > c.
         Similarly, "(a,b)>(c,d)" means that H1: a > c AND b > c AND b > c AND
         b > d.
 }
}

\value{
A pair of named matrices for every hypothesis specified in the
\code{hyp} argument; one matrix named ERr, specifying equality constraints,
and one matrix named IRr, specifying inequality constraints.
}


\author{
  Caspar van Lissa
}
%%\note{
  %%  ~~further notes~~
    %%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

  %%  \seealso{
    %% ~~objects to See Also as \code{\link{help}}, ~~~
      %%  }
\examples{

varnames <- c("a","b","c","d","e","f")
hyp <- "e=f<a=b=c"
create_matrices(varnames, hyp)

hyp <- ".5*f>e>0&c>0&d=0&a=b=0"
create_matrices(varnames, hyp)

hyp <- "f>e<0&c>0&d=0&a=b=0"
create_matrices(varnames, hyp)

hyp <- "f<e>0&c>0&d=0&a=b=0"
create_matrices(varnames, hyp)

hyp <- "a>b<c>d"
create_matrices(varnames, hyp)

hyp <- "d<c>b<a"
create_matrices(varnames, hyp)

hyp <- "a>2"
create_matrices(varnames, hyp)

hyp <- "a>2& b+3=0& b>c=d&e=f=4"
create_matrices(varnames, hyp)

hyp <- "(a,b)>c&d=0"
create_matrices(varnames, hyp)

hyp <- "a>-1&a<1"
create_matrices(varnames, hyp)

hyp <- "a>b&b>c&c>a"
create_matrices(varnames, hyp)

hyp <- "a+b>2"
create_matrices(varnames, hyp)

hyp <- "(a+b)>2"
create_matrices(varnames, hyp)

hyp <- "a&b>c&d"
create_matrices(varnames, hyp)

hyp <- "(.5*a,b)>c"
create_matrices(varnames, hyp)

hyp <- "1/2*a>c"
create_matrices(varnames, hyp)

hyp <- "a < -2"
create_matrices(varnames, hyp)
hyp <- "-2 > a"
create_matrices(varnames, hyp)

hyp <- "a>2; b=0; c< -0.5; d>e=f"
create_matrices(varnames, hyp)

hyp1 <- "a>2; b=0; c< -.5; d>e=f"
create_matrices(varnames, hyp1)

hyp <- "(a+a+a,c+c+c)>b+b"
create_matrices(varnames, hyp)

hyp1 <- "a+a+b>b+b+a"
create_matrices(varnames, hyp1)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


\alias{create_matrices}



