\name{plot.Bain}
\alias{plot.Bain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Pie chart for Bain result.
}
\description{
  This function plots pie charts for posterior model probabilities of candidate hypotheses.
}
\usage{
\method{plot}{Bain}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Result from Bain, Bain_ttest, Bain_ttestData, Bain_anova, Bain_anova_core,
    Bain_ancova, Bain_ancova_core or Bain_regression function.
  }
  \item{...}{
    Optional arguments passed to the plot function
  }
}

\value{
Plot a pie chart of posterior model probabilities if only one informative hypothesis is evaluated.
Plot two pie charts of posterior model probabilities if more than one informative hypotheses are evaluated: one for posterior model probabilities excluding the unconstrained hypothesis; one for posterior model probabilities including the unconstrained hypothesis.

}

  \seealso{
    \code{\link{Bain},\link{plot}}
  }

\alias{plot.Bain}



