% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/surveyIdxPlots.R
\name{factorplot}
\alias{factorplot}
\title{Plot estimates of factor levels (or random effects) from a GAM model.}
\usage{
factorplot(x, name, invlink = exp, levs = NULL, ylim = NULL, ...)
}
\arguments{
\item{x}{A model of class 'gam' or 'glm'}

\item{name}{name of the factor}

\item{invlink}{inverse link function}

\item{levs}{optional custom factor level names}

\item{...}{extra arguments to 'plot'}
}
\description{
Plot estimates of factor levels (or random effects) from a GAM model.
}
