% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sde.R
\name{sde_milstein}
\alias{sde_milstein}
\title{Simulation of gbm}
\usage{
sde_milstein(a, b, x0, t, delta_t)
}
\arguments{
\item{a}{the paramter of gbm}

\item{b}{the paramter of gbm}

\item{x0}{the start point of gbm}

\item{t}{the time of gbm}

\item{delta_t}{the step of simulation}
}
\value{
a path of gbm
}
\description{
Simulation of gbm with milstein method.
}
