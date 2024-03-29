% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSOFunc_eval_js_jac.R
\name{eval_cs_jac}
\alias{eval_cs_jac}
\title{Jacobian (partial derivative matrix) for the inequalities}
\usage{
eval_cs_jac(x, spp)
}
\arguments{
\item{x}{a list/column, e.g. run$solution,}

\item{spp}{a data.frame, usually data or option data.frame}
}
\value{
returns a vector
}
\description{
Jacobian (partial derivative matrix) for the inequalities
}
\details{
Internal GSO function. Used within gso(), an input for nloptr()

Since each constraint is linear, the Jacobian is the coefficients of each
variable in eval_cs().  Calculates the coefficient for each variable in eval_cs()
}
