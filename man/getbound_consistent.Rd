% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation.R
\name{getbound_consistent}
\alias{getbound_consistent}
\title{getbound_consistent}
\usage{
getbound_consistent(obj, constraints, slack, verb = 1)
}
\arguments{
\item{obj}{defines the linear coefficients of the linear objective function}

\item{constraints}{output from make_constraints() - defines the linear constraints in the form of A%*%x - b <= b_n}

\item{slack}{defines the amount to slack b_n in the equality constraints to get consistent plug-in estimator of the identified set.}

\item{verb}{defines how much information of the optimization routine to be printed. Default is set at 1 for no print, change to 5 for full output print}
}
\value{
bound is the lower and upper bound of the value function, phat is the optimized solution for lower and upper bound, status records the status of the two optimization problem
}
\description{
calls Rmosek to solve for two linear programming problem for the identified set with a small slack introduced on the constraints for consistent plug-in estimator of the identified set
}
\examples{
 ## min/max_x obj'x
 s.t. Acons \%*\% x - bcons <= slack
      Acons \%*\% x - bcons >= -slack
}
