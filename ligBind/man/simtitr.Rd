% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_binding.R
\name{simtitr}
\alias{simtitr}
\title{Simulate ligand binding data.
 
This function simulates protein-ligand binding data given}
\usage{
simtitr(kd, prot_conc, npoints, Lmin, Lmax)
}
\arguments{
\item{kd}{Theoretical 'true' binding affinity between the protein and the ligand}

\item{prot_conc}{The concentration of the free receptor in the titration}

\item{npoints}{The number of points in the ligand titration.}

\item{minimum}{ligand concentration}

\item{maximum}{ligand concentration}
}
\value{
fit.nls.summary List containing a summary of the fitted data.
}
\description{
The output of the function is a list containing the fitted data as well as the estimated binding affinity.
}
