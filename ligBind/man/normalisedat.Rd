% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_binding.R
\name{normalisedat}
\alias{normalisedat}
\title{Fit a binding curve to ligand binding data assuming a 1:1 binding stoichiometry.
 
This function reads a data-frame of the concentration of the ligand vs. the response (eg. fluorescence signal)
of the binding event. The user also supplies a predicted binding affinity and the concentration of the receptor.
used in the experiment.}
\usage{
normalisedat(datavec)
}
\arguments{
\item{bindingdat}{Dataframe with two equal-length columns including 1. the concentration of the titrant and 2. the binding response at that concentration of titrant}

\item{kd_pred}{Initial estimate of the binding affinity}

\item{prot_conc}{Concentration of receptor used in the titration.}
}
\value{
fit.nls.summary List containing a summary of the fitted data.
}
\description{
The output of the function is a list containing the fitted data as well as the estimated binding affinity.
}
