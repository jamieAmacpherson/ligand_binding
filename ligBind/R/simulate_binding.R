#' Simulate ligand binding data.
#'  
#' This function simulates protein-ligand binding data given 
#'
#' The output of the function is a list containing the fitted data as well as the estimated binding affinity.
#'
#' @param kd Theoretical 'true' binding affinity between the protein and the ligand
#' @param prot_conc The concentration of the free receptor in the titration
#' @param npoints The number of points in the ligand titration.
#' @return fit.nls.summary List containing a summary of the fitted data.
#' @export 


simtitr = function(kd, prot_conc, npoints) {
	
	## initialise spectra signals for protein and protein-ligand complex
	Sp = 1
    Spl = 4

    ## simulated data using the Martin equation
    Po = prot_conc
    Lmin = 0.01
    Lmax=Po*100
    Lo = emdbook::lseq(Lmin, Lmax, npoints)
    kd.pred = kd

    ## gradient of the quadratic
    b = ((((kd.pred+Lo+Po)-(sqrt((kd.pred + Po + Lo)^2 - (4 * Po * Lo))))/2))

    ## calculate the observed spectral signal
    obs = (Sp*Po)+((Spl-Sp)*b)

    ## normalise the observed spectra signal
    normobs = (obs-min(obs))/(max(obs)-min(obs))

    ## noise added to the binding data
    noisey.nobs = normobs + rnorm(length(normobs), sd=0.025)

    titr = as.data.frame(cbind(Lo, noisey.nobs))
        
    ## fit the simulated binding data
    ligBind::fit_binding(bindingdat = titr, kd_pred = kd, prot_conc = prot_conc)
}

