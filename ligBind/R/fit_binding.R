#' Fit a binding curve to ligand binding data assuming a 1:1 binding stoichiometry.
#'  
#' This function reads a data-frame of the concentration of the ligand vs. the response (eg. fluorescence signal)
#' of the binding event. The user also supplies a predicted binding affinity and the concentration of the receptor.
#' used in the experiment.
#'
#' The output of the function is a list containing the fitted data as well as the estimated binding affinity.
#'
#' @param bindingdat Dataframe with two equal-length columns including 1. the concentration of the titrant and 2. the binding response at that concentration of titrant 
#' @param kd_pred Initial estimate of the binding affinity
#' @param prot_conc Concentration of receptor used in the titration.
#' @return fit.nls.summary List containing a summary of the fitted data.
#' @export 



## normalise obs between 0 and 1
normalisedat = function(datavec) {
        tmp = (datavec - min(datavec))/(max(datavec)-min(datavec))
        return(tmp)
}



## Fiting function
fit_binding = function(bindingdat, kd_pred, prot_conc){

	bindingdat = as.data.frame(cbind(bindingdat[,1], normalisedat(bindingdat[,2])))
	
	# rename the dataframe columns
	names(bindingdat) = c("L", "obs")
		
	# Constrain the free protein concentration to the input value 	
	Po = prot_conc

	# Quadratic function used to fit binding curve	
	fiteq = as.formula(obs ~
             	(Sp*Po) +
                ((Spl-Sp) *
                ((((kd.cal+L+Po)-(sqrt((kd.cal + Po + L)^2 - (4 * Po * L))))/2))))

	# Non-linear least squares fitting of above equation
    fit.nls = nls(fiteq,
        data=bindingdat,
        start=list(kd.cal=kd_pred, Spl=1, Sp=0), 
        	
		# Constrain the kd, Spl and Sp
		lower = c(kd_pred * 0.1, 0, 0),   
        upper = c(kd_pred * 10, 1, 1), 
        	
		# Verbose output
		algorithm="port",
        trace=T)
	    
	fit.nls.summary = summary(fit.nls)
	print(fit.nls.summary)
    fit.nls.Kd      = fit.nls.summary$param[1]
    fit.nls.predict = predict(fit.nls)
    results 		= as.data.frame(cbind(bindingdat, fit.nls.predict))
	names(results) 	= c("l", "i", "fit")
	residuals 		= residuals(fit.nls)	

	# predict binding model
	mml = data.frame(L = seq(0, max(bindingdat$L), length.out = 100))
	mml$obs = predict(fit.nls, newdata = mml)
	    
    # Return summary of single-site fitting model
    outlist = list(fit.nls.summary, results, mml, residuals)

    return(outlist)

}


plt_binding_curve = function(fitting_model_list){

	## raw binding data as an element in the list returned by fit_binding()
	results = fitting_model_list[[2]]

	## the binding model is reassigned to a vector
	mml = fitting_model_list[[3]]

	## the residuals between the raw data and the fit are reassigned to a vector
	residuals = fitting_model_list[[4]]

	
	# arrange layout for main plot and subplot of the residuals
	layout(matrix(1:2, ncol=1),
		widths=1,
		heights=c(2.5, 1.5),
		respect=FALSE)

	# set bespoke margins
	par(mar = c(0, 5, 1, 1.2))

	# Plot the raw binding signal vs. ligand concentration
    plot(results$l, results$i, 
        xlab="[Ligand]",
        ylab="Fraction bound",
        cex.lab=1.5,
        cex.axis=1.5,
        cex=1.5,
        panel.first=grid(),
        xaxt = 'n',
        ylim = c(0,1.4))
        	
    lines(mml$L, mml$obs, col="red")

	# subplot of the residuals
	par(mar = c(5, 5, 0, 1.2))
	plot(results$l, residuals,
	        pch=16,
	        cex=1.5,
	        cex.lab=1.5,
	        cex.axis=1.5,
	        xlab="[Ligand]",
	        ylab="Residuals",
	        ylim = c((min(residuals) + (min(residuals)*0.4)), max(residuals) + (max(residuals)*0.4)),
	        panel.first=grid())

	# line going through y=0 
	abline(0,0, lty=2, lwd=1)

}











