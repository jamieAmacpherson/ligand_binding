#' Determine the error of the Kd over a range of protein concentrations.
#'
#' This function loops over a range of protein concentrations and determine
#' the error of the binding fit for a defined theoretical Kd.
#' 
#' This function returns a plot of the kd error vs. the protein concentration.
#' 
#'
#' @param P0_min Minimum total protein concentration.
#' @param P0_max Maximum total protein concentration.
#' @param L0_min Minimum total ligand concentration.
#' @param L0_max Maximum total ligand concentration.
#' @param kd Theoretical Kd 
#' @return plt A plot of the fit error vs. P0
#' @export 

simulate_error = function(P0_min, P0_max, L0_min, L0_max, kd){
	
	## initialise a vector containing the total protein concentrations
	prot_conc = seq(from = P0_min, to = P0_max)

	## initialise empty vector to contain error
	errors = c()

	success_prot_conc = c()

	## loop through protein concentrations
	for(i in prot_conc){

		## simulate binding for a given P0, ligand concentration range and theoretical kd
		
		tmp = try(
			  simtitr(kd, i, 12, L0_min, i*5)[[1]]$coefficients[1],
			  silent = TRUE)
		
		## if fitting is successful then calculate the percentage error
		if(class(tmp) == "numeric"){
			
			## calculate the percentage error
			p.error = abs(tmp - 1) * 100
	
			## append the error to errors vector
			errors = append(errors, p.error)

			## append protein conc of successful fit 
			success_prot_conc = append(success_prot_conc, i)
		}
	}
	
	## return 
	dat = as.data.frame(cbind(success_prot_conc, errors))
	return(dat)
}



plt_error = function(dat){

	p0 = seq(from = 1, to = nrow(dat))

	## linear model for data
	lmdat = lm(dat[,2] ~ dat[,1])	

	pltdat = as.data.frame(cbind(p0,
				     lmdat$fitted.values*1.5,
				     lmdat$residuals))

	names(pltdat) = c('p0', 'means', 'sds')

	## extract error of the linear model
	err = summary(lmdat)$coefficients[ , 2][2]*2

	## plot data
	plt = ggplot2::ggplot(pltdat, aes(x=p0, y=means)) + 
		geom_ribbon(aes(ymin = means - means*err, ymax = means + means*err),
			    fill = "grey70", alpha = 0.5) + 
		geom_line(aes(y = means)) +
		theme_bw() +
		theme(text = element_text(size = 26)) +
		xlab('Protein concentration') + 
		ylab('Error')
	
	print(plt)
}

#	par(mar = c(5,5,2,2))

#	plot(dat,
#	     xlab = '[Protein]',
#	     ylab = 'Error',
#	     panel.first=grid(),
#	     cex = 1,
#	     cex.lab = 2,
#	     cex.axis = 2)
#
#	abline(lm(errors ~ success_prot_conc), col = 'red', lwd=2)
#}
