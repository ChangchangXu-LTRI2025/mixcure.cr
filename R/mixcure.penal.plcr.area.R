##############################################
### Mixcure PLCR area by polar             ###
### Input: mixcure.penal.polar.plcr output ###
##############################################

mixcure.penal.plcr.area <- function(obj) {
  require(splines)
  require(survival)
  require(abind)
  library(reshape)
  library(reshape2)


  r.phi.mcvar = obj$plcr
  dat_r_phi_mcvar0 <- melt(r.phi.mcvar[,-1], id=c("phi"))
  colnames(dat_r_phi_mcvar0) <- c("phi", "variable","value")
  varname <- unique(dat_r_phi_mcvar0$variable)

  area.vec <- rep(NA, length(varname))
  for (i in 1: length(varname) ) {

  datafm2 <-dat_r_phi_mcvar0[dat_r_phi_mcvar0$variable==varname[i],]
  radius.val2 <- datafm2$value
  area.vec[i] <- sum(0.0635*radius.val2^2)
  }

  out.area <- cbind.data.frame(
    var.name = varname,
    cr.area = area.vec
  )
  return(out.area)
}
