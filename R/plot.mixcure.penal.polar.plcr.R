#########################################
### plotting mixcure.penal.polar.plcr ###
#########################################


plot.mixcure.penal.polar.plcr <- function(obj_cr, cr_type="PLCR", title_name =NULL, plci=NULL) {
  
  require(ggplot2)
  require(reshape)
  require(reshape2)
  require(ellipse)
  library(dplyr)
  library(purrr)
  library(tidyr)
  
  plcr.dat <- obj_cr$plcr
  est <- obj_cr$coefficients
  n.parm <- length(rownames(est$cure)[-1])
  vcov.max <- obj_cr$cov[-c(1, (2+n.parm)), -c(1, (2+n.parm))]
  
  dat_r_phi_5var0 <- melt(plcr.dat[,-1], id=c("phi"))
  colnames(dat_r_phi_5var0) <- c("phi", "variable","value")
  
  a.value.tab <- cbind.data.frame(variable = rownames(est$cure)[-1], 
                                  axmle = est$cure[,1][-1])
  b.value.tab <- cbind.data.frame(variable = rownames(est$surv)[-1],
                                  bymle = est$surv[,1][-1])
  
  merged_df <- merge(dat_r_phi_5var0, a.value.tab, by = "variable")
  dat_r_phi_5var <- merge(merged_df, b.value.tab, by = "variable")
  dat_r_phi_5var$ax <- with(dat_r_phi_5var, axmle+value*cos(phi)) 
  dat_r_phi_5var$by <- with(dat_r_phi_5var, bymle+value*sin(phi)) 
  


  if (cr_type=="Wald") {
    
    means <- data.frame(
      group = rownames(est$cure)[-1], ax = est$cure[,1][-1], by = est$surv[,1][-1],
      lo95.ax=est$cure[,6][-1], hi95.ax=est$cure[,7][-1],
      lo95.by=est$surv[,6][-1], hi95.by=est$surv[,7][-1]
    )
    
    covs <- list()
    for (i in 1: n.parm) {
      covs[[i]] = vcov.max[c(i,(n.parm+i)),c(i,(n.parm+i))]
     }
    names(covs) <- rownames(est$cure)[-1]
    
    wald_ellipse_df <- function(mu, V, group, level = 0.95, pts = 181) {
      c2 <- qchisq(level, df = 2)
      R <- chol(V)
      t <- seq(0, 2*pi, length.out = pts)
      xy <- t(sqrt(c2) * t(R) %*% rbind(cos(t), sin(t)))
      tibble(
        x = mu[1] + xy[,1],
        y = mu[2] + xy[,2]
      ) |> mutate(group = group)
    }
    
    ellipses <- means %>%
      mutate(ellipse = pmap(
        list(mu = map2(ax, by, c), group = group),
        \(mu, group) wald_ellipse_df(mu, covs[[group]], group)
      )) %>%
      select(group, ellipse) %>%
      unnest(ellipse, names_sep = "_")
    
 
    ggplot() +
      geom_path(data = ellipses, aes(x = ellipse_x, y = ellipse_y)) +
      geom_point(data = means, aes(x = ax, y = by), size = 2) +
      geom_errorbarh(data = means, aes(y = by, xmin = lo95.ax, xmax = hi95.ax), height = 0.5) +
      geom_errorbar(data = means, aes(x = ax, ymin = lo95.by, ymax = hi95.by), width = 0.5) +
      coord_equal() +
      facet_wrap(~ group) +
      theme_bw()
    
    # wald_data <- cbind.data.frame(ax = est$cure[,1][-1], by= est$surv[,1][-1], lo95.ax=est$cure[,6][-1], hi95.ax=est$cure[,7][-1],
    #                               lo95.by=est$surv[,6][-1], hi95.by=est$surv[,7][-1],
    #                               variable = rownames(est$cure)[-1])
    # wald.cr.dat <- data.frame()
    # for (i in 1: n.parm) {
    #   corr.ML = vcov.max[c(i,(n.parm+i)),c(i,(n.parm+i))]
    #   wald.cr.dat <- rbind.data.frame(cbind(rownames(est$cure)[-1][i], 
    #                                         ellipse::ellipse(x=corr.ML, scale = c(1,1), 
    #                                                 centre = c(est$cure[,1][i+1], est$surv[,1][i+1]))), wald.cr.dat)
    # }
    # colnames(wald.cr.dat) <- c("variable", "ax", "by")
    # 
    # ggplot(wald.cr.dat, aes(x=ax, y=by), cex = 3, pch = 19) +
    #   geom_point() + ggtitle(paste0(title_name)) +
    #   geom_errorbar(wald_data, mapping = aes(ymin = lo95.by, ymax = hi95.by, width=.1)) +
    #   geom_errorbar(wald_data, mapping = aes(xmin = lo95.ax, xmax = hi95.ax, width=.1)) +
    #   geom_point(wald_data, mapping = aes(x=ax, y=by)) +
    # #  geom_point(ref.dat.alt, mapping=aes(x=alpha, y=beta),color="blue", size =2, shape=1) +
    #   facet_wrap(~variable) +
    #   theme(legend.position = "top")
    
  } else {
  
  if (is.null(plci)==TRUE) {
 
  ggplot(dat_r_phi_5var, aes(x=ax, y=by), cex = 3, pch = 19) +
    geom_point() + xlab("log OR") + ylab("log HR") +
    ggtitle(paste0(title_name)) +
    facet_wrap(~variable) +
    theme(axis.text.y=element_text(face="bold", size =10),
          #      axis.ticks.y=element_blank(),
          axis.text.x=element_text(face="bold", size =10),
          axis.title=element_text(size=16,face="bold"),
          legend.text = element_text(size=16),
          legend.title = element_blank(),
          #      strip.text.y.left = element_text(hjust=0,vjust = 1, angle=0,face="bold"),
          strip.text = element_text(size=16),
          legend.position = "top") + 
    geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.3) + 
    geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.3)
  }
    
 else {
    
    plci_data <- cbind.data.frame(ax = plci$coefficients$cure[,1][-1], by= plci$coefficients$surv[,1][-1], lo95.ax=plci$coefficient$cure[,3][-1], 
                                  hi95.ax=plci$coefficient$cure[,4][-1],
                                  lo95.by=plci$coefficient$surv[,3][-1], hi95.by=plci$coefficient$surv[,4][-1],
                                  variable = rownames(plci$coefficient$cure)[-1]) 
    
    ggplot(dat_r_phi_5var, aes(x=ax, y=by), cex = 3, pch = 19) +
      geom_point() + xlab("log OR") + ylab("log HR") +
      ggtitle(paste0(title_name)) +
        geom_errorbar(plci_data, mapping = aes(ymin = lo95.by, ymax = hi95.by, width=.1)) +
        geom_errorbar(plci_data, mapping = aes(xmin = lo95.ax, xmax = hi95.ax, width=.1)) +
        geom_point(plci_data, mapping = aes(x=ax, y=by)) +
      facet_wrap(~variable) +
      theme(axis.text.y=element_text(face="bold", size =10),
            #      axis.ticks.y=element_blank(),
            axis.text.x=element_text(face="bold", size =10),
            axis.title=element_text(size=16,face="bold"),
            legend.text = element_text(size=16),
            legend.title = element_blank(),
            #      strip.text.y.left = element_text(hjust=0,vjust = 1, angle=0,face="bold"),
            strip.text = element_text(size=16),
            legend.position = "top") + 
      geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.3) + 
      geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.3)
    
  }
  }
}  
    
    