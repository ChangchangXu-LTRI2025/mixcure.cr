################################
### mixcure.penal.polar.plcr ###
################################


mixcure.penal.polar.plcr <- function(formula, data, apct = 0.05,init, pl){
  
  require(splines)
  require(survival)
  require(abind)
  
  mat.inv <- function(matx) {
    
    detm = det(matx)
    #2x2 matrix inverse;
    if (ncol(matx) == 2) {
      inv.matx = (1/detm) * matrix(c(matx[2,2],-matx[2,1],-matx[1,2],matx[1,1]), nrow = 2)
    }
    
    else {
      #For any n>2 dimension square matrix;
      adjug.matx <- matrix(rep(0, ncol(matx)^2), nrow = nrow(matx))
      for (i in 1:nrow(matx)) {
        for (j in 1:ncol(matx)) {
          adjug.matx[i,j] <- (-1)^(i+j)*det(matx[-i,][,-j])
        }
      }
      inv.matx <- t(adjug.matx/detm)
    }
    
    return(inv.matx)
  }
  
  
  design.matrix <- model.frame(formula, data = data, na.action = na.omit);
  survt <- design.matrix[,1];
  
  design.matrix <- model.matrix(formula, data = design.matrix);
  
  # index ranges of coefficients of glm and cox models
  index.cure.v <- 1 : ncol(design.matrix);
  index.surv.v <- (ncol(design.matrix) + 1) : (2*length(index.cure.v))
  # index of alpha,the shape parameter
  index.gamma <- 2*length(index.cure.v)+1;
  
  samp.s <- nrow(design.matrix)
  
  # p is a vector of values for parameter vector r (radius) for all covariates X and gamma (shape)
  
  loglik.mixture <- function(p, survt, design.matrix, 
                             index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl) {
    
    ####  parameter and variable dep parameters;
    #####
    theta = 1/(1+exp(-design.matrix%*%p[index.cure.var]))
    eps = survt[,1]^(p[index.gamma])*exp(design.matrix%*%p[index.surv.var])
    eta = 1/((exp(eps)-1)*theta+1)
    delta = 1/(theta/(1-theta)*exp(eps)+1)
    #kap = theta*(1-theta)*(1-eta)-(1-theta)^2*eta*(1-eta)
    kap= (1-eta)*(1-theta)*(theta + eta)
    pic = exp(eps)*eps*eta^2
    lambda = (1-theta)^2*eta*(1-eta)*((2*eta-1)*(1-theta)+3)
    phi = theta*(1-theta)*((2*eta-1)*(1-theta)+theta)*pic
    
    
    #calculate loglikelihood for the unpenalized;
    cure.par <- p[1 : ncol(design.matrix) ];
    surv.par <- p[ (ncol(design.matrix) + 1) : (2*length(cure.par)) ];
    p.gamma <- p[ 2*length(cure.par) + 1 ];  #use original shape parameter instead of exp();
    
    # loglikelihood is defined as the negative of the actual loglikelihood for feeding nlm() minimizer;
    loglikelihood <- -sum( ( log(1-theta) + log(p.gamma)-log(survt[,1])
                             +log(eps)-eps )[survt[, 2] == 1] ) -
      sum( (log(theta + (1-theta)*exp(-eps)))[survt[, 2] == 0] );
    
    if (pl == F) {
      loglik = loglikelihood
    } else {
      
      
      n.elema = length(index.cure.var)^2
      a.sub1 <- matrix(rep(0,n.elema), nrow = length(index.cure.var))
      a.sub2 <- matrix(rep(0,n.elema), nrow = length(index.cure.var))
      
      for (i in c(index.cure.var)) {
        for (j in c(index.cure.var)) {
          a.sub1[i,j] <- sum((design.matrix[,i]*design.matrix[,j]*theta*(1-theta))[survt[, 2] == 1])
          a.sub2[i,j] <- sum((design.matrix[,i]*design.matrix[,j]*kap)[survt[, 2] == 0])
        }
      }
      info.a = a.sub1 + a.sub2
      
      
      design.xt <- cbind(design.matrix, log(survt[,1]))
      n.elemb <- length(index.cure.var)*(length(index.cure.var)+1)
      b.sub <- matrix(rep(0,n.elemb), nrow = length(index.surv.var))
      
      for (i in c(index.cure.var)) {
        for (j in c(index.cure.var,length(index.surv.var)+1)) {
          b.sub[i,j] <- -sum((design.matrix[,i]*design.xt[,j]*eps*(1-delta)*delta)[survt[, 2] == 0]) #alternative expression for est
        }
      }
      info.b = b.sub  #Upper right block of fisher.info;
      
      
      n.elemd <- (length(index.surv.var)+1)^2
      d.sub1 <- matrix(rep(0,n.elemd), nrow = (length(index.surv.var)+1))
      d.sub2 <- matrix(rep(0,n.elemd), nrow = (length(index.surv.var)+1))
      for (i in c(index.cure.var,length(index.surv.var)+1)) {
        for (j in c(index.cure.var,length(index.surv.var)+1)) {
          d.sub1[i,j] <- sum((design.xt[,i]*design.xt[,j]*eps)[survt[, 2] == 1])
          d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
          # d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta^2))[survt[, 2] == 0])
          
        }
      }
      info.d = d.sub1 + d.sub2 +
        matrix(c(rep(0, (n.elemd-1)),sum(survt[, 2] == 1)/(p[index.gamma]^2)),nrow = (length(index.surv.var)+1))
      
      
      
      info.d.inv = mat.inv(info.d)
      
      fisher.info = rbind(cbind(info.a,info.b),cbind(t(info.b),info.d))
      #hessian.mat = -fisher.info
      
      # #info.set0 is (A-BD^-1B^T), dif than used in modified score;
      info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)
      
      #determinant of hessian matrix;
      det.info = det(info.set0)*det(info.d)
      #det.info = matrix.det(fisher.info)
      
      
      loglik = loglikelihood - 0.5*log(det.info)
    }
    #loglik = loglikelihood
    return(loglik)
    
  }
  
  
  maximizer0 <- nlm(
    f = loglik.mixture, p = init, survt=survt, design.matrix=design.matrix,
    pl = pl,
    iterlim = 100, hessian=T);
  
  hessmat <- maximizer0$hessian

  var.mat <- solve(hessmat)
  
  alpha.hat <- maximizer0$estimate[index.gamma];
  
  loglik <- -maximizer0$minimum  #in loglik function loglik was calculated as minus of actual loglik value
  
  
  loglik.mixture.profile <- function(p, survt, k=k, design.matrix,  
                                     ka.hat, kb.hat, r.est, phi.est, 
                                     index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl) {
    
    
    ik = k-length(index.cure.var);
    
    #parameter and variable dep parameters;
    if (length(index.cure.var)<3) {
      theta = 1/(1+exp(-design.matrix[,index.cure.var[-k]]*as.matrix(p[index.cure.var[-length(index.cure.var)]])-design.matrix[,k]*(ka.hat+r.est*cos(phi.est))))
      eps = survt[,1]^(p[index.gamma-2])*exp(design.matrix[,index.cure.var[-k]]*as.matrix(p[(index.surv.var[-length(index.cure.var)]-1)])+design.matrix[,k]*(kb.hat+r.est*sin(phi.est)))
    } else {
      theta = 1/(1+exp(-design.matrix[,index.cure.var[-k]]%*%as.matrix(p[index.cure.var[-length(index.cure.var)]])-design.matrix[,k]*(ka.hat+r.est*cos(phi.est))))
      eps = survt[,1]^(p[index.gamma-2])*exp(design.matrix[,index.cure.var[-k]]%*%as.matrix(p[(index.surv.var[-length(index.cure.var)]-1)])+design.matrix[,k]*(kb.hat+r.est*sin(phi.est)))
    }
    eta = 1/((exp(eps)-1)*theta+1)
    delta = 1/(theta/(1-theta)*exp(eps)+1)
    #kap = theta*(1-theta)*(1-eta)-(1-theta)^2*eta*(1-eta) #for LRT
    kap= (1-eta)*(1-theta)*(theta + eta)  # for est,PLCI
    pic = exp(eps)*eps*eta^2
    lambda = (1-theta)^2*eta*(1-eta)*((2*eta-1)*(1-theta)+3)
    phi = theta*(1-theta)*((2*eta-1)*(1-theta)+theta)*pic
    
    ####################################################################################################
    # Note: below constructs fisher info matrix; steps are divide into 4 blocks, 2 square blocks (A&D) #
    # on upper left and lower right, 2 identical transposed blocks (B) on upper right and lower left;  #
    # the idential B blocks are not identical in reduced models unless it's a global LRT, needs to be C#
    ####################################################################################################
    
    #calculate loglikelihood for the unpenalized;
    p.gamma <- p[index.gamma-2];  #use original shape parameter instead of exp();
    
    # loglikelihood is defined as the negative of the actual loglikelihood for feeding nlm() minimizer;
    loglikelihood <- -sum( ( log(1-theta) + log(p.gamma)-log(survt[,1])
                             +log(eps)-eps )[survt[, 2] == 1] ) -
      sum( (log(theta + (1-theta)*exp(-eps)))[survt[, 2] == 0] );
    
    if (pl==T) {
      max.len = max(length(index.cure.var),length(index.surv.var))
      n.elema = max.len^2
      a.sub1 <- matrix(rep(0,n.elema), nrow = max.len)
      a.sub2 <- matrix(rep(0,n.elema), nrow = max.len)
      
      for (i in c(index.cure.var)) {
        for (j in c(index.cure.var)) {
          a.sub1[i,j] <- sum((as.matrix(design.matrix)[,i]*as.matrix(design.matrix)[,j]*theta*(1-theta))[survt[, 2] == 1])
          a.sub2[i,j] <- sum((as.matrix(design.matrix)[,i]*as.matrix(design.matrix)[,j]*kap)[survt[, 2] == 0])
        }
      }
      # if (k <= length(index.cure.var)) {info.a = (a.sub1 + a.sub2)[index.cure.var[-k],index.cure.var[-k]]} else
      info.a = (a.sub1 + a.sub2)[index.cure.var,index.cure.var]
      
      ##info matrix block B
      design.xt <- cbind(design.matrix, log(survt[,1]))
      n.elemb <- max.len*(max.len+1)
      b.sub <- matrix(rep(0,n.elemb), nrow = max.len)
      
      for (i in c(index.cure.var)) {
        for (j in c(1:length(index.surv.var), max.len+1)) {
          b.sub[i,j] <- -sum((as.matrix(design.matrix)[,i]*design.xt[,j]*eps*(1-delta)*delta)[survt[, 2] == 0])
          #b.sub[i,j] <- -sum((design.matrix1[,i]*design.xt0[,j]*eps*(1-eta)*eta*(1-theta))[survt[, 2] == 0])
          
        }
      }
      # if (k <= length(index.cure.var)) {info.b = b.sub[index.cure.var[-k],c(index.surv.var-max.len,index.gamma-max.len)]} else
      # {info.b = b.sub[index.cure.var,c(index.surv.var[-ik]-max.len,index.gamma-max.len)]}
      info.b = b.sub[index.cure.var,c(index.surv.var-max.len,index.gamma-max.len)]
      
      ###info matrix block d
     # design.xt1 <- cbind(design.matrix1, log(survt[,1]))
      
      n.elemd <- (max.len+1)^2
      d.sub1 <- matrix(rep(0,n.elemd), nrow = (max.len+1))
      d.sub2 <- matrix(rep(0,n.elemd), nrow = (max.len+1))
      
      for (i in c(index.surv.var-max.len, max.len +1)) {
        for (j in c(index.surv.var-max.len, max.len +1)) {
          d.sub1[i,j] <- sum((design.xt[,i]*design.xt[,j]*eps)[survt[, 2] == 1])
          d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
          #d.sub2[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*(eps*delta^2))[survt[, 2] == 0])
        }
      }
      d.sub = d.sub1 + d.sub2 +
        matrix(c(rep(0, (n.elemd - 1)),sum(survt[, 2] == 1)/(p[index.gamma-2]^2)),
               nrow = (max.len + 1))
      
      # if (k <= length(index.cure.var))
      #   {info.d = d.sub[c(index.surv.var-max.len,index.gamma-max.len),c(index.surv.var-max.len,index.gamma-max.len)]} else
      #   {info.d = d.sub[c(index.surv.var[-ik]-max.len,index.gamma-max.len),c(index.surv.var[-ik]-max.len,index.gamma-max.len)]}
      
      info.d = d.sub[c(index.surv.var-max.len,index.gamma-max.len),c(index.surv.var-max.len,index.gamma-max.len)]
      
      info.d.inv = mat.inv(info.d)
      
      fisher.info = rbind(cbind(info.a,info.b),cbind(t(info.b),info.d))
      #hessian.mat = -fisher.info
      
      # #info.set0 is (A-BD^-1B^T), dif than used in modified score;
      info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)
      
      #determinant of hessian matrix;
      det.info = det(info.set0)*det(info.d)
      #det.info = matrix.det(fisher.info)
      
      loglik.part = loglikelihood - 0.5*log(det.info)
    } else
      if (pl == FALSE)
      {
        loglik.part = loglikelihood
      }
    
    return(loglik.part)
  }
  
  l.null = loglik - 0.5 * qchisq(1-apct,df=2,ncp = 0,lower.tail=T)
  
  phi_value = seq(0, pi*2, length=100)
  
  ni = 1
  tol = 0.4
  r_phi <- matrix(0, nrow = 100, ncol = (length(index.cure.v)+1))
  
  #ncores=detectCores()-1
 # cl<-makeCluster(10) #change the 2 to your number of CPU cores
 # registerDoParallel(cl)
  
  #r_phi <-rep(0, (length(index.cure.v)+1))    
  
  #  for(i.phi in 1:100){
  for (i.phi in 1:100){
    require(R.utils)
    
      for (k in index.cure.v[-1]) {
        ik = k+length(index.cure.v);
        max.est = maximizer0$estimate 
        n = ni + 1
        r.est=0.25
        param.est.up = max.est 
        sign.delta =1
        
      converge <- FALSE; iter1 <- 1; EXIT1 <-FALSE; delta.up = 0
      while (!converge & iter1 <= 25 & !EXIT1 & !is.nan(delta.up)) {
        
        maximizer.temp <-  nlm(
          f = loglik.mixture, p = max.est, design.matrix=design.matrix,
          survt=survt, pl=pl,
          iterlim = 60, hessian=TRUE)
        score.temp = maximizer.temp$gradient
        hessian.temp = maximizer.temp$hessian
        if (det(hessian.temp) < 1e-05) diag(hessian.temp) <- diag(hessian.temp) + 1e-06
        
        l0.b.up <- -maximizer.temp$minimum
        inv.hessian.temp <- solve(hessian.temp)
        
        
        lambda <- (2*(l0.b.up - l.null + score.temp %*% inv.hessian.temp %*% score.temp)/inv.hessian.temp[k,k])^0.5
        #define increment for estimated value: delta=-A^-1(U-lambda*e)
        delta.up <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda)        
        
        # maximizing loop for unpenalized estimates;
        #if (pl == F) {
        inside <- FALSE; iter2 <- 1;
        while (!inside & iter2 <= 100 & !is.nan(delta.up)) {
          
          # add increment to stepwise parameter value;
          param.est.temp.up <- param.est.up[-c(k,ik)]
          r.est <- r.est + delta.up *sign.delta
          # if (k==2) {param.est.temp.up[1] <- -1;param.est.temp.up[9] <- 0.1}
          
          #compute loglikelihood function using updated parameter values;
          
          
          maximizer.temp1 <- nlm( f = loglik.mixture.profile, p = param.est.temp.up,
                                  survt=survt, k = k, design.matrix=design.matrix, 
                                  ka.hat=max.est[k], kb.hat=max.est[ik], r.est=r.est, phi.est=phi_value[i.phi], 
                                  pl = pl, iterlim = 100, hessian=F)
          
          l.temp.up = -maximizer.temp1$minimum
          
          #if (!is.nan(l.temp.up))
          
          #compare to see if updated l is still
          inside <- (l.temp.up < (l.null - 0.05)) #l.null - 0.05 for all others, 0.2 for k=3 of high rate H0
          alevel.up <- pchisq(2*(loglik - l.temp.up), df=2, ncp=0,lower.tail = T)
          #print(c(delta.up, alevel.up, n,l.temp.up,k,iter1,iter2, i.phi))
          if (!inside) {delta.up <- delta.up/((n+4)/n); sign.delta=1; iter2 <- iter2 + 1}  #(n+0.1)/n for low rate H0;
          # if (is.nan(delta.up)) {param.est.temp.up[k] <- NA}
        } #for iter2
        
        #}
        #Using converged increment for parameter to get corresponding score and variance expressions;
       
        param.est.up <- insert(maximizer.temp1$estimate, ats=k, values=r.est[,1]*cos(phi_value[i.phi]))
        param.est.up <- insert(param.est.up, ats=ik, values=r.est[,1]*sin(phi_value[i.phi]))
        
        l0.b.up = l.temp.up
        
        diff.up = l0.b.up - l.null
        converge <- (abs(diff.up) <= tol)
        if ((!converge| is.nan(l0.b.up)) & !is.nan(delta.up)) {iter1 <- iter1 + 1; n = n + 1; sign.delta=-1} else {EXIT1 = T}
        if (is.nan(delta.up)==T) {param.est.up[k] <- NA}
      } #for iter1
      
      r_phi[i.phi, k] <- r.est
     # print(c(i.phi, k))
      }
      r_phi[i.phi, 1] <- i.phi
      r_phi[i.phi, (length(index.cure.v)+1)] <- phi_value[i.phi]
  #    print(r_phi)
  #    return(r_phi)
  }
  
 # r_phi[, (length(index.cure.v)+1)] <- phi_value
  colnames(r_phi) <- c(colnames(design.matrix), "phi")
  
  z.score <- maximizer0$estimate / sqrt(diag(var.mat));
  
  coef.table.cure <- cbind(
    'coef'        = maximizer0$estimate[index.cure.v],
    'exp(coef)'   = exp(maximizer0$estimate[index.cure.v]),
    'se(coef)'    = sqrt(diag(var.mat)[index.cure.v]),
    'z'           = z.score[index.cure.v],
    'Pr(>|z|)'    = 2 * (1 - pnorm(abs(z.score[index.cure.v]))),
    'LCI.95%' = maximizer0$estimate[index.cure.v] - 1.96 * sqrt(diag(var.mat)[index.cure.v]),
    'UCI.95%' = maximizer0$estimate[index.cure.v] + 1.96 * sqrt(diag(var.mat)[index.cure.v])
  );
  rownames(coef.table.cure) <- colnames(design.matrix);
  
  coef.table.surv <- cbind(
    'coef'        = maximizer0$estimate[index.surv.v],
    'exp(coef)'   = exp(maximizer0$estimate[index.surv.v]),
    'se(coef)'    = sqrt(diag(var.mat)[index.surv.v]),
    'z'           = z.score[index.surv.v],
    'Pr(>|z|)'    = 2 * (1 - pnorm(abs(z.score[index.surv.v]))),
    'LCI.95%' = maximizer0$estimate[index.surv.v] - 1.96 * sqrt(diag(var.mat)[index.surv.v]),
    'UCI.95%' = maximizer0$estimate[index.surv.v] + 1.96 * sqrt(diag(var.mat)[index.surv.v])
  );
  rownames(coef.table.surv) <- colnames(design.matrix);
  
  coef.table.alpha <- cbind(
    'coef'     = alpha.hat,
    'se(coef)' = sqrt(diag(var.mat)[index.gamma]),
    'z'        = z.score[index.gamma],
    'Pr(>|z|)' = 2 * (1 - pnorm(abs(z.score[index.gamma]))),
    'LCI.95%'  = maximizer0$estimate[index.gamma] - 1.96 * sqrt(diag(var.mat)[index.gamma]),
    'UCI.95%'  = maximizer0$estimate[index.gamma] + 1.96 * sqrt(diag(var.mat)[index.gamma]),
    'loglik' = -maximizer0$minimum 
  );
  rownames(coef.table.alpha) <- 'alpha';
  
  out <- list(
    coefficients = list(
      cure = coef.table.cure, 
      surv = coef.table.surv, 
      alpha = coef.table.alpha
      #   run.time
    ),
    plcr = as.data.frame(r_phi),
    cov = var.mat
  );

  return(out)
}
