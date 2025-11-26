#######################################
### Function for mixcure model under ##
#### penalized loglikelihoods or ML  ##
#######################################################################
#### PROFILE LIKELIHOOD CONFIDENCE REGION FOR BIVARIATE PARAMETERS ####
#### REGION ESTIMATION FOR CURE PART OF MC                         ####
#######################################################################

########### CREATED ON 2021-10-26;
########### LAST MODIFIED 2021-10-26:

mixcure.penal.lrt.2d <- function(formula, data, init, pl, est.tru, iterlim = 200) {
require(splines)
require(survival)
require(abind)

  #########################################################################################
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
  loglik.mixture <- function(p, survt, design.matrix, index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl, PLCI=F) {
    
    ####  parameter and variable dep parameters;
    #####
    theta = 1/(1+exp(-design.matrix%*%p[index.cure.var]))
    eps = survt[,1]^(p[index.gamma])*exp(design.matrix%*%p[index.surv.var])
    eta = 1/((exp(eps)-1)*theta+1)
    delta = 1/(theta/(1-theta)*exp(eps)+1)
    # kap = theta*(1-theta)*(1-eta)-(1-theta)^2*eta*(1-eta)   ##for mixcure.penal.mi.R
    kap= (1-eta)*(1-theta)*(theta + eta)
    pi = exp(eps)*eps*eta^2
    lambda = (1-theta)^2*eta*(1-eta)*((2*eta-1)*(1-theta)+3)
    phi = theta*(1-theta)*((2*eta-1)*(1-theta)+theta)*pi
    
    #calculate loglikelihood for the unpenalized;
    cure.par <- p[1 : ncol(design.matrix) ];
    surv.par <- p[ (ncol(design.matrix) + 1) : (2*length(cure.par)) ];
    p.gamma <- p[ 2*length(cure.par) + 1 ];  #use original shape parameter instead of exp();
    
    # loglikelihood is defined as the negative of the actual loglikelihood for feeding nlm() minimizer; 
    loglikelihood <- -sum( ( log(1-theta) + log(p.gamma)-log(survt[,1])
                             +log(eps)-eps )[survt[, 2] == 1] ) - 
      sum( (log(theta + (1-theta)*exp(-eps)))[survt[, 2] == 0] );
    
    if (pl==T) {
      ####calculate inverse of info matrix by block matrix;
      
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
          b.sub[i,j] <- -sum((design.matrix[,i]*design.xt[,j]*theta*(1-theta)*pi)[survt[, 2] == 0])
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
          
        }
      }
      info.d = d.sub1 + d.sub2 + 
        matrix(c(rep(0, (n.elemd-1)),sum(survt[, 2] == 1)/(p[index.gamma]^2)),nrow = (length(index.surv.var)+1))
      
      
      info.d.inv = mat.inv(info.d)
      
      #    fisher.info = rbind(cbind(info.a,info.b),cbind(t(info.b),info.d))
      #hessian.mat = -fisher.info
      
      # #info.set0 is (A-BD^-1B^T), dif than used in modified score;
      info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)
      
      #determinant of hessian matrix;
      det.info = det(info.set0)*det(info.d)
      #   det.info = matrix.det(fisher.info)
      
      if (PLCI==T)  
      {loglik = loglikelihood - 0.5*log(abs(det.info))} else {
        
        loglik = loglikelihood - 0.5*log(det.info)
      }
    }
    
    else if (pl == FALSE)
    {
      loglik = loglikelihood
    }
    
    #loglik = loglikelihood
    return(loglik)
    
  }
  
  ######END of loglik.mixture####################################
  
  
  # Parameter estimation under Ha (non-restricted likelihood)
  # maximize penalized or unpenalized loglikelihood by nlm; 
  maximizer0 <- nlm(
    f = loglik.mixture, p = init, survt=survt, design.matrix=design.matrix, 
    pl = pl, 
    iterlim = iterlim, hessian=TRUE);

  loglik0 <- -maximizer0$minimum
  
  loglik.mixture.profile <- function(p, survt, k.cur, k.sur, param.est.cur, param.est.sur, design.matrix1=design.matrix, design.matrix0=design.matrix, index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl) {
    
    design.mtx.comb = cbind(design.matrix0,design.matrix1)

    theta = 1/(1+exp(-design.matrix1[,index.cure.var[-k.cur]]%*%as.matrix(p[index.cure.var[-length(index.cure.var)]]) - design.mtx.comb[,k.cur]*param.est.cur))
    eps = survt[,1]^(p[index.gamma-2])*exp(design.mtx.comb[,index.surv.var[-k.cur]]%*%as.matrix(p[index.surv.var[-length(index.surv.var)]-1]) + design.mtx.comb[,k.sur]*param.est.sur)
    
    eta = 1/((exp(eps)-1)*theta+1)
    delta = 1/(theta/(1-theta)*exp(eps)+1)
    kap= (1-eta)*(1-theta)*(theta + eta)  # for est,PLCI
    pi = exp(eps)*eps*eta^2
    lambda = (1-theta)^2*eta*(1-eta)*((2*eta-1)*(1-theta)+3)
    phi = theta*(1-theta)*((2*eta-1)*(1-theta)+theta)*pi
    
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
          a.sub1[i,j] <- sum((as.matrix(design.matrix0)[,i]*as.matrix(design.matrix0)[,j]*theta*(1-theta))[survt[, 2] == 1])
          a.sub2[i,j] <- sum((as.matrix(design.matrix0)[,i]*as.matrix(design.matrix0)[,j]*kap)[survt[, 2] == 0])
        }
      }
      # if (k <= length(index.cure.var)) {info.a = (a.sub1 + a.sub2)[index.cure.var[-k],index.cure.var[-k]]} else
      info.a = (a.sub1 + a.sub2)[index.cure.var,index.cure.var]
      
      ##info matrix block B
      design.xt0 <- cbind(design.matrix0, log(survt[,1]))
      n.elemb <- max.len*(max.len+1)
      b.sub <- matrix(rep(0,n.elemb), nrow = max.len)
      
      for (i in c(index.cure.var)) {
        for (j in c(1:length(index.surv.var), max.len+1)) {
          b.sub[i,j] <- -sum((as.matrix(design.matrix1)[,i]*design.xt0[,j]*eps*(1-delta)*delta)[survt[, 2] == 0])
          #b.sub[i,j] <- -sum((design.matrix1[,i]*design.xt0[,j]*eps*(1-eta)*eta*(1-theta))[survt[, 2] == 0])
          
        }
      }
      
      info.b = b.sub[index.cure.var,c(index.surv.var-max.len,index.gamma-max.len)]
      
      ###info matrix block d
      design.xt1 <- cbind(design.matrix1, log(survt[,1]))
      
      n.elemd <- (max.len+1)^2
      d.sub1 <- matrix(rep(0,n.elemd), nrow = (max.len+1))
      d.sub2 <- matrix(rep(0,n.elemd), nrow = (max.len+1))
      
      for (i in c(index.surv.var-max.len, max.len +1)) {
        for (j in c(index.surv.var-max.len, max.len +1)) {
          d.sub1[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*eps)[survt[, 2] == 1])
          d.sub2[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
          #d.sub2[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*(eps*delta^2))[survt[, 2] == 0])
        }
      }
      d.sub = d.sub1 + d.sub2 +
        matrix(c(rep(0, (n.elemd - 1)),sum(survt[, 2] == 1)/(p[index.gamma-2]^2)),
               nrow = (max.len + 1))
      
      
      info.d = d.sub[c(index.surv.var-max.len,index.gamma-max.len),c(index.surv.var-max.len,index.gamma-max.len)]
      
      info.d.inv = mat.inv(info.d)
      
      fisher.info = rbind(cbind(info.a,info.b),cbind(t(info.b),info.d))
      
      # #info.set0 is (A-BD^-1B^T), dif than used in modified score;
      info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)
      
      #determinant of hessian matrix;
      det.info = det(info.set0)*det(info.d)
      #det.info = matrix.det(fisher.info)
      
      loglik.part = loglikelihood - 0.5*log(det.info)
    } else if (pl == FALSE)
    {
      loglik.part = loglikelihood
    }
    
    return(loglik.part)
  }
  

  #################################################################
  #### parameter estimation under H0 for individual parameter
  #### loglikelihood ratio test statistics for parameters of a single variable;

  dim.v <- ncol(design.matrix)

    ll.compl3 <- matrix(0, nrow=dim.v,ncol = 1)
    ll.pval <- matrix(0, nrow=dim.v,ncol = 1)
    
    for (k in index.cure.v[-1]) {

      ik = k + length(index.cure.v)
      maximizer <- nlm(
        f = loglik.mixture.profile, p = init[-c(k,ik)], k.cur = k, k.sur=ik,
        param.est.cur=est.tru[k], param.est.sur=est.tru[ik],
        survt = survt, design.matrix0 = design.matrix,
        design.matrix1=design.matrix,
        pl=pl, iterlim = iterlim, hessian=F);

      loglik.part = -maximizer$minimum;
      ll.compl3[k,1]<- loglik.part
      ll.pval[k,1] <- pchisq(2*abs(loglik0 - loglik.part), df=2, lower.tail = F)
    }

  ll.table <- cbind.data.frame(llr=ll.compl3,pval=ll.pval);
  rownames(ll.table) <- colnames(design.matrix);


#run.time = proc.time() - init.time


#######################################
## Output tables from either method; ##
#######################################


out <- list(ll.table);
class(out) <- c('mixcure', 'list');

return(out);

}

