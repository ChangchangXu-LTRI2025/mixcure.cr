################################
### Function for mixcure model##
#### penalized loglikelihoods ##
###################################################
#### Last modified Dec31 2018 for one x variable ##
###################################################

########### NOTE on Oct22:
########## Error checking for cure uppper points iter2;
########## maximizer$temp1 is not consistent in and out of iter2 loop;

mixcure.penal.profile.point <- function(formula, data, init, k, kval, pl, iterlim = 200) {
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



  loglik.mixture <- function(p, survt, design.matrix, index.cure.v, index.surv.v, pl) {

    ####  parameter and variable dep parameters;
    #####
    p[k] <- kval[1]
    p[k+length(index.cure.v)] <- kval[2]

    theta = 1/(1+exp(-design.matrix%*%p[index.cure.v]))
    eps = survt[,1]^(p[index.gamma])*exp(design.matrix%*%p[index.surv.v])
    #calculate loglikelihood for the unpenalized;
    # cure.par <- p[1 : ncol(design.matrix) ];
    # surv.par <- p[ (ncol(design.matrix) + 1) : (2*length(cure.par)) ];
    # p.gamma <- p[ 2*length(cure.par) + 1 ];  #use original shape parameter instead of exp();
    p.gamma <- p[index.gamma];  #use original shape parameter instead of exp();

    # loglikelihood is defined as the negative of the actual loglikelihood for feeding nlm() minimizer;
    loglikelihood <- -sum( ( log(1-theta) + log(p.gamma)-log(survt[,1])
                             +log(eps)-eps )[survt[, 2] == 1] ) -
      sum( (log(theta + (1-theta)*exp(-eps)))[survt[, 2] == 0] );

    if (pl==T) {

      theta = 1/(1+exp(-design.matrix%*%p[index.cure.v]))
      eps = survt[,1]^(p[index.gamma])*exp(design.matrix%*%p[index.surv.v])
      eta = 1/((exp(eps)-1)*theta+1)
      delta = 1/(theta/(1-theta)*exp(eps)+1)
      #kap = -theta*(1-theta)*(1-eta)+(1-theta)^2*eta*(1-eta)  # for est and PLCI
      kap= (1-eta)*(1-theta)*(theta + eta)    # exp for est and PLCI
      pi = exp(eps)*eps*eta^2

      ####calculate inverse of info matrix by block matrix;

      n.elema = length(index.cure.v)^2
      a.sub1 <- matrix(rep(0,n.elema), nrow = length(index.cure.v))
      a.sub2 <- matrix(rep(0,n.elema), nrow = length(index.cure.v))

      for (i in c(index.cure.v)) {
        for (j in c(index.cure.v)) {
          a.sub1[i,j] <- sum((design.matrix[,i]*design.matrix[,j]*theta*(1-theta))[survt[, 2] == 1])
          a.sub2[i,j] <- sum((design.matrix[,i]*design.matrix[,j]*kap)[survt[, 2] == 0])
        }
      }
      info.a = a.sub1 + a.sub2

      design.xt <- cbind(design.matrix, log(survt[,1]))
      n.elemb <- length(index.cure.v)*(length(index.cure.v)+1)
      b.sub <- matrix(rep(0,n.elemb), nrow = length(index.surv.v))

      for (i in c(index.cure.v)) {
        for (j in c(index.cure.v,length(index.surv.v)+1)) {
          #b.sub[i,j] <- -sum((design.matrix[,i]*design.xt[,j]*theta*(1-theta)*pi)[survt[, 2] == 0]) #for est
          #b.sub[i,j] <- -sum((design.matrix[,i]*design.xt[,j]*eps*(1-eta)*eta*(1-theta))[survt[, 2] == 0])  #for LRT
          b.sub[i,j] <- -sum((design.matrix[,i]*design.xt[,j]*eps*(1-delta)*delta)[survt[, 2] == 0]) #alternative expression for est

        }
      }
      info.b = b.sub  #Upper right block of fisher.info;


      n.elemd <- (length(index.surv.v)+1)^2
      d.sub1 <- matrix(rep(0,n.elemd), nrow = (length(index.surv.v)+1))
      d.sub2 <- matrix(rep(0,n.elemd), nrow = (length(index.surv.v)+1))

      for (i in c(index.cure.v,length(index.surv.v)+1)) {
        for (j in c(index.cure.v,length(index.surv.v)+1)) {
          d.sub1[i,j] <- sum((design.xt[,i]*design.xt[,j]*eps)[survt[, 2] == 1])
          #d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta^2))[survt[, 2] == 0]) #for est, PLCI
          #d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta-eps^2*(delta*(1-delta))))[survt[, 2] == 0]) #for LRT, same as below
          d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0]) #for est, PLCI

        }
      }
      info.d = d.sub1 + d.sub2 +
        matrix(c(rep(0, (n.elemd-1)),sum(survt[, 2] == 1)/(p[index.gamma]^2)),nrow = (length(index.surv.v)+1))


      info.d.inv = mat.inv(info.d)

      #    fisher.info = rbind(cbind(info.a,info.b),cbind(t(info.b),info.d))
      #hessian.mat = -fisher.info

      # #info.set0 is (A-BD^-1B^T), dif than used in modified score;
      info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)

      #determinant of hessian matrix;
      det.info = det(info.set0)*det(info.d)
      #   det.info = matrix.det(fisher.info)


      loglik = loglikelihood - 0.5*log(det.info)
    }

    else {
      loglik = loglikelihood
    }

    #loglik = loglikelihood
    return(loglik)

  }

  ######END of loglik.mixture####################################


  maximizer <- nlm(
    f = loglik.mixture,
    p = init,
    survt=survt, design.matrix=design.matrix,
    index.cure.v=index.cure.v,
    index.surv.v=index.surv.v,
    pl = pl,
    iterlim = iterlim, hessian=F);


  loglik.part <- -maximizer$minimum  #in loglik function loglik was calculated as minus of actual loglik value
  parm.est <- maximizer$estimate

  coef.vec <- c(k, kval, loglik.part, parm.est)


  out <- list(
    coefficients = list(coef.vec)
    #cov = var.mat
  );
  #class(out) <- c('mixcure', 'list');

  return(out);
}


