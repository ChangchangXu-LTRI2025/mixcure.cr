##################################
### mixcure.penal.wald.cr      ###
### Creates:  estimates and se ###
### Wald 2d test p, and CR area ##
##################################


mixcure.penal.wald.cr.area <- function(formula, data, apct = 0.05,init, true.val, pl){

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
    iterlim = 100, hessian=F);

  var.est <- maximizer0$estimate
  var.mat <- solve(maximizer0$hessmat)

  wald.2d.test <- rep(NA, length(ncol(design.matrix)))
  wald.2d.test.p <- rep(NA, length(ncol(design.matrix)))
  wald.2d.cr.area <- rep(NA, length(ncol(design.matrix)))
  corr.var <- rep(NA, length(ncol(design.matrix)))

  for (i in 1: ncol(design.matrix) ) {
    mat.index = c(i, (i+ncol(design.matrix)) )
    test.2d = (var.est-true.val)[mat.index]%*%solve(vcov.mat[mat.index, mat.index])%*%(var.est-true.val)[mat.index]
    wald.2d.test.p[i] = pchisq(test.2d, df=2,lower.tail = F)
    wald.2d.cr.area[i] = sqrt(diag(var.mat[i]))*sqrt(diag(var.mat[i+ncol(design.matrix)]))*pi
    wald.2d.test[i] <- test.2d
    corr.var[i] <- vcov.mat[mat.index, mat.index][2,1]
      }

  coef.table <- cbind(
    'coef.cure'   = var.est[index.cure.v],
    'coef.surv'   = exp(var.est[index.surv.v]),
    'se.cure'    = sqrt(diag(var.mat)[index.cure.v]),
    'se.surv'    = sqrt(diag(var.mat)[index.surv.v]),
    'corr'        = corr.var[index.surv.v],
    '2dchi'       = wald.2d.test[index.cure.v],
    'Pval'    = wald.2d.test.p[index.cure.v],
    'area_cr' = wald.2d.cr.area[index.cure.v]
   );
  rownames(coef.table) <- colnames(design.matrix);

  return(coef.table);
}
