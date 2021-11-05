

#' Score test for TCR repertorire and phenotypes
#'
#' @param Y a vector, \code{Y} should be continous or binary variable.
#' @param X a covariate matrix, each row represents a subject.
#' @param fR a variant matrix, each row represents a subject, each column represents an extracted feature of TCR repertorire.
#' @param W a feature character matrix, each row represents an extracted feature, each column represents a variant. The row number of \code{W} should be consistent with the column number of \code{fR}.
#' @param S a homology matrix, which represents the correlation between subjects. It can be cauculated by using function \code{seqhom}.
#'
#' @return
#' \describe{
#' \item{fix.effect}{p-value for testing the fixed variant effect.}
#' \item{random.effect}{p-value for testing the random variant effect.}
#' \item{Overall.pval}{overall p-value for testing the assocation between TCR repertorire and phenotype. It combines \code{fix.effect} and \code{rand.effect} by using Fisher's procedure.}
#' }
#'
#' @export
#'
#' @examples
#' data("Example.data")

#' ## extract features
#' fR <- AAfreq(TCRdat[,1],TCRdat[,2],as.numeric(TCRdat[,3]))
#'
#' ## homology matrix
#' S <- seqhom(TCRdat[,1],TCRdat[,2],as.numeric(TCRdat[,3]),'BLOSUM62')
#'
#' ## include both fixed and random effect
#' TCRL.cont(Y,X,fR,W,S)


TCRL_bin <-function(Y,X,fR,W,S){

    n<-nrow(Y)
    fRW <- fR%*%W
    r <- ncol(fRW)
    d = ncol(X)

    ## score test for fixed part
    fit0 <-glm(Y~0+X,family=binomial(link=logit))
    mu0<-fit0$fitted.values
    Omega0<-diag(mu0*(1-mu0),n,n)
    tXOmega0 <- t(X)%*%Omega0
    Gamma1<-Omega0-t(tXOmega0)%*%solve(tXOmega0%*%X)%*%tXOmega0
    U1<-matrix(t(Y-mu0)%*%fRW,r,1)
    Sigma1<-t(fRW)%*%Gamma1%*%fRW
    T1 <-t(U1)%*%solve(Sigma1)%*%U1
    pval_fix <-pchisq(T1, df=r, lower.tail=FALSE)

    ## score test for random part
    Z <- cbind(X,fRW)
    d <- ncol(Z)

    fit<-glm(Y~0+Z,family=binomial(link=logit))
    mu<-fit$fitted.values
    Omega<-diag(mu*(1-mu),n,n)
    Omega_half = diag(sqrt(mu*(1-mu)),n,n)
    ywork = Z%*%coef(fit) + (Y-mu)/(mu*(1-mu)) # working response
    newY =  Omega_half%*%ywork
    newZ <- Omega_half%*%Z
    P0 <- diag(n)- newZ%*%solve(t(newZ)%*%newZ)%*%t(newZ)
    sig2 = t(newY)%*%P0%*%newY/(n-d)
    newmu = newZ%*%coef(fit)
    SP0 = S%*%P0
    Tmu = tr(SP0)*sig2
    newT <-t(newY-newmu)%*%S%*%(newY-newmu)
    Tvar = 2*sig2^2*(tr(SP0%*%SP0) - (tr(SP0))^2/(n-d))
    TYnorm <- (newT - Tmu)/sqrt(Tvar)
    pval_rand <- pnorm(abs(TYnorm),0,1,lower.tail=FALSE)*2


    ## fisher combination
    pval.overall  <- fisher(c(pval_fix,pval_rand))

    return(c(fixed.effect=pval_fix,rand.effect=pval_rand,Overall.pval=pval.overall))

}

#' @rdname TCRL_bin
TCRL_cont <-function(Y,X,fR,W,S){

    n<-nrow(Y)
    fRW <- fR%*%W
    r <- ncol(fRW)

    ## score test for fixed part
    fit1<-lm(Y~0+X)
    mu1<-fit1$fitted.value
    sig1<-summary(fit1)$sigma**2
    Omega1<-diag(sig1,n,n)
    tXOmega1 <- t(X)%*%Omega1
    Gamma1<-Omega1-t(tXOmega1)%*%solve(tXOmega1%*%X)%*%tXOmega1
    U1<-matrix(t(Y-mu1)%*%fRW,r,1)
    Sigma1<-t(fRW)%*%Gamma1%*%fRW
    T1 <-t(U1)%*%solve(Sigma1)%*%U1
    pval.fix<-pchisq(T1, df=r, lower.tail=FALSE)

    ## score test for random part
    Z <- cbind(X,fRW)
    d <- ncol(Z)

    fit <- lm(Y~0+Z)
    resid <- Y - fit$fitted.value
    sig2 <- summary(fit)$sigma**2
    T <-t(resid)%*%S%*%(resid)
    P0 = diag(n) - Z%*%solve(t(Z)%*%Z)%*%t(Z)

    SP0 = S%*%P0
    Tmu = tr(SP0)*sig2
    Tvar = 2*sig2^2*(tr(SP0%*%SP0) - (tr(SP0))^2/(n-d))
    TYnorm <- (T - Tmu)/sqrt(Tvar)
    pval.rand <- pnorm(abs(TYnorm),0,1,lower.tail=FALSE)*2

    ## fisher combination
    pval.overall <- fisher(c(pval.fix,pval.rand))

    return(c(fixed.effect=pval1,rand.effect=pval.rand,Overall.pval = pval.overall))
}


