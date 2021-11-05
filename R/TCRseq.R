## only random effect is considered

#' Score test for TCR repertorire and phenotypes
#'
#' @param Y a vector, \code{Y} should be continous variable.
#' @param X a covariate matrix, each row represents a subject.
#' @param S a homology matrix, which represents the correlation between subjects. It can be cauculated by using function \code{seqhom}.
#'
#' @return \item{Overall.pval}{p-value for testing the assocation between TCR repertorire and phenotype.}
#' @export
#'
#' @examples
#' data("Example.data")
#'
#' ## extract features
#' fR <- AAfreq(TCRdat[,1],TCRdat[,2],as.numeric(TCRdat[,3]))
#'
#' ## homology matrix
#' S <- seqhom(TCRdat[,1],TCRdat[,2],as.numeric(TCRdat[,3]),'BLOSUM62')
#'
#' ## include random effect only
#' TCRseq.cont(Y,X,S)


TCRseq_bin <-function(Y,X,S){

    n<-nrow(Y)
    d = ncol(X)

    ## score test for random part
    fit<-glm(Y~0+X,family=binomial(link=logit))
    mu<-fit$fitted.values
    Omega<-diag(mu*(1-mu),n,n)
    Omega_half = diag(sqrt(mu*(1-mu)),n,n)

    ywork = X%*%coef(fit) + (Y-mu)/(mu*(1-mu)) # working response
    newY =  Omega_half%*%ywork
    newX <- Omega_half%*%X
    P0 <- diag(n)- newX%*%solve(t(newX)%*%newX)%*%t(newX)
    sig2 = t(newY)%*%P0%*%newY/(n-d)
    newmu = newX%*%coef(fit)

    SP0 = S%*%P0
    Tmu = tr(SP0)*sig2
    newT <-t(newY-newmu)%*%S%*%(newY-newmu)
    Tvar = 2*sig2^2*(tr(SP0%*%SP0) - (tr(SP0))^2/(n-d))
    TYnorm <- (newT - Tmu)/sqrt(Tvar)
    pval <- pnorm(abs(TYnorm),0,1,lower.tail=FALSE)*2


    return(c(TCR.Seq = pval))
}

#' @rdname TCRseq_bin
TCRseq_cont <-function(Y,X,S){

    n <- nrow(Y)
    d = ncol(X)

    ## score test for random part, chisq approx
    fit <- lm(Y~0+X)
    resid <- Y - fit$fitted.value
    sig2 <- summary(fit)$sigma**2
    T <-t(resid)%*%S%*%(resid)
    P0 = diag(n) - X%*%solve(t(X)%*%X)%*%t(X)

    SP0 = S%*%P0
    Tmu = tr(SP0)*sig2
    Tvar = 2*sig2^2*(tr(SP0%*%SP0) - (tr(SP0))^2/(n-d))
    TYnorm <- (T - Tmu)/sqrt(Tvar)
    pval <- pnorm(abs(TYnorm),0,1,lower.tail=FALSE)*2

    return(c(TCR.seq=pval))
}




