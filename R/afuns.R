#' Compute TCR repertoire homology between subjects.

#' @param Subject.ID a vector of subject IDs.
#' @param AAseq A vector of amino acid sequences. The order of sequences should be consistent with \code{Subject.ID}.
#' @param Abundance A numeric vector of abundance. The order of abundances should be consistent with  \code{Subject.ID} and \code{AAseq}.
#' @param substitutionMatrix A character, representing the fixe substitution scores for an alignment. It can take value from following options: BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM100, PAM30, PAM40, PAM70, PAM120, and PAM250.
#'
#' @return \describe{\item{S}{a homology matrix.}}
#' @export
#'
#' @examples
#' data("Example.data")
#'
#' ## homology matrix
#' S <- seqhom(TCRdat[,1],TCRdat[,2],as.numeric(TCRdat[,3]),'BLOSUM62')
seqhom <- function(Subject.ID,AAseq,Abundance,substitutionMatrix){
    sID <- unique(Subject.ID)
    n = length(sID)
    S <- matrix(NA,n,n)
    for(i in 1:(n-1)){
        for(j in (i+1):n){
            iid <- which(Subject.ID==sID[i])
            jid <- which(Subject.ID==sID[j])
            ni = length(iid)
            nj = length(jid)
            seqmat <- matrix(NA,ni,nj)
            for(k in 1:ni){
                iseq <- AAseq[iid][k]
                sk <- Biostrings::pairwiseAlignment(iseq,iseq, substitutionMatrix = substitutionMatrix)@score
                for(m in 1:nj){
                    jseq <- AAseq[jid][m]
                    skm <- Biostrings::pairwiseAlignment(iseq, jseq, substitutionMatrix = substitutionMatrix)@score
                    sm <- Biostrings::pairwiseAlignment(jseq, jseq, substitutionMatrix = substitutionMatrix)@score
                    seqmat[k,m] <- skm/sqrt(sk*sm)
                }
            }
            abund1 <- Abundance[iid]
            abund2 <- Abundance[jid]
            s1 = sum(apply(seqmat,1,max)*abund1)
            s2 = sum(apply(seqmat,2,max)*abund2)
            S[i,j] <- S[j,i] <- (s1+s2)/sum(c(abund1,abund2))
        }
    }
    diag(S) = 1
    if(!is.positive.semi.definite(S)) S <- psd(S)
    return(S)
}

psd <- function(S){
  if (isSymmetric(S)==FALSE){
    stop("S should be symmetric.")
  }
  es = eigen(S)
  S <- as.matrix(S)
  Q <- es$values
  m <- sum(Q>=0)
  P <- es$vectors
  n <- dim(S)[1]
  my.psd <- matrix(0, n, n)
  rownames(my.psd) <- rownames(S)
  colnames(my.psd) <- colnames(S)
  for (j in 1: m){
    my.psd <- my.psd + Q[j]*P[,j]%*%t(P[,j])
  }
  return(my.psd)
}


#' Compute amino acid frequence in TCR repertoire.
#'
#' @param Subject.ID a vector of subject IDs.
#' @param AAseq A vector of amino acid sequences. The order of sequences should be consistent with \code{Subject.ID}.
#' @param Abundance A numeric vector of abundance. The order of abundances should be consistent with  \code{Subject.ID} and \code{AAseq}.
#'
#' @return \describe{\item{out}{an amino acid matrix with each row represents a subject and each column represents an animo acid.}}
#' @export
#'
#' @examples
#' data("Example.data")
#'
#' ## extract features
#' fR <- AAfreq(TCRdat[,1],TCRdat[,2],as.numeric(TCRdat[,3]))

AAfreq <- function(Subject.ID,AAseq,Abundance){ ## extract features
    sID <- unique(Subject.ID)
    n = length(sID)
    out = matrix(NA,n,20)
    for(i in 1:n){
        iid <- which(Subject.ID==sID[i])
        sAA = strsplit(paste(rep(AAseq[iid],Abundance[iid]),collapse=""),"")[[1]]
        out[i,] <- table(factor(sAA,levels=AA_STANDARD))/length(sAA)
    }
    return(out)
}


fisher <- function(pvec){
  n = length(pvec)
  q = -2*(sum(log(pvec)))
  1- pchisq(q,2*n)
}

tr <- function(x) sum(diag(x))
