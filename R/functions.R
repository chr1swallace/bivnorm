##' @importFrom mvtnorm dmvnorm rmvnorm
NULL

##' Closed form Cholesky of 2x2 variance matrices
##'
##' Let Sigma=((v1 v12), (v12 v2))
##' Calculate A st A=((a 0),(b c)) and AA'=S
##' @title 2D Cholesky decomposition
##' @param v1 element of Sigma
##' @param v2 element of Sigma
##' @param v12 element of Sigma
##' @return elements a, b, c of Cholesky decomposition
##' @author Chris Wallace
chol.abc <- function(v1,v2,v12) {
  A <- sqrt(v1)
  B <- v12/A
  C <- sqrt(v2-B^2)
  cbind(A,C,B)
}

checkv <- function(n,mean,v1,v2,v12) {
  if(length(v1)==1)
    v1 <- rep(v1,n)
  if(length(v2)==1)
    v2 <- rep(v2,n)
  if(length(v12)==1)
    v12 <- rep(v12,n)
  stopifnot(length(v1)==n &&
            length(v2)==n && 
            length(v12)==n)
  if(is.vector(mean)) {
    stopifnot(length(mean)==2)
    mean <- matrix(mean,n,2,byrow=TRUE)
  }
  stopifnot(is.matrix(mean) & nrow(mean)==n)
  list(mean=mean,v1=v1,v2=v2,v12=v12)
}

ltri.inv <- function(A)
  cbind(A[,2],A[,1],-A[,3])/(A[,1]*A[,2])

ltri.mult <- function(A,Z) {
  if(is.vector(Z))
    Z <- matrix(Z,nrow(A),2,byrow=TRUE)
  cbind(A[,1]*Z[,1],
        A[,3]*Z[,1] + A[,2] * Z[,2])
}
  
##' rmvnorm in bivariate case
##'
##' Let Sigma=((v1 v12)', (v12 v2)'), mu=(m1 m2)'
##'
##' If Z ~ N(0,I2) is the bivariate standard normal, then
##' Y = A %*% Z + mu is a bivariate normal with mean mu and variance
##' Sigma = AA'
##' @title bivariate random variables with non-constant mean and var
##' @param n number of replicates
##' @param mean mean of the ditribution
##' @param v1 element of Sigma
##' @param v2 element of Sigma
##' @param v12 element of Sigma
##' @param ... passed to rmvnorm
##' @return matrix of random samples
##' @export
##' @author Chris Wallace
rbivnorm <- function(n, mean=rep(0,2), v1=1, v2=1, v12=0,...) {
  ##check
  V <- checkv(n,mean,v1,v2,v12)
  ## calc A
  A <- chol.abc(V$v1,V$v2,V$v12)
  Z <- rmvnorm(n,sigma=diag(2),...)
  V$mean + ltri.mult(A, Z)
}

## helpful for debugging
## create symmetric matrix
sym <- function(a1,a2,a12)
  matrix(c(a1,a12,a12,a2),2,2)
## create lower triangular matrix
ltri <- function(a1,a2,a12)
  matrix(c(a1,a12,0,a2),2,2)

##' @param x matrix of quantiles
##' @inheritParams rbivnorm
##' @param ... passed to dmvnorm
##' @return
##' @export
##' @rdname rbivnorm
##' @author Chris Wallace
dbivnorm <- function(x,mean=rep(0,2), v1=1, v2=1, v12=0,...) {
  ##check
  V <- checkv(nrow(x),mean,v1,v2,v12)
  ## calc A
  A <- chol.abc(V$v1,V$v2,V$v12) # lower triangular matrix
  A1 <- ltri.inv(A) # also lower triangular
  s=sym(V$v1[1],V$v2[1],V$v12[1])
  s.A=t(chol(s))
  s.A1=solve(s.A)
  m.A=ltri(A[1,1],A[1,2],A[1,3])
  m.A1=ltri(A1[1,1],A1[1,2],A1[1,3])
  ## create Z
  Z <- ltri.mult(A1, x - V$mean)
  ## likelihood
  dmvnorm(Z,mean=rep(0,2),sigma=diag(2),...)/sqrt(V$v1*V$v2 - V$v12^2)
}
