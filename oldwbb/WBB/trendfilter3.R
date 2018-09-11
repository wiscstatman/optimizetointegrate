trendfilter3 = function(y,w,lambda,ord=1)
{
  L = length(w)
  n = length(y)
  W = sparseMatrix(i=1:n,j=1:n,x=sqrt(w[-L]/w[L]))
  
  y = as.numeric(W%*%y)
  Winv = sparseMatrix(i=1:n,j=1:n,x=sqrt(w[L]/w[-L]))
  
  
  D = getDtfSparse(n, ord)%*%Winv
    
  k = n-ord-1
  ui = matrix(0,nrow = 2*k,ncol = k)
  for(i in 1:k)
  {
    ui[2*i-1,i] = 1
    ui[2*i,i] = -1
  }
  ci = rep(-lambda,2*k)
  
  dvec <- crossprod(t(D), y)
  Dmat <- crossprod(t(D), t(D))
  diag(Dmat) <- diag(Dmat) + 1e-08
  Amat <- t(rbind(ui))
  bvec <- c(ci)
  u <- solve.QP(Dmat, dvec, Amat, bvec, meq = 0)$solution
  
  out= Winv%*%as.matrix(y - t(D) %*% u)
  return(out)
}

getDtfSparse <- function(n,ord) {
  D = bandSparse(n, m=n, c(0,1), diagonals=list(rep(-1,n),rep(1,n-1)))
  D0 = D
  for (i in Seq(1,ord)) D = D0 %*% D
  return(D[Seq(1,n-ord-1),,drop=FALSE])
}

Seq <- function(a, b, ...) {
  if (a<=b) return(seq(a,b,...))
  else return(numeric(0))
}
