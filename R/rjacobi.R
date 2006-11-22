


jacobi <- function(x,n,a=0.0, b=0.0){
  .Call("R_jacobi", x, n, a, b, PACKAGE='rjacobi')
}

djacobi <- function(x,n,a=0.0, b=0.0){
  .Call("R_djacobi", x, n, a, b, PACKAGE='rjacobi')
}


jacobiZeros <- function(Q, a=0, b=0)
  .Call("R_jacobi_zeros", Q, a, b, PACKAGE='rjacobi')


zerosGJ <- function(Q, a=0, b=0)
  .Call("R_zeros_gj", Q, a, b, PACKAGE='rjacobi')

zerosGLJ <- function(Q, a=0, b=0)
  .Call("R_zeros_glj", Q, a, b, PACKAGE='rjacobi')

zerosGRJM <- function(Q, a=0, b=0)
  .Call("R_zeros_grjm", Q, a, b, PACKAGE='rjacobi')

zerosGRJP <- function(Q, a=0, b=0)
  .Call("R_zeros_grjp", Q, a, b, PACKAGE='rjacobi')


weightsGJ <- function(x, a=0,b=0)
  .Call("R_weights_gj", x, a, b, PACKAGE='rjacobi')

weightsGLJ <- function(x,a=0,b=0)
  .Call("R_weights_glj", x, a, b, PACKAGE='rjacobi')

weightsGRJM <- function(x,a=0,b=0)
  .Call("R_weights_grjm", x, a, b, PACKAGE='rjacobi')

weightsGRJP <- function(x,a=0,b=0)
  .Call("R_weights_grjp", x, a, b, PACKAGE='rjacobi')


diffmatGJ <- function(x,a=0, b=0)
  t(.Call("R_diffmat_gj", x, a, b, PACKAGE='rjacobi'))

diffmatGLJ <- function(x,a=0, b=0)
  t(.Call("R_diffmat_glj", x, a, b, PACKAGE='rjacobi'))

diffmatGRJM <- function(x,a=0, b=0)
  t(.Call("R_diffmat_grjm", x, a, b, PACKAGE='rjacobi'))

diffmatGRJP <- function(x,a=0, b=0)
  t(.Call("R_diffmat_grjp", x, a, b, PACKAGE='rjacobi'))



interpmatGJ <- function(zp, x, a=0, b=0)
  t(.Call("R_interpmat_gj", zp, x, a, b, PACKAGE='rjacobi'))

interpmatGLJ <- function(zp, x, a=0, b=0)
  t(.Call("R_interpmat_glj", zp, x, a, b, PACKAGE='rjacobi'))

interpmatGRJM <- function(zp, x, a=0, b=0)
  t(.Call("R_interpmat_grjm", zp, x, a, b, PACKAGE='rjacobi'))

interpmatGRJP <- function(zp, x, a=0, b=0)
  t(.Call("R_interpmat_grjp", zp, x, a, b, PACKAGE='rjacobi'))





# Wrappers of the functions above to ease on the use of the library
quadGJ <- function(Q, a=0, b=0){
  x <- zerosGJ(Q, a, b)
  w <- weightsGJ(x, a, b)
  D <- diffmatGJ(x, a, b)

  # Interpolation matrix, build an abstraction
  interp <- function(zp)
    interpmatGJ(zp, x, a, b)

  return(list(x=x, w=w, D=D, interp=interp))
}

quadGLJ <- function(Q, a=0, b=0){
  x <- zerosGLJ(Q, a, b)
  w <- weightsGLJ(x, a, b)
  D <- diffmatGLJ(x, a, b)

  # Interpolation matrix, build an abstraction
  interp <- function(zp)
    interpmatGLJ(zp, x, a, b)

  return(list(x=x, w=w, D=D, interp=interp))
}


quadGRJM <- function(Q, a=0, b=0){
  x <- zerosGRJM(Q, a, b)
  w <- weightsGRJM(x, a, b)
  D <- diffmatGRJM(x, a, b)

  # Interpolation matrix, build an abstraction
  interp <- function(zp)
    interpmatGRJM(zp, x, a, b)

  return(list(x=x, w=w, D=D, interp=interp))
}

         
quadGRJP <- function(Q, a=0, b=0){
  x <- zerosGRJP(Q, a, b)
  w <- weightsGRJP(x, a, b)
  D <- diffmatGRJP(x, a, b)

  # Interpolation matrix, build an abstraction
  interp <- function(zp)
    interpmatGRJP(zp, x, a, b)

  return(list(x=x, w=w, D=D, interp=interp))
}

quadrature <- function(Q, type='GJ', a=0, b=0, xinterp=NULL){
  type <- match.arg(toupper(type), c('GJ', 'GLJ', 'GRJM', 'GRJP'))

  quadfun <- get(paste('quad', type, sep=''))
  q <- quadfun(Q, a, b)
  q$a <- a
  q$b <- b
  q$type <- type
  
  lagr <- get(paste('lagrange', type, sep=''))
  lagrange <- function(x,i)
    lagr(x, i, q$x, q$a, q$b)

  q$lagrange <- lagrange
  if (!is.null(xinterp)){
    q$xinterp <- xinterp
    q$imat <- q$interp(xinterp)
  }else{
    q$xinterp <- NULL
    q$imat <- NULL
  }

  class(q) <- "quadrature"
  return(q)
}


xinterp <- function(q){
  return(q$xinterp)
}

'xinterp<-' <- function(q, value){
  q$xinterp <- value
  q$imat <- q$interp(value)
  return(q)
}

quadIntegr <- function(f,quad){
  if (is.function(f)) f <- f(quad$x)

  return(sum(f*quad$w))
}

quadDeriv <- function(f, quad){
  if (is.function(f)) f <- f(quad$x)

  return(quad$D %*% f)
}

quadInterp <- function(imat, f){
  if (class(imat) == 'quadrature'){
    if (is.function(f)) f <- f(imat$x)
    imat <- imat$imat
  }
  
  return(imat %*% f)
}


lagrangeGJ <- function(z, i, x, a=0.0, b=0.0){
  .Call("R_lagrange_gj", i, z, x, a, b, PACKAGE='rjacobi')
}


lagrangeGLJ <- function(z, i, x, a=0.0, b=0.0){
  .Call("R_lagrange_glj", i, z, x, a, b, PACKAGE='rjacobi')
}


lagrangeGRJM <- function(z, i, x, a=0.0, b=0.0){
  .Call("R_lagrange_grjm", i, z, x, a, b, PACKAGE='rjacobi')
}


lagrangeGRJP <- function(z, i, x, a=0.0, b=0.0){
  .Call("R_lagrange_grjp", i, z, x, a, b, PACKAGE='rjacobi')
}

