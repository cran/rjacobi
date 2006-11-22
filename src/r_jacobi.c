/* jacobi.h
 * Copyright (C) 2006 Paulo José Saiz Jabardo
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  Paulo J. Saiz Jabardo */



/** 
   This file presents the R interface for the jacobi library.
   Originally this library was supposed to be used in conjunction
   with the GNU Scientific Library. To remove this dependency when
   using the R interface I took some (mostly) declarations from GSL
   and removed the dependence.
*/


#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "jacobi.h"

SEXP 
R_jacobi(SEXP x, SEXP n, SEXP a, SEXP b)
{
  int nn, np;
  SEXP ans;
  double *aux_mem;

  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(a = coerceVector(a, REALSXP));
  PROTECT(b = coerceVector(b, REALSXP));
  PROTECT(n = coerceVector(n, INTSXP));
  np = length(x);
  nn = INTEGER(n)[0];
  PROTECT(ans = allocVector(REALSXP, np));

  aux_mem = (double *) R_alloc(2*np, sizeof(double));

 jac_jacobi_array(np, REAL(x), nn, REAL(ans),
		  REAL(a)[0], REAL(b)[0], aux_mem);
 UNPROTECT(5);
 return ans;
}


SEXP 
R_djacobi(SEXP x, SEXP n, SEXP a, SEXP b)
{
  int nn, np;
  SEXP ans;
  double *aux_mem;

  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(a = coerceVector(a, REALSXP));
  PROTECT(b = coerceVector(b, REALSXP));
  PROTECT(n = coerceVector(n, INTSXP));
  np = length(x);
  nn = INTEGER(n)[0];
  PROTECT(ans = allocVector(REALSXP, np));

  aux_mem = (double *) R_alloc(2*np, sizeof(double));

 jac_djacobi_array(np, REAL(x), nn, REAL(ans),
		  REAL(a)[0], REAL(b)[0], aux_mem);
 UNPROTECT(5);
 return ans;
}


SEXP 
R_jacobi_zeros(SEXP Q, SEXP alfa, SEXP beta)
{

  int q;
  SEXP ans;

  PROTECT(Q = coerceVector(Q, INTSXP));
  PROTECT(alfa = coerceVector(alfa, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  q = INTEGER(Q)[0];
  PROTECT(ans = allocVector(REALSXP, q));

  jac_jacobi_zeros(REAL(ans), q, REAL(alfa)[0], REAL(beta)[0]);

  UNPROTECT(4);

  return ans;
}




/** Calculates the Gauss-Jacobi quadrature points*/

SEXP R_zeros_gj(SEXP Q, SEXP alpha, SEXP beta)
{
  int q;
  SEXP ans;

  PROTECT(Q = coerceVector(Q, INTSXP));
  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  q = INTEGER(Q)[0];
  PROTECT(ans = allocVector(REALSXP, q));

  jac_zeros_gj(REAL(ans), q, REAL(alpha)[0], REAL(beta)[0]);

  UNPROTECT(4);

  return ans;
}


SEXP R_zeros_glj(SEXP Q, SEXP alpha, SEXP beta)
{
  int q;
  SEXP ans;

  PROTECT(Q = coerceVector(Q, INTSXP));
  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  q = INTEGER(Q)[0];
  PROTECT(ans = allocVector(REALSXP, q));

  jac_zeros_glj(REAL(ans), q, REAL(alpha)[0], REAL(beta)[0]);

  UNPROTECT(4);

  return ans;
}
  
  


SEXP R_zeros_grjm(SEXP Q, SEXP alpha, SEXP beta)
{
  int q;
  SEXP ans;

  PROTECT(Q = coerceVector(Q, INTSXP));
  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  q = INTEGER(Q)[0];
  PROTECT(ans = allocVector(REALSXP, q));

  jac_zeros_grjm(REAL(ans), q, REAL(alpha)[0], REAL(beta)[0]);

  UNPROTECT(4);

  return ans;
}


SEXP R_zeros_grjp(SEXP Q, SEXP alpha, SEXP beta)
{
  int q;
  SEXP ans;

  PROTECT(Q = coerceVector(Q, INTSXP));
  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  q = INTEGER(Q)[0];
  PROTECT(ans = allocVector(REALSXP, q));

  jac_zeros_grjp(REAL(ans), q, REAL(alpha)[0], REAL(beta)[0]);

  UNPROTECT(4);

  return ans;
}



SEXP R_weights_gj(SEXP x, SEXP alpha, SEXP beta)
{
  int q;
  SEXP ans;
  double *aux_mem;
  
  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  PROTECT(x = coerceVector(x, REALSXP));

  q = length(x);
  aux_mem = (double *) R_alloc(2*q, sizeof(double));
  

  PROTECT(ans = allocVector(REALSXP, q));

  jac_weights_gj(REAL(x), REAL(ans), q, REAL(alpha)[0], REAL(beta)[0], aux_mem);

  UNPROTECT(4);

  return ans;
}


SEXP R_weights_glj(SEXP x, SEXP alpha, SEXP beta)
{
  int q;
  SEXP ans;
  double *aux_mem;
  
  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  PROTECT(x = coerceVector(x, REALSXP));

  q = length(x);
  aux_mem = (double *) R_alloc(2*q, sizeof(double));
  

  PROTECT(ans = allocVector(REALSXP, q));

  jac_weights_glj(REAL(x), REAL(ans), q, REAL(alpha)[0], REAL(beta)[0], aux_mem);

  UNPROTECT(4);

  return ans;
}



SEXP R_weights_grjm(SEXP x, SEXP alpha, SEXP beta)
{
  int q;
  SEXP ans;
  double *aux_mem;
  
  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  PROTECT(x = coerceVector(x, REALSXP));

  q = length(x);
  aux_mem = (double *) R_alloc(2*q, sizeof(double));
  

  PROTECT(ans = allocVector(REALSXP, q));

  jac_weights_grjm(REAL(x), REAL(ans), q, REAL(alpha)[0], REAL(beta)[0], aux_mem);

  UNPROTECT(4);

  return ans;
}



SEXP R_weights_grjp(SEXP x, SEXP alpha, SEXP beta)
{
  int q;
  SEXP ans;
  double *aux_mem;
  
  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  PROTECT(x = coerceVector(x, REALSXP));

  q = length(x);
  aux_mem = (double *) R_alloc(2*q, sizeof(double));
  

  PROTECT(ans = allocVector(REALSXP, q));

  jac_weights_grjp(REAL(x), REAL(ans), q, REAL(alpha)[0], REAL(beta)[0], aux_mem);

  UNPROTECT(4);

  return ans;
}






SEXP R_diffmat_gj(SEXP x, SEXP alpha, SEXP beta)
{
  int q;
  SEXP ans;
  double *aux_mem;
  
  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  PROTECT(x = coerceVector(x, REALSXP));

  q = length(x);
  aux_mem = (double *) R_alloc(3*q, sizeof(double));
  

  PROTECT(ans = allocMatrix(REALSXP, q, q));
  
  jac_diffmat_gj(REAL(x), REAL(ans), q, REAL(alpha)[0], REAL(beta)[0], aux_mem);
  

  UNPROTECT(4);

  return ans;
}



SEXP R_diffmat_glj(SEXP x, SEXP alpha, SEXP beta)
{
  int q;
  SEXP ans;
  double *aux_mem;
  
  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  PROTECT(x = coerceVector(x, REALSXP));

  q = length(x);
  aux_mem = (double *) R_alloc(3*q, sizeof(double));
  

  PROTECT(ans = allocMatrix(REALSXP, q, q));
  
  jac_diffmat_glj(REAL(x), REAL(ans), q, REAL(alpha)[0], REAL(beta)[0], aux_mem);
  

  UNPROTECT(4);

  return ans;
}


SEXP R_diffmat_grjm(SEXP x, SEXP alpha, SEXP beta)
{
  int q;
  SEXP ans;
  double *aux_mem;
  
  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  PROTECT(x = coerceVector(x, REALSXP));

  q = length(x);
  aux_mem = (double *) R_alloc(3*q, sizeof(double));
  

  PROTECT(ans = allocMatrix(REALSXP, q, q));
  
  jac_diffmat_grjm(REAL(x), REAL(ans), q, REAL(alpha)[0], REAL(beta)[0], aux_mem);
  

  UNPROTECT(4);

  return ans;
}


SEXP R_diffmat_grjp(SEXP x, SEXP alpha, SEXP beta)
{
  int q;
  SEXP ans;
  double *aux_mem;
  
  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  PROTECT(x = coerceVector(x, REALSXP));

  q = length(x);
  aux_mem = (double *) R_alloc(3*q, sizeof(double));
  

  PROTECT(ans = allocMatrix(REALSXP, q, q));
  
  jac_diffmat_grjp(REAL(x), REAL(ans), q, REAL(alpha)[0], REAL(beta)[0], aux_mem);
  

  UNPROTECT(4);

  return ans;
}





SEXP R_interpmat_glj(SEXP zp, SEXP x, SEXP alpha, SEXP beta)
{

  int q;
  int np;

  SEXP ans;

  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(zp = coerceVector(zp, REALSXP));
  
  q = length(x);
  np = length(zp);
  
  PROTECT(ans = allocMatrix(REALSXP, q, np));

  jac_interpmat_glj(REAL(ans), REAL(zp), np, REAL(x), q, REAL(alpha)[0], REAL(beta)[0]);

  
  
  UNPROTECT(5);
  return ans;
}


SEXP R_interpmat_gj(SEXP zp, SEXP x, SEXP alpha, SEXP beta)
{

  int q;
  int np;

  SEXP ans;

  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(zp = coerceVector(zp, REALSXP));
  
  q = length(x);
  np = length(zp);

  PROTECT(ans = allocMatrix(REALSXP, q, np));

  jac_interpmat_gj(REAL(ans), REAL(zp), np, REAL(x), q, REAL(alpha)[0], REAL(beta)[0]);

  
  
  UNPROTECT(5);
  return ans;
}



SEXP R_interpmat_grjm(SEXP zp, SEXP x, SEXP alpha, SEXP beta)
{

  int q;
  int np;

  SEXP ans;

  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(zp = coerceVector(zp, REALSXP));
  
  q = length(x);
  np = length(zp);
  
  PROTECT(ans = allocMatrix(REALSXP, q, np));

  jac_interpmat_grjm(REAL(ans), REAL(zp), np, REAL(x), q, REAL(alpha)[0], REAL(beta)[0]);

  
  
  UNPROTECT(5);
  return ans;
}




SEXP R_interpmat_grjp(SEXP zp, SEXP x, SEXP alpha, SEXP beta)
{

  int q;
  int np;

  SEXP ans;

  PROTECT(alpha = coerceVector(alpha, REALSXP));
  PROTECT(beta = coerceVector(beta, REALSXP));
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(zp = coerceVector(zp, REALSXP));
  
  q = length(x);
  np = length(zp);
  
  PROTECT(ans = allocMatrix(REALSXP, q, np));

  jac_interpmat_grjp(REAL(ans), REAL(zp), np, REAL(x), q, REAL(alpha)[0], REAL(beta)[0]);

  
  
  UNPROTECT(5);
  return ans;
}




SEXP R_lagrange_glj(SEXP i, SEXP z, SEXP x, SEXP a, SEXP b)
{
  


  int np;  /* Number of points*/
  int Q;   /* Number of quadrature points*/
  int indx, j;
  SEXP ans;

  PROTECT(z = coerceVector(z, REALSXP));
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(a = coerceVector(a, REALSXP));
  PROTECT(b = coerceVector(b, REALSXP));
  PROTECT(i = coerceVector(i, INTSXP));

  Q = length(x);
  np = length(z);
  indx = INTEGER(i)[0]-1;

  if (indx < 0 && indx >= Q)
    {
      error("i should be 1 < i <= Q");
    }
  
  
  PROTECT(ans = allocVector(REALSXP, np));

  for (j = 0; j < np; ++j)
    REAL(ans)[j] = jac_lagrange_glj(indx, REAL(z)[j], Q, REAL(x), REAL(a)[0], REAL(b)[0]);
  
  
  
  UNPROTECT(6);
  return ans;
}


SEXP R_lagrange_gj(SEXP i, SEXP z, SEXP x, SEXP a, SEXP b)
{
  


  int np;  /* Number of points*/
  int Q;   /* Number of quadrature points*/
  int indx, j;
  SEXP ans;

  PROTECT(z = coerceVector(z, REALSXP));
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(a = coerceVector(a, REALSXP));
  PROTECT(b = coerceVector(b, REALSXP));
  PROTECT(i = coerceVector(i, INTSXP));

  Q = length(x);
  np = length(z);
  indx = INTEGER(i)[0]-1;

  if (indx < 0 && indx >= Q)
    {
      error("i should be 1 < i <= Q");
    }
  
  
  PROTECT(ans = allocVector(REALSXP, np));

  for (j = 0; j < np; ++j)
    REAL(ans)[j] = jac_lagrange_gj(indx, REAL(z)[j], Q, REAL(x), REAL(a)[0], REAL(b)[0]);
  
  
  
  UNPROTECT(6);
  return ans;
}


SEXP R_lagrange_grjm(SEXP i, SEXP z, SEXP x, SEXP a, SEXP b)
{
  


  int np;  /* Number of points1*/
  int Q;   /* Number of quadrature points*/
  int indx, j;
  SEXP ans;

  PROTECT(z = coerceVector(z, REALSXP));
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(a = coerceVector(a, REALSXP));
  PROTECT(b = coerceVector(b, REALSXP));
  PROTECT(i = coerceVector(i, INTSXP));

  Q = length(x);
  np = length(z);
  indx = INTEGER(i)[0]-1;

  if (indx < 0 && indx >= Q)
    {
      error("i should be 1 < i <= Q");
    }
  
  
  PROTECT(ans = allocVector(REALSXP, np));

  for (j = 0; j < np; ++j)
    REAL(ans)[j] = jac_lagrange_grjm(indx, REAL(z)[j], Q, REAL(x), REAL(a)[0], REAL(b)[0]);
  
  
  
  UNPROTECT(6);
  return ans;
}



SEXP R_lagrange_grjp(SEXP i, SEXP z, SEXP x, SEXP a, SEXP b)
{
  


  int np;  /* Number of points*/
  int Q;   /* Number of quadrature points*/
  int indx, j;
  SEXP ans;

  PROTECT(z = coerceVector(z, REALSXP));
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(a = coerceVector(a, REALSXP));
  PROTECT(b = coerceVector(b, REALSXP));
  PROTECT(i = coerceVector(i, INTSXP));

  Q = length(x);
  np = length(z);
  indx = INTEGER(i)[0]-1;

  if (indx < 0 && indx >= Q)
    {
      error("i should be 1 < i <= Q");
    }
  
  
  PROTECT(ans = allocVector(REALSXP, np));

  for (j = 0; j < np; ++j)
    REAL(ans)[j] = jac_lagrange_grjp(indx, REAL(z)[j], Q, REAL(x), REAL(a)[0], REAL(b)[0]);
  
  
  
  UNPROTECT(6);
  return ans;
}


R_CallMethodDef callMethods[]={
  {"R_jacobi", (void *) &R_jacobi, 4},
  {"R_djacobi", (void *)&R_jacobi, 4},

  {"R_jacobi_zeros", (void *)&R_jacobi_zeros, 3},

  {"R_zeros_gj",   (void *)&R_zeros_gj, 3},
  {"R_zeros_glj",  (void *)&R_zeros_glj, 3},
  {"R_zeros_grjm", (void *)&R_zeros_grjm, 3},
  {"R_zeros_grjp", (void *)&R_zeros_grjp, 3},

  {"R_weights_gj",   (void *)&R_weights_gj, 3},
  {"R_weights_glj",  (void *)&R_weights_glj, 3},
  {"R_weights_grjm", (void *)&R_weights_grjm, 3},
  {"R_weights_grjp", (void *)&R_weights_grjp, 3},

  {"R_diffmat_gj",   (void *)&R_diffmat_gj, 3},
  {"R_diffmat_glj",  (void *)&R_diffmat_glj, 3},
  {"R_diffmat_grjm", (void *)&R_diffmat_grjm, 3},
  {"R_diffmat_grjp", (void *)&R_diffmat_grjp, 3},
  
  {"R_interpmat_gj",   (void *)&R_interpmat_gj, 4},
  {"R_interpmat_glj",  (void *)&R_interpmat_glj, 4},
  {"R_interpmat_grjm", (void *)&R_interpmat_grjm, 4},
  {"R_interpmat_grjp", (void *)&R_interpmat_grjp, 4},

  {"R_lagrange_gj",   (void *)&R_lagrange_gj, 5},
  {"R_lagrange_glj",  (void *)&R_lagrange_glj, 5},
  {"R_lagrange_grjm", (void *)&R_lagrange_grjm, 5},
  {"R_lagrange_grjp", (void *)&R_lagrange_grjp, 5},

  {NULL, NULL, 0}
};


void
R_init_rjacobi(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}


