//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file weno5_simple.cpp
//  \brief  WENO5 interpolation

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "interpolation.hpp"
#include "interp_weno5.hpp"
#include "interp_weno3.hpp"

#ifdef UNIFORM_GRID
  #define interp_weno5x1 interp_weno5
#endif

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::Weno5X1()
//  \brief 

void Reconstruction::Weno5X1(const int k, const int j, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  for (int n=0; n<=NVAPOR; ++n) 
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real scale = 0.;
      for (int m = -2; m <= 2; ++m) scale += (w(n,k,j,i+m) + 1.E-16)/5.;
      wl(n,i+1) = interp_weno5x1(w(n,k,j,i+2)/scale, w(n,k,j,i+1)/scale, w(n,k,j,i)/scale,
                                 w(n,k,j,i-1)/scale, w(n,k,j,i-2)/scale);
      wl(n,i+1) *= scale;
      wr(n,i) = interp_weno5x1(w(n,k,j,i-2)/scale, w(n,k,j,i-1)/scale, w(n,k,j,i)/scale,
                               w(n,k,j,i+1)/scale, w(n,k,j,i+2)/scale);
      wr(n,i) *= scale;
    }

  for (int n=NVAPOR+1; n<NMASS; ++n) 
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real scale = 0.;
      for (int m = -1; m <= 1; ++m) scale += (w(n,k,j,i+m) + 1.E-16)/3.;
      wl(n,i+1) = interp_weno3x1(w(n,k,j,i+1)/scale, w(n,k,j,i)/scale, w(n,k,j,i-1)/scale);
      wl(n,i+1) *= scale;

      wr(n,i) = interp_weno3x1(w(n,k,j,i-1)/scale, w(n,k,j,i)/scale, w(n,k,j,i+1)/scale);
      wr(n,i) *= scale;
    }

  for (int n=NMASS; n<NHYDRO; ++n) {
#pragma omp simd
    for (int i=il; i<=il+2; ++i) {
      wl(n,i+1) = interp_weno5x1(w(n,k,j,i+2),w(n,k,j,i+1),w(n,k,j,i),w(n,k,j,i-1),w(n,k,j,i-2));
      wr(n,i) = interp_weno5x1(w(n,k,j,i-2),w(n,k,j,i-1),w(n,k,j,i),w(n,k,j,i+1),w(n,k,j,i+2));
    }

#pragma omp simd
    for (int i=il+3; i<=iu-3; ++i) {
      wl(n,i+1) = interp_cp5x1(w(n,k,j,i+2),w(n,k,j,i+1),w(n,k,j,i),w(n,k,j,i-1),w(n,k,j,i-2));
      wr(n,i) = interp_cp5x1(w(n,k,j,i-2),w(n,k,j,i-1),w(n,k,j,i),w(n,k,j,i+1),w(n,k,j,i+2));
    }

#pragma omp simd
    for (int i=iu-2; i<=iu; ++i) {
      wl(n,i+1) = interp_weno5x1(w(n,k,j,i+2),w(n,k,j,i+1),w(n,k,j,i),w(n,k,j,i-1),w(n,k,j,i-2));
      wr(n,i) = interp_weno5x1(w(n,k,j,i-2),w(n,k,j,i-1),w(n,k,j,i),w(n,k,j,i+1),w(n,k,j,i+2));
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::Weno5X2()
//  \brief 

void Reconstruction::Weno5X2(const int k, const int j, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  for (int n=0; n<=NVAPOR; ++n)
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real scale = 0.;
      for (int m = -2; m <= 2; ++m) scale += (w(n,k,j+m,i) + 1.E-16)/5.;
      wl(n,i) = interp_weno5(w(n,k,j+2,i)/scale, w(n,k,j+1,i)/scale, w(n,k,j,i)/scale,
                             w(n,k,j-1,i)/scale, w(n,k,j-2,i)/scale);
      wl(n,i) *= scale;
      wr(n,i) = interp_weno5(w(n,k,j-2,i)/scale, w(n,k,j-1,i)/scale, w(n,k,j,i)/scale,
                             w(n,k,j+1,i)/scale, w(n,k,j+2,i)/scale);
      wr(n,i) *= scale;
    }

  for (int n=NVAPOR+1; n<NMASS; ++n) 
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real scale = 0.;
      for (int m = -1; m <= 1; ++m) scale += (w(n,k,j+m,i) + 1.E-16)/3.;
      wl(n,i) = interp_weno3(w(n,k,j+1,i)/scale, w(n,k,j,i)/scale, w(n,k,j-1,i)/scale);
      wl(n,i) *= scale;
      wr(n,i) = interp_weno3(w(n,k,j-1,i)/scale, w(n,k,j,i)/scale, w(n,k,j+1,i)/scale);
      wr(n,i) *= scale;
    }

  for (int n=NMASS; n<NHYDRO; ++n) {
#pragma omp simd
    for (int i=il; i<=il+2; ++i) {
      wl(n,i) = interp_weno5(w(n,k,j+2,i),w(n,k,j+1,i),w(n,k,j,i),w(n,k,j-1,i),w(n,k,j-2,i));
      wr(n,i) = interp_weno5(w(n,k,j-2,i),w(n,k,j-1,i),w(n,k,j,i),w(n,k,j+1,i),w(n,k,j+2,i));
    }
#pragma omp simd
    for (int i=il+3; i<=iu-3; ++i) {
      wl(n,i) = interp_cp5(w(n,k,j+2,i),w(n,k,j+1,i),w(n,k,j,i),w(n,k,j-1,i),w(n,k,j-2,i));
      wr(n,i) = interp_cp5(w(n,k,j-2,i),w(n,k,j-1,i),w(n,k,j,i),w(n,k,j+1,i),w(n,k,j+2,i));
    }
#pragma omp simd
    for (int i=iu-2; i<=iu; ++i) {
      wl(n,i) = interp_weno5(w(n,k,j+2,i),w(n,k,j+1,i),w(n,k,j,i),w(n,k,j-1,i),w(n,k,j-2,i));
      wr(n,i) = interp_weno5(w(n,k,j-2,i),w(n,k,j-1,i),w(n,k,j,i),w(n,k,j+1,i),w(n,k,j+2,i));
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::Weno5X3()
//  \brief 

void Reconstruction::Weno5X3(const int k, const int j, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  for (int n=0; n<=NVAPOR; ++n)
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real scale = 0.;
      for (int m = -2; m <= 2; ++m) scale += (w(n,k+m,j,i) + 1.E-16)/5.;
      wl(n,i) = interp_weno5(w(n,k+2,j,i)/scale, w(n,k+1,j,i)/scale, w(n,k,j,i)/scale,
                             w(n,k-1,j,i)/scale, w(n,k-2,j,i)/scale);
      wl(n,i) *= scale;
      wr(n,i) = interp_weno5(w(n,k-2,j,i)/scale, w(n,k-1,j,i)/scale, w(n,k,j,i)/scale,
                             w(n,k+1,j,i)/scale, w(n,k+2,j,i)/scale);
      wr(n,i) *= scale;
    }

  for (int n=NVAPOR+1; n<NMASS; ++n) 
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      Real scale = 0.;
      for (int m = -1; m <= 1; ++m) scale += (w(n,k+m,j,i) + 1.E-16)/3.;
      wl(n,i) = interp_weno3(w(n,k+1,j,i)/scale, w(n,k,j,i)/scale, w(n,k-1,j,i)/scale);
      wl(n,i) *= scale;
      wr(n,i) = interp_weno3(w(n,k-1,j,i)/scale, w(n,k,j,i)/scale, w(n,k+1,j,i)/scale);
      wr(n,i) *= scale;
    }

  for (int n=NMASS; n<NHYDRO; ++n) {
#pragma omp simd
    for (int i=il; i<=il+2; ++i) {
      wl(n,i) = interp_weno5(w(n,k+2,j,i),w(n,k+1,j,i),w(n,k,j,i),w(n,k-1,j,i),w(n,k-2,j,i));
      wr(n,i) = interp_weno5(w(n,k-2,j,i),w(n,k-1,j,i),w(n,k,j,i),w(n,k+1,j,i),w(n,k+2,j,i));
    }
#pragma omp simd
    for (int i=il+3; i<=iu-3; ++i) {
      wl(n,i) = interp_cp5(w(n,k+2,j,i),w(n,k+1,j,i),w(n,k,j,i),w(n,k-1,j,i),w(n,k-2,j,i));
      wr(n,i) = interp_cp5(w(n,k-2,j,i),w(n,k-1,j,i),w(n,k,j,i),w(n,k+1,j,i),w(n,k+2,j,i));
    }
#pragma omp simd
    for (int i=iu-2; i<=iu; ++i) {
      wl(n,i) = interp_weno5(w(n,k+2,j,i),w(n,k+1,j,i),w(n,k,j,i),w(n,k-1,j,i),w(n,k-2,j,i));
      wr(n,i) = interp_weno5(w(n,k-2,j,i),w(n,k-1,j,i),w(n,k,j,i),w(n,k+1,j,i),w(n,k+2,j,i));
    }
  }

  return;
}
