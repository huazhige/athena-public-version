// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../coordinates/coordinates.hpp"
#include "../../../thermodynamics/thermodynamics.hpp"
#include "../../../utils/utils.hpp"
#include "../../../hydro/hydro.hpp"
#include "bvals_hydro.hpp"

void HydroBoundaryVariable::OutflowInnerX1(
    Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh) {
  // extend velocities
  for (int n = IVX; n <= IVZ; ++n) {
    for (int k=kl; k<=ku; ++k)
      for (int j=jl; j<=ju; ++j)
#pragma omp simd
        for (int i=1; i<=ngh; ++i)
          (*var_cc)(n,k,j,il-i) = (*var_cc)(n,k,j,il);
  }

  // conforms to the case without gravity
  if (pmy_block_->phydro->hsrc.GetG1() == 0.) {
    for (int n = 0; n < NMASS; ++n) {
      for (int k=kl; k<=ku; ++k)
        for (int j=jl; j<=ju; ++j)
#pragma omp simd
          for (int i=1; i<=ngh; ++i)
            (*var_cc)(n,k,j,il-i) = (*var_cc)(n,k,j,il);
    }
    for (int k=kl; k<=ku; ++k)
      for (int j=jl; j<=ju; ++j)
#pragma omp simd
        for (int i=1; i<=ngh; ++i)
          (*var_cc)(IPR,k,j,il-i) = (*var_cc)(IPR,k,j,il);
  }
}

void HydroBoundaryVariable::OutflowOuterX1(
    Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // extend velocities
  for (int n = IVX; n <= IVZ; ++n) {
    for (int k=kl; k<=ku; ++k)
      for (int j=jl; j<=ju; ++j)
#pragma omp simd
        for (int i=1; i<=ngh; ++i)
          (*var_cc)(n,k,j,iu+i) = (*var_cc)(n,k,j,iu);
  }

  Real **w1;
  NewCArray(w1, 1+ngh, NHYDRO);
  Coordinates *pcoord = pmy_block_->pcoord;
  Thermodynamics *pthermo = pmy_block_->pthermo;
  Real grav = - pmy_block_->phydro->hsrc.GetG1();

  // extrapolation based on temperature gradient
  for (int k=kl; k<=ku; ++k)
    for (int j=jl; j<=ju; ++j) {
      Real P1 = (*var_cc)(IPR,k,j,iu);
      Real T1 = pthermo->Temp(var_cc->at(k,j,iu));
      Real T2 = pthermo->Temp(var_cc->at(k,j,iu-1));
      Real dz = pcoord->dx1f(iu);
      Real dTdz = (T1 - T2)/dz;
      for (int n = 0; n < NMASS; ++n)
        w1[0][n] = (*var_cc)(n,k,j,iu);
      pthermo->ConstructAdiabat(w1, T1, P1, grav, dz, 1+ngh, Adiabat::reversible, dTdz);
#pragma omp simd
      for (int i = 1; i <= ngh; ++i) {
        for (int n = 0; n < NMASS; ++i)
          (*var_cc)(n,k,j,iu+i) = w1[i][n];
        (*var_cc)(IPR,k,j,iu+i) = w1[i][IPR];
      }
    }

  FreeCArray(w1);
}
