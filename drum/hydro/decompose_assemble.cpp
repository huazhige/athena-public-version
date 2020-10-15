// Athena++ headers
#include "hydro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../coordinates/coordinates.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../globals.hpp"

inline void IntegrateUpwards(AthenaArray<Real>& psf, AthenaArray<Real> const& w, Coordinates *pco,
  Real grav, int il, int iu, int jl, int ju, int kl, int ku)
{
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i) {
        psf(k,j,i+1) = psf(k,j,i) + grav*w(IDN,k,j,i)*pco->dx1f(i);
        if (psf(k,j,i+1) < 0.)
          psf(k,j,i+1) = psf(k,j,i)*exp(grav*w(IDN,k,j,i)*pco->dx1f(i)/w(IPR,k,j,i));
      }

}

inline void IntegrateDownwards(AthenaArray<Real>& psf, AthenaArray<Real> const& w, Coordinates *pco,
  Real grav, int il, int iu, int jl, int ju, int kl, int ku)
{
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = iu; i >= il; --i)
        psf(k,j,i) = psf(k,j,i+1) - grav*w(IDN,k,j,i)*pco->dx1f(i);
}

void Hydro::DecomposePressure(AthenaArray<Real> &w, int kl, int ku, int jl, int ju)
{
  MeshBlock *pmb = pmy_block;
  Coordinates *pco = pmb->pcoord;
  Thermodynamics *pthermo = pmb->pthermo;

  Real grav = hsrc.GetG1();  // positive in the x1-increasing direction
  Real Rd = pthermo->GetRd();

  int is = pmb->is, ie = pmb->ie;
  
  if (grav == 0.) return;

  // find top and bot neighbor
  NeighborBlock ntop, nbot;
  bool has_top_neighbor = false;
  bool has_bot_neighbor = false;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock& nb = pmb->pbval->neighbor[n];
    if ((nb.ni.ox1 == -1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0)) {
      nbot = nb;
      has_bot_neighbor = true;
    } if ((nb.ni.ox1 == 1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0)) {
      ntop = nb;
      has_top_neighbor = true;
    }
  }

  // calculate local polytropic index
  pthermo->PolytropicIndex(gamma_, w);

  if (has_top_neighbor)
    RecvTopPressure(psf_, psbuf_, ntop);
  else {  // isothermal extrapolation to find the pressure at top boundary
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j) {
        Real Tv = w(IPR,k,j,ie)/(w(IDN,k,j,ie)*Rd);
        psf_(k,j,ie+1) = w(IPR,k,j,ie)*exp(grav*pco->dx1f(ie)/(2.*Rd*Tv));
      }
  }
  IntegrateDownwards(psf_, w, pco, grav, is, ie, jl, ju, kl, ku);
  
  if (has_bot_neighbor)
    SendBotPressure(psf_, psbuf_, nbot);

  // integrate ghost cells
  if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::reflect)
    IntegrateDownwards(psf_, w, pco, -grav, is - NGHOST, is-1, jl, ju, kl, ku);
  else  // block boundary
    IntegrateDownwards(psf_, w, pco,  grav, is - NGHOST, is-1, jl, ju, kl, ku);

  if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::reflect)
    IntegrateUpwards(psf_, w, pco, -grav, ie + 1, ie + NGHOST, jl, ju, kl, ku);
  else  // block boundary
    IntegrateUpwards(psf_, w, pco,  grav, ie + 1, ie + NGHOST, jl, ju, kl, ku);

  // decompose pressure
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j) {
      // 1. change density and pressure (including ghost cells)
      for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
        // interpolate hydrostatic pressure, prevent divided by zero
        if (fabs(psf_(k,j,i) - psf_(k,j,i+1)) < 1.E-6)
          psv_(k,j,i) = 0.5*(psf_(k,j,i) + psf_(k,j,i+1));
        else
          psv_(k,j,i) = (psf_(k,j,i) - psf_(k,j,i+1))/log(psf_(k,j,i)/psf_(k,j,i+1));

        // change density to pseudo entropy
        w(IDN,k,j,i) = log(w(IPR,k,j,i)) - gamma_(k,j,i)*log(w(IDN,k,j,i));
        
        // change pressure to pertubation pressure
        w(IPR,k,j,i) -= psv_(k,j,i);  
      }

      /* 2. adjust bottom boundary condition
      if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
        for (int i = is - NGHOST; i < is; ++i) {
          w(IDN,k,j,i) = w(IDN,k,j,is);
          w(IPR,k,j,i) = w(IPR,k,j,is);
        }
      }

      // 3. adjust top boundary condition
      if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::outflow) {
        for (int i = ie + 1; i <= ie + NGHOST; ++i) {
          w(IDN,k,j,i) = w(IDN,k,j,ie);
          w(IPR,k,j,i) = w(IPR,k,j,ie);
        }
      }*/
    }

  // finish send bottom pressure
  if (has_bot_neighbor && (nbot.snb.rank != Globals::my_rank))
    WaitBotPressure();
}

void Hydro::AssemblePressure(AthenaArray<Real> &w,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr, int k, int j, int il, int iu)
{
  MeshBlock *pmb = pmy_block;
  int is = pmb->is, ie = pmb->ie;
  Real grav = hsrc.GetG1();
  if (grav == 0.) return;
  
  for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
    w(IPR,k,j,i) += psv_(k,j,i);
    if (w(IPR,k,j,i) < 0.) w(IPR,k,j,i) = psv_(k,j,i);
    w(IDN,k,j,i) = exp((log(w(IPR,k,j,i)) - w(IDN,k,j,i))/gamma_(k,j,i));
  }

  for (int i = il; i <= iu; ++i) {
    wr(IPR,i) += psf_(k,j,i);
    wl(IPR,i+1) += psf_(k,j,i+1);
    wr(IDN,i) = exp((log(wr(IPR,i)) - wr(IDN,i))/gamma_(k,j,i));
    wl(IDN,i+1) = exp((log(wl(IPR,i+1)) - wl(IDN,i+1))/gamma_(k,j,i));
  }
}
