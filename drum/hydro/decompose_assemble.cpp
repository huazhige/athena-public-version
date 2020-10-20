// C/C++ headers
#include <sstream>

// Athena++ headers
#include "hydro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../coordinates/coordinates.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../globals.hpp"
#include "../reconstruct/interpolation.hpp"
#include "../utils/utils.hpp"

inline void IntegrateUpwards(AthenaArray<Real>& psf, AthenaArray<Real> const& w, Coordinates *pco,
  Real grav, int kl, int ku, int jl, int ju, int il, int iu)
{
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i) {
        psf(k,j,i+1) = psf(k,j,i) - grav*w(IDN,k,j,i)*pco->dx1f(i);
        if (psf(k,j,i+1) < 0.)
          psf(k,j,i+1) = psf(k,j,i)*exp(-grav*w(IDN,k,j,i)*pco->dx1f(i)/w(IPR,k,j,i));
      }

}

inline void IntegrateDownwards(AthenaArray<Real>& psf, AthenaArray<Real> const& w, Coordinates *pco,
  Real grav, int kl, int ku, int jl, int ju, int il, int iu)
{
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = iu; i >= il; --i)
        psf(k,j,i) = psf(k,j,i+1) + grav*w(IDN,k,j,i)*pco->dx1f(i);
}

void Hydro::DecomposePressure(AthenaArray<Real> &w, int kl, int ku, int jl, int ju)
{
  // Need to integrate upwards
  MeshBlock *pmb = pmy_block;
  Coordinates *pco = pmb->pcoord;
  Thermodynamics *pthermo = pmb->pthermo;

  Real grav = -hsrc.GetG1();  // positive downward pointing
  Real Rd = pthermo->GetRd();

  int is = pmb->is, ie = pmb->ie;
  
  if (grav == 0.) return;

  std::stringstream msg;
  if (NGHOST < 3) {
    msg << "### FATAL ERROR in function [Hydro::DecomposePressure]"
        << std::endl << "number of ghost cells (NGHOST) must be at least 3" <<std::endl;
    ATHENA_ERROR(msg);
  }


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

  Real **w1;
  NewCArray(w1, 2, NHYDRO);

  if (has_bot_neighbor) {
    RecvBotPressure(psf_, entropy_, gamma_, nbot, kl, ku, jl, ju);
  } else {
    // bottom layer polytropic index
    pthermo->PolytropicIndex(gamma_, w, kl, ku, jl, ju, is, is);
    // bottom layer entropy
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        entropy_(k,j,is) = log(w(IPR,k,j,is)) - gamma_(k,j,is)*log(w(IDN,k,j,is));

    if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
      // adiabatic extrapolation to bottom boundary
      for (int k = kl; k <= ku; ++k)
        for (int j = jl; j <= ju; ++j) {
          Real P1 = w(IPR,k,j,is-1);
          Real T1 = pthermo->Temp(w.at(k,j,is-1));
          Real dz = pco->dx1f(is-1);
          for (int n = 0; n < NMASS; ++n)
            w1[0][n] = w(n,k,j,is-1);
          pthermo->ConstructAdiabat(w1, T1, P1, grav, dz/2., 2, Adiabat::reversible);
          psf_(k,j,is) = w1[1][IPR];
          for (int n = 0; n < NHYDRO; ++n)
            hydro_face_(n,k,j,is) = w1[1][n];
        }
    } else {  // reflecing boundary condition or else
      // polynomical interpolation to find the pressure at bottom boundary
      for (int k = kl; k <= ku; ++k)
        for (int j = jl; j <= ju; ++j)
          psf_(k,j,is) = exp(interp_cp6(log(w(IPR,k,j,is-3)), log(w(IPR,k,j,is-2)), 
            log(w(IPR,k,j,is-1)), log(w(IPR,k,j,is)), log(w(IPR,k,j,is+1)), log(w(IPR,k,j,is+2))));
    }
  }
  IntegrateUpwards(psf_, w, pco, grav, kl, ku, jl, ju, is, ie);
  
  if (has_top_neighbor)
    SendBotPressure(psf_, entropy_, gamma_, ntop, kl, ku, jl, ju);

  // integrate ghost cells
  if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::reflect)
    IntegrateDownwards(psf_, w, pco, -grav, kl, ku, jl, ju, is - NGHOST, is - 1);
  else  // block boundary
    IntegrateDownwards(psf_, w, pco,  grav, kl, ku, jl, ju, is - NGHOST, is - 1);

  if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::reflect)
    IntegrateUpwards(psf_, w, pco, -grav, kl, ku, jl, ju, ie + 1, ie + NGHOST);
  else  // block boundary
    IntegrateUpwards(psf_, w, pco,  grav, kl, ku, jl, ju, ie + 1, ie + NGHOST);

  // decompose pressure
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j) {
      // 1. change density and pressure (including ghost cells)
      for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
        // interpolate hydrostatic pressure, prevent divided by zero
        if (fabs(psf_(k,j,i) - psf_(k,j,i+1)) < 1.E-6)
          psv_(k,j,i) = (psf_(k,j,i) + psf_(k,j,i+1))/2.;
        else
          psv_(k,j,i) = (psf_(k,j,i) - psf_(k,j,i+1))/log(psf_(k,j,i)/psf_(k,j,i+1));

        // save reference adiabatic density profile
        dsv_(k,j,i) = pow(psv_(k,j,i), 1./gamma_(k,j,is))*exp(-entropy_(k,j,is)/gamma_(k,j,is));

        // change pressure/density to pertubation quantities
        w(IPR,k,j,i) -= psv_(k,j,i);
        w(IDN,k,j,i) -= dsv_(k,j,i);
      }

      // 2. override outflow bottom boundary
      if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
        Real p1 = w(IPR,k,j,is);
        Real p2 = w(IPR,k,j,is+1);
        Real p3 = w(IPR,k,j,is+2);
        // is-1
        psv_(k,j,is-1) += w(IPR,k,j,is-1);
        w(IPR,k,j,is-1) = inflection3_cell1(p1, p2, p3);
        psv_(k,j,is-1) -= w(IPR,k,j,is-1);

        // is-2
        psv_(k,j,is-2) += w(IPR,k,j,is-2);
        w(IPR,k,j,is-2) = inflection3_cell2(p1, p2, p3);
        psv_(k,j,is-2) -= w(IPR,k,j,is-2);

        // is-3
        psv_(k,j,is-3) += w(IPR,k,j,is-3);
        w(IPR,k,j,is-3) = inflection3_cell3(p1, p2, p3);
        psv_(k,j,is-3) -= w(IPR,k,j,is-3);
      }

      // 3. override outflow top boundary
      if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::outflow) {
        Real p1 = w(IPR,k,j,ie);
        Real p2 = w(IPR,k,j,ie-1);
        Real p3 = w(IPR,k,j,ie-2);
        // ie+1
        psv_(k,j,ie+1) += w(IPR,k,j,ie+1);
        w(IPR,k,j,ie+1) = inflection3_cell1(p1, p2, p3);
        psv_(k,j,ie+1) -= w(IPR,k,j,ie+1);

        // ie+2
        psv_(k,j,ie+2) += w(IPR,k,j,ie+2);
        w(IPR,k,j,ie+2) = inflection3_cell2(p1, p2, p3);
        psv_(k,j,ie+2) -= w(IPR,k,j,ie+2);
        
        // ie+3
        psv_(k,j,ie+3) += w(IPR,k,j,ie+3);
        w(IPR,k,j,ie+3) = inflection3_cell3(p1, p2, p3);
        psv_(k,j,ie+3) -= w(IPR,k,j,ie+3);
      }
    }

  /* debug
  if (Globals::my_rank == 0) {
    std::cout << "========== " << jl << std::endl;
    for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
      if (i == is)
        std::cout << "-------- ";
      if (i == 0)
        std::cout << "i = " << "-1/2 ";
      else if (i == 1)
        std::cout << "i = " << "+1/2 ";
      else
        std::cout << "i = " << i-1 << "+1/2 ";
      std::cout << "psf = " << psf_(kl,jl,i) << std::endl;
      std::cout << "i = " << i  << "    ";
      //std::cout << "pre = " << w(IPR,kl,jl,i) << std::endl;
      std::cout << "psv = " << psv_(kl,jl,i) + w(IPR,kl,jl,i) << " pre = " << w(IPR,kl,jl,i) 
                << " dsv = " << dsv_(kl,jl,i) << " den = " << w(IDN,kl,jl,i) << std::endl;
      if (i == ie)
        std::cout << "-------- ";
      if (i == ie + NGHOST) {
        std::cout << "i = " << i+1 << "+1/2 ";
        std::cout << "psf = " << psf_(kl,jl,i+1) << std::endl;
      }
    }
    std::cout << "==========" << std::endl;
  }*/

  // finish send top pressure
  if (has_top_neighbor && (ntop.snb.rank != Globals::my_rank))
    WaitBotPressure();

  FreeCArray(w1);
}

void Hydro::AssemblePressure(AthenaArray<Real> &w,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr, int k, int j, int il, int iu)
{
  MeshBlock *pmb = pmy_block;
  int is = pmb->is, ie = pmb->ie;
  if (hsrc.GetG1() == 0.) return;
  
  for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
    w(IPR,k,j,i) += psv_(k,j,i);
    w(IDN,k,j,i) += dsv_(k,j,i);
  }

  for (int i = il; i <= iu; ++i) {
    wr(IPR,i) += psf_(k,j,i);
    wl(IPR,i+1) += psf_(k,j,i+1);
    wr(IDN,i) += pow(psf_(k,j,i), 1./gamma_(k,j,is))*exp(-entropy_(k,j,is)/gamma_(k,j,is));
    wl(IDN,i+1) += pow(psf_(k,j,i+1), 1./gamma_(k,j,is))*exp(-entropy_(k,j,is)/gamma_(k,j,is));
  }

  // override outflow bottom boundary
  if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
    wl(IPR,is) = psf_(k,j,is);
    wl(IDN,is) = hydro_face_(IDN,k,j,is);
  }
}
