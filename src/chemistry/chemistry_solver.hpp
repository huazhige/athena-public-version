#ifndef CHEMISTRY_SOLVER_HPP
#define CHEMISTRY_SOLVER_HPP

#include "../math/linalg.h"

class ChemistrySolver {
public:
  ChemistrySolver() {}

  // solve independent chemistry: A*dq = b
  // the second argument Real *dq is both input (b) and output (dq)
  void SolveIndependent(Real **A, Real *dq) {
    for (int i = 0; i < NVAPOR; ++i) {
      // assemble matrix for each species
      #pragma ivdep
      for (int j = 0; j < NPHASE; ++j) {
        for (int k = 0; k < NPHASE; ++k)
          a_[j][k] = A[1+j*NVAPOR+i][1+k*NVAPOR+i];
        sol_[j] = dq[1+j*NVAPOR+i];
      }

      ludcmp<NPHASE>(a_, indx_, vv_);
      lubksb<NPHASE>(a_, indx_, sol_);

      for (int j = 0; j < NPHASE; ++j)
        dq[1+j*NVAPOR+i] = sol_[j];
    }
    dq[0] = 0.;
  }

private:
  // scratch arrays
  Real a_[NPHASE][NPHASE];
  Real sol_[NPHASE], vv_[NPHASE];
  int indx_[NPHASE];
};

#endif
