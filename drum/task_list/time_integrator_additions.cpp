// C/C++ headers
#include <iostream>

// Athena++ headers
#include "../athena.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../chemistry/chemistry.hpp"
#include "task_list.hpp"

//----------------------------------------------------------------------------------------
// Functions for implicit correction
enum TaskStatus TimeIntegratorTaskList::UpdateHydro(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  Real dt = pmb->pmy_mesh->dt;

  if (stage <= nstages) {
    if (ph->implicit_flag)
      ph->ImplicitCorrection(ph->du, ph->w, stage_wghts[stage-1].beta*dt);
    Real wghts[3];
    wghts[0] = 1.;
    wghts[1] = 1.;
    wghts[2] = 0.;
    pmb->WeightedAve(ph->u, ph->du, ph->u2, wghts);
  }

  return TaskStatus::next;
}

//----------------------------------------------------------------------------------------
// Functions to integrate chemistry
TaskStatus TimeIntegratorTaskList::IntegrateChemistry(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;

  if (stage != nstages) return TaskStatus::next;

  // frictional heating
  pmb->pchem->AddFrictionalHeating(ph->u, ph->w, pmb->pmy_mesh->dt);

  // do slow chemistry
  pmb->pchem->EvolveOneStep(ph->u, pmb->pmy_mesh->time, pmb->pmy_mesh->dt);

  // do fast chemistry
  pmb->pthermo->SaturationAdjustment(ph->u);

  return TaskStatus::next;
}
