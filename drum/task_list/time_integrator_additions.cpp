// C/C++ headers
#include <iostream>

// Athena++ headers
#include "../athena.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../chemistry/chemistry.hpp"
#include "task_list.hpp"

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
