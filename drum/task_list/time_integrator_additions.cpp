// Athena++ headers
#include "../athena.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "task_list.hpp"

//----------------------------------------------------------------------------------------
// Functions to integrate chemistry
TaskStatus TimeIntegratorTaskList::IntegrateChemistry(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;

  if (stage != nstages) return TaskStatus::next;

  // do fast chemistry
  pmb->pthermo->SaturationAdjustment(ph->u);

  return TaskStatus::next;
}
