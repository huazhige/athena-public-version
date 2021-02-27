// C/C++ header
#include <ctime>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../math/interpolation.h"
#include "../globals.hpp"
#include "../utils/utils.hpp"
#include "../thermodynamics/thermodynamics.hpp"

// molecules
enum {iH2Ov = 1,  iH2Oc = 2,  iH2Op = 3}; 

Real P0, T0;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(5);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "thetav");
  SetUserOutputVariableName(3, "mse");
  SetUserOutputVariableName(4, "rh_h2o");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  Real dq[1+NVAPOR], rh;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0,k,j,i) = pthermo->Temp(phydro->w.at(k,j,i));
        user_out_var(1,k,j,i) = pthermo->Theta(phydro->w.at(k,j,i), P0);
        // theta_v
        user_out_var(2,k,j,i) = user_out_var(1,k,j,i)*pthermo->Qeps(phydro->w.at(k,j,i));
        user_out_var(3,k,j,i) = pthermo->MSE(phydro->w.at(k,j,i), grav*pcoord->x1v(i));
        pthermo->SaturationSurplus(dq, phydro->w.at(k,j,i));
        rh = phydro->w(iH2Ov,k,j,i)/(phydro->w(iH2Ov,k,j,i) - dq[iH2Ov]);
        user_out_var(4,k,j,i) = rh;
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real qH2O = pin->GetReal("problem", "qH2O");
  Real gamma = pin->GetReal("hydro", "gamma");

  // estimate surface temperature and pressure
  Real Rd = pthermo->GetRd();
  Real cp = gamma/(gamma - 1.)*Rd;
  Real temp = T0;
  Real R = Rd*( qH2O/pthermo->GetMassRatio(iH2Ov) +
                (1. - qH2O) );
  Real rho0 = P0/(R*temp);
  // initial condition
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        phydro->w(IPR,k,j,i) = P0;
        phydro->w(iH2Ov,k,j,i) = qH2O;
        phydro->w(IDN,k,j,i) = rho0;
        phydro->w(IV1,k,j,i) = 0.;
        phydro->w(IV2,k,j,i) = 0.;
        phydro->w(IV3,k,j,i) = 0.;

      }
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
