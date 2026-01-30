//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file PwP_Gamma_thermal.cpp
//! \brief implements functions in class EquationOfState for a piecewise polytropic EoS 
//!  with a Gamma law thermal component for temperature. Start with one polytrope.
//======================================================================================

// C headers

// C++ headers

// Athena++ headers
#include "../eos.hpp"
# include <cmath>
//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::PresFromRhoEg(Real rho, Real egas)
//! \brief Return gas pressure
Real EquationOfState::PresFromRhoEg(Real rho, Real egas) {
  Real K_c = 1.;
  Real Gamma_c = 1.4;
  Real Gamma_th = 1.8;
  return (K_c * pow(rho, Gamma_c) * (1. - (Gamma_th - 1.)/(Gamma_c - 1.))) + (egas * (Gamma_th - 1.));
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::EgasFromRhoP(Real rho, Real pres)
//! \brief Return internal energy density
Real EquationOfState::EgasFromRhoP(Real rho, Real pres) {
  Real K_c = 1.;
  Real Gamma_c = 1.4;
  Real Gamma_th = 1.8;
  return (pres - (K_c * pow(rho, Gamma_c) * (1. - (Gamma_th - 1.)/(Gamma_c - 1.))))/(Gamma_th - 1.);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::AsqFromRhoP(Real rho, Real pres)
//! \brief Return adiabatic sound speed squared
Real EquationOfState::AsqFromRhoP(Real rho, Real pres) {
  Real K_c = 1.;
  Real Gamma_c = 1.4;
  Real Gamma_th = 1.8;
  return Gamma_th * pres / rho + K_c * (Gamma_c - Gamma_th) * pow(rho, Gamma_c - 1.);
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::InitEosConstants(ParameterInput* pin)
//! \brief Initialize constants for EOS
void EquationOfState::InitEosConstants(ParameterInput *pin) {
  return;
}
