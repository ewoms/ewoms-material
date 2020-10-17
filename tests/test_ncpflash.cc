// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the eWoms project.

  eWoms is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  eWoms is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with eWoms.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 *
 * \brief This is a program to test the flash calculation which uses
 *        non-linear complementarity problems (NCP)
 *
 * A flash calculation determines the pressures, saturations and
 * composition of all phases given the total mass (or, as in this case
 * the total number of moles) in a given amount of pore space.
 */
#include "config.h"

#include <ewoms/material/constraintsolvers/ncpflash.hh>
#include <ewoms/material/constraintsolvers/misciblemultiphasecomposition.hh>
#include <ewoms/material/constraintsolvers/computefromreferencephase.hh>

#include <ewoms/material/fluidstates/compositionalfluidstate.hh>

#include <ewoms/material/fluidsystems/h2on2fluidsystem.hh>

#include <ewoms/material/fluidmatrixinteractions/linearmaterial.hh>
#include <ewoms/material/fluidmatrixinteractions/regularizedbrookscorey.hh>
#include <ewoms/material/fluidmatrixinteractions/efftoabslaw.hh>
#include <ewoms/material/fluidmatrixinteractions/materialtraits.hh>

#include <dune/common/parallel/mpihelper.hh>

template <class Scalar, class FluidState>
void checkSame(const FluidState& fsRef, const FluidState& fsFlash)
{
    enum { numPhases = FluidState::numPhases };
    enum { numComponents = FluidState::numComponents };

    Scalar tol = std::max(std::numeric_limits<Scalar>::epsilon()*1e4, 1e-6);

    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        Scalar error;

        // check the pressures
        error = 1 - fsRef.pressure(phaseIdx)/fsFlash.pressure(phaseIdx);
        if (std::abs(error) > tol) {
            std::ostringstream oss;
            oss << "pressure error for phase " << phaseIdx << " exceeds tolerance"
                << " (" << fsFlash.pressure(phaseIdx)  << " flash vs "
                << fsRef.pressure(phaseIdx) << " reference,"
                << " error=" << error << ")";
            throw std::runtime_error(oss.str());
        }

        // check the saturations
        error = fsRef.saturation(phaseIdx) - fsFlash.saturation(phaseIdx);
        if (std::abs(error) > tol) {
            std::ostringstream oss;
            oss << "saturation error for phase " << phaseIdx << " exceeds tolerance"
                << " (" << fsFlash.saturation(phaseIdx) << " flash vs "
                << fsRef.saturation(phaseIdx) << " reference,"
                << " error=" << error << ")";
            throw std::runtime_error(oss.str());
        }

        // check the compositions
        for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx) {
            error = fsRef.moleFraction(phaseIdx, compIdx) - fsFlash.moleFraction(phaseIdx, compIdx);
            if (std::abs(error) > tol) {
                std::ostringstream oss;
                oss << "composition error phase " << phaseIdx << ", component " << compIdx << " exceeds tolerance"
                    << " (" << fsFlash.moleFraction(phaseIdx, compIdx) << " flash vs "
                    << fsRef.moleFraction(phaseIdx, compIdx) << " reference,"
                    << " error=" << error << ")";
                throw std::runtime_error(oss.str());
            }
        }
    }
}

template <class Scalar, class FluidSystem, class MaterialLaw, class FluidState>
void checkNcpFlash(const FluidState& fsRef,
                   typename MaterialLaw::Params& matParams)
{
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;
    typedef typename FluidSystem::template ParameterCache<typename FluidState::Scalar> ParameterCache;

    // calculate the total amount of stuff in the reference fluid
    // phase
    ComponentVector globalMolarities(0.0);
    for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            globalMolarities[compIdx] +=
                fsRef.saturation(phaseIdx)*fsRef.molarity(phaseIdx, compIdx);
        }
    }

    // initialize the fluid state for the flash calculation
    typedef Ewoms::NcpFlash<Scalar, FluidSystem> NcpFlash;
    FluidState fsFlash;

    fsFlash.setTemperature(fsRef.temperature(/*phaseIdx=*/0));

    // run the flash calculation
    ParameterCache paramCache;
    paramCache.updateAll(fsFlash);
    NcpFlash::guessInitial(fsFlash, globalMolarities);
    NcpFlash::template solve<MaterialLaw>(fsFlash, matParams, paramCache, globalMolarities);

    // compare the "flashed" fluid state with the reference one
    checkSame<Scalar>(fsRef, fsFlash);
}

template <class Scalar, class FluidSystem, class MaterialLaw, class FluidState>
void completeReferenceFluidState(FluidState& fs,
                                 typename MaterialLaw::Params& matParams,
                                 unsigned refPhaseIdx)
{
    enum { numPhases = FluidSystem::numPhases };

    typedef Ewoms::ComputeFromReferencePhase<Scalar, FluidSystem> ComputeFromReferencePhase;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

    unsigned otherPhaseIdx = 1 - refPhaseIdx;

    // calculate the other saturation
    fs.setSaturation(otherPhaseIdx, 1.0 - fs.saturation(refPhaseIdx));

    // calulate the capillary pressure
    PhaseVector pC;
    MaterialLaw::capillaryPressures(pC, matParams, fs);
    fs.setPressure(otherPhaseIdx,
                   fs.pressure(refPhaseIdx)
                   + (pC[otherPhaseIdx] - pC[refPhaseIdx]));

    // make the fluid state consistent with local thermodynamic
    // equilibrium
    typename FluidSystem::template ParameterCache<typename FluidState::Scalar> paramCache;
    ComputeFromReferencePhase::solve(fs,
                                     paramCache,
                                     refPhaseIdx,
                                     /*setViscosity=*/false,
                                     /*setEnthalpy=*/false);
}

template <class Scalar>
inline void testAll()
{
    typedef Ewoms::H2ON2FluidSystem<Scalar> FluidSystem;
    typedef Ewoms::CompositionalFluidState<Scalar, FluidSystem> CompositionalFluidState;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

    enum { H2OIdx = FluidSystem::H2OIdx };
    enum { N2Idx = FluidSystem::N2Idx };

    typedef Ewoms::TwoPhaseMaterialTraits<Scalar, liquidPhaseIdx, gasPhaseIdx> MaterialTraits;
    typedef Ewoms::RegularizedBrooksCorey<MaterialTraits> EffMaterialLaw;
    typedef Ewoms::EffToAbsLaw<EffMaterialLaw> MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    std::cout << "---- using " << Dune::className<Scalar>() << " as scalar ----\n";
    Scalar T = 273.15 + 25;

    // initialize the tables of the fluid system
    Scalar Tmin = T - 1.0;
    Scalar Tmax = T + 1.0;
    unsigned nT = 3;

    Scalar pmin = 0.0;
    Scalar pmax = 1.25 * 2e6;
    unsigned np = 100;

    FluidSystem::init(Tmin, Tmax, nT, pmin, pmax, np);

    // set the parameters for the capillary pressure law
    MaterialLawParams matParams;
    matParams.setResidualSaturation(MaterialLaw::wettingPhaseIdx, 0.0);
    matParams.setResidualSaturation(MaterialLaw::nonWettingPhaseIdx, 0.0);
    matParams.setEntryPressure(0);
    matParams.setLambda(2.0);
    matParams.finalize();

    CompositionalFluidState fsRef;

    // create an fluid state which is consistent

    // set the fluid temperatures
    fsRef.setTemperature(T);

    ////////////////
    // only liquid
    ////////////////
    std::cout << "testing single-phase liquid\n";

    // set liquid saturation
    fsRef.setSaturation(liquidPhaseIdx, 1.0);

    // set pressure of the liquid phase
    fsRef.setPressure(liquidPhaseIdx, 2e5);

    // set the liquid composition to pure water
    fsRef.setMoleFraction(liquidPhaseIdx, N2Idx, 0.0);
    fsRef.setMoleFraction(liquidPhaseIdx, H2OIdx, 1.0 - fsRef.moleFraction(liquidPhaseIdx, N2Idx));

    // "complete" the fluid state
    completeReferenceFluidState<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams, liquidPhaseIdx);

    // check the flash calculation
    checkNcpFlash<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams);

    ////////////////
    // only gas
    ////////////////
    std::cout << "testing single-phase gas\n";
    // set gas saturation
    fsRef.setSaturation(gasPhaseIdx, 1.0);

    // set pressure of the gas phase
    fsRef.setPressure(gasPhaseIdx, 1e6);

    // set the gas composition to 99.9% nitrogen and 0.1% water
    fsRef.setMoleFraction(gasPhaseIdx, N2Idx, 0.999);
    fsRef.setMoleFraction(gasPhaseIdx, H2OIdx, 0.001);

    // "complete" the fluid state
    completeReferenceFluidState<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams, gasPhaseIdx);

    // check the flash calculation
    checkNcpFlash<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams);

    ////////////////
    // both phases
    ////////////////
    std::cout << "testing two-phase\n";

    // set saturations
    fsRef.setSaturation(liquidPhaseIdx, 0.5);
    fsRef.setSaturation(gasPhaseIdx, 0.5);

    // set pressures
    fsRef.setPressure(liquidPhaseIdx, 1e6);
    fsRef.setPressure(gasPhaseIdx, 1e6);

    typename FluidSystem::template ParameterCache<Scalar> paramCache;
    typedef Ewoms::MiscibleMultiPhaseComposition<Scalar, FluidSystem> MiscibleMultiPhaseComposition;
    MiscibleMultiPhaseComposition::solve(fsRef, paramCache,
                                         /*setViscosity=*/false,
                                         /*setEnthalpy=*/false);

    // check the flash calculation
    checkNcpFlash<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams);

    ////////////////
    // with capillary pressure
    ////////////////
    MaterialLawParams matParams2;
    matParams2.setResidualSaturation(MaterialLaw::wettingPhaseIdx, 0.0);
    matParams2.setResidualSaturation(MaterialLaw::nonWettingPhaseIdx, 0.0);
    matParams2.setEntryPressure(1e3);
    matParams2.setLambda(2.0);
    matParams2.finalize();

    // set gas saturation
    fsRef.setSaturation(gasPhaseIdx, 0.5);
    fsRef.setSaturation(liquidPhaseIdx, 0.5);

    // set pressure of the liquid phase
    fsRef.setPressure(liquidPhaseIdx, 1e6);

    // calulate the capillary pressure
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;
    PhaseVector pC;
    MaterialLaw::capillaryPressures(pC, matParams2, fsRef);
    fsRef.setPressure(gasPhaseIdx,
                      fsRef.pressure(liquidPhaseIdx)
                      + (pC[gasPhaseIdx] - pC[liquidPhaseIdx]));

    typedef Ewoms::MiscibleMultiPhaseComposition<Scalar, FluidSystem> MiscibleMultiPhaseComposition;
    MiscibleMultiPhaseComposition::solve(fsRef, paramCache,
                                         /*setViscosity=*/false,
                                         /*setEnthalpy=*/false);

    // check the flash calculation
    checkNcpFlash<Scalar, FluidSystem, MaterialLaw>(fsRef, matParams2);
}

int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);

    testAll<double>();
    testAll<float>();

    return 0;
}
