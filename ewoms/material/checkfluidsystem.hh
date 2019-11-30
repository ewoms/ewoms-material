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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc checkFluidSystem
 */
#ifndef EWOMS_CHECK_FLUIDSYSTEM_HH
#define EWOMS_CHECK_FLUIDSYSTEM_HH

// include all fluid systems in ewoms-material
#include <ewoms/material/fluidsystems/singlephasefluidsystem.hh>
#include <ewoms/material/fluidsystems/twophaseimmisciblefluidsystem.hh>
#include <ewoms/material/fluidsystems/h2on2fluidsystem.hh>
#include <ewoms/material/fluidsystems/h2on2liquidphasefluidsystem.hh>
#include <ewoms/material/fluidsystems/h2oairfluidsystem.hh>
#include <ewoms/material/fluidsystems/h2oairmesitylenefluidsystem.hh>
#include <ewoms/material/fluidsystems/h2oairxylenefluidsystem.hh>
#include <ewoms/material/fluidsystems/spe5fluidsystem.hh>

// include all fluid states
#include <ewoms/material/fluidstates/pressureoverlayfluidstate.hh>
#include <ewoms/material/fluidstates/saturationoverlayfluidstate.hh>
#include <ewoms/material/fluidstates/temperatureoverlayfluidstate.hh>
#include <ewoms/material/fluidstates/compositionalfluidstate.hh>
#include <ewoms/material/fluidstates/nonequilibriumfluidstate.hh>
#include <ewoms/material/fluidstates/immisciblefluidstate.hh>

#include <ewoms/common/unused.hh>
#include <dune/common/classname.hh>

#include <iostream>
#include <string>

/*!
 * \brief This is a fluid state which makes sure that only the quantities
 *        allowed are accessed.
 */
template <class ScalarT,
          class FluidSystem,
          class BaseFluidState = Ewoms::CompositionalFluidState<ScalarT, FluidSystem> >
class HairSplittingFluidState
    : protected BaseFluidState
{
public:
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    typedef ScalarT Scalar;

    HairSplittingFluidState()
    {
        // initially, do not allow anything
        allowTemperature(false);
        allowPressure(false);
        allowComposition(false);
        allowDensity(false);

        // do not allow accessing any phase
        restrictToPhase(1000);
    }

    void allowTemperature(bool yesno)
    { allowTemperature_ = yesno; }

    void allowPressure(bool yesno)
    { allowPressure_ = yesno; }

    void allowComposition(bool yesno)
    { allowComposition_ = yesno; }

    void allowDensity(bool yesno)
    { allowDensity_ = yesno; }

    void restrictToPhase(int phaseIdx)
    { restrictPhaseIdx_ = phaseIdx; }

    BaseFluidState& base()
    { return *static_cast<BaseFluidState*>(this); }

    const BaseFluidState& base() const
    { return *static_cast<const BaseFluidState*>(this); }

    auto temperature(unsigned phaseIdx) const
        -> decltype(this->base().temperature(phaseIdx))
    {
        assert(allowTemperature_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        return this->base().temperature(phaseIdx);
    }

    auto pressure(unsigned phaseIdx) const
        -> decltype(this->base().pressure(phaseIdx))
    {
        assert(allowPressure_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        return this->base().pressure(phaseIdx);
    }

    auto moleFraction(unsigned phaseIdx, unsigned compIdx) const
        -> decltype(this->base().moleFraction(phaseIdx, compIdx))
    {
        assert(allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        return this->base().moleFraction(phaseIdx, compIdx);
    }

    auto massFraction(unsigned phaseIdx, unsigned compIdx) const
        -> decltype(this->base().massFraction(phaseIdx, compIdx))
    {
        assert(allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        return this->base().massFraction(phaseIdx, compIdx);
    }

    auto averageMolarMass(unsigned phaseIdx) const
        -> decltype(this->base().averageMolarMass(phaseIdx))
    {
        assert(allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        return this->base().averageMolarMass(phaseIdx);
    }

    auto molarity(unsigned phaseIdx, unsigned compIdx) const
        -> decltype(this->base().molarity(phaseIdx, compIdx))
    {
        assert(allowDensity_ && allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        return this->base().molarity(phaseIdx, compIdx);
    }

    auto molarDensity(unsigned phaseIdx) const
        -> decltype(this->base().molarDensity(phaseIdx))
    {
        assert(allowDensity_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        return this->base().molarDensity(phaseIdx);
    }

    auto molarVolume(unsigned phaseIdx) const
        -> decltype(this->base().molarVolume(phaseIdx))
    {
        assert(allowDensity_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        return this->base().molarVolume(phaseIdx);
    }

    auto density(unsigned phaseIdx) const
        -> decltype(this->base().density(phaseIdx))
    {
        assert(allowDensity_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        return this->base().density(phaseIdx);
    }

    auto saturation(unsigned phaseIdx) const
        -> decltype(this->base().saturation(phaseIdx))
    {
        assert(false);
        return  this->base().saturation(phaseIdx);
    }

    auto fugacity(unsigned phaseIdx, unsigned compIdx) const
        -> decltype(this->base().fugacity(phaseIdx, compIdx))
    {
        assert(false);
        return this->base().fugacity(phaseIdx, compIdx);
    }

    auto fugacityCoefficient(unsigned phaseIdx, unsigned compIdx) const
        -> decltype(this->base().fugacityCoefficient(phaseIdx, compIdx))
    {
        assert(false);
        return this->base().fugacityCoefficient(phaseIdx, compIdx);
    }

    auto enthalpy(unsigned phaseIdx) const
        -> decltype(this->base().enthalpy(phaseIdx))
    {
        assert(false);
        return this->base().enthalpy(phaseIdx);
    }

    auto internalEnergy(unsigned phaseIdx) const
        -> decltype(this->base().internalEnergy(phaseIdx))
    {
        assert(false);
        return this->base().internalEnergy(phaseIdx);
    }

    auto viscosity(unsigned phaseIdx) const
        -> decltype(this->base().viscosity(phaseIdx))
    {
        assert(false);
        return this->base().viscosity(phaseIdx);
    }

private:
    bool allowSaturation_;
    bool allowTemperature_;
    bool allowPressure_;
    bool allowComposition_;
    bool allowDensity_;
    int restrictPhaseIdx_;
};

template <class Scalar, class BaseFluidState>
void checkFluidState(const BaseFluidState& fs)
{
    // fluid states must be copy-able
    BaseFluidState EWOMS_UNUSED tmpFs(fs);
    tmpFs = fs;

    // a fluid state must provide a checkDefined() method
    fs.checkDefined();

    // fluid states must export the types which they use as Scalars
    typedef typename BaseFluidState::Scalar FsScalar;
    static_assert(std::is_same<FsScalar, Scalar>::value,
                  "Fluid states must export the type they are given as scalar in an unmodified way");

    // make sure the fluid state provides all mandatory methods
    while (false) {
        Scalar val = 1.0;

        val = 2*val; // get rid of GCC warning (only occurs with paranoid warning flags)

        val = fs.temperature(/*phaseIdx=*/0);
        val = fs.pressure(/*phaseIdx=*/0);
        val = fs.moleFraction(/*phaseIdx=*/0, /*compIdx=*/0);
        val = fs.massFraction(/*phaseIdx=*/0, /*compIdx=*/0);
        val = fs.averageMolarMass(/*phaseIdx=*/0);
        val = fs.molarity(/*phaseIdx=*/0, /*compIdx=*/0);
        val = fs.molarDensity(/*phaseIdx=*/0);
        val = fs.molarVolume(/*phaseIdx=*/0);
        val = fs.density(/*phaseIdx=*/0);
        val = fs.saturation(/*phaseIdx=*/0);
        val = fs.fugacity(/*phaseIdx=*/0, /*compIdx=*/0);
        val = fs.fugacityCoefficient(/*phaseIdx=*/0, /*compIdx=*/0);
        val = fs.enthalpy(/*phaseIdx=*/0);
        val = fs.internalEnergy(/*phaseIdx=*/0);
        val = fs.viscosity(/*phaseIdx=*/0);
    };
}

/*!
 * \brief Checks whether a fluid system adheres to the specification.
 */
template <class Scalar, class FluidSystem, class RhsEval, class LhsEval>
void checkFluidSystem()
{
    std::cout << "Testing fluid system '" << Dune::className<FluidSystem>() << "'\n";

    // make sure the fluid system provides the number of phases and
    // the number of components
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    typedef HairSplittingFluidState<RhsEval, FluidSystem> FluidState;
    FluidState fs;
    fs.allowTemperature(true);
    fs.allowPressure(true);
    fs.allowComposition(true);
    fs.restrictToPhase(-1);

    // initialize memory the fluid state
    fs.base().setTemperature(273.15 + 20.0);
    for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
        fs.base().setPressure(phaseIdx, 1e5);
        fs.base().setSaturation(phaseIdx, 1.0/numPhases);
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
            fs.base().setMoleFraction(phaseIdx, compIdx, 1.0/numComponents);
        }
    }

    static_assert(std::is_same<typename FluidSystem::Scalar, Scalar>::value,
                  "The type used for floating point used by the fluid system must be the same"
                  " as the one passed to the checkFluidSystem() function");

    // check whether the parameter cache adheres to the API
    typedef typename FluidSystem::template ParameterCache<LhsEval> ParameterCache;

    ParameterCache paramCache;
    try { paramCache.updateAll(fs); } catch (...) {};
    try { paramCache.updateAll(fs, /*except=*/ParameterCache::None); } catch (...) {};
    try { paramCache.updateAll(fs, /*except=*/ParameterCache::Temperature | ParameterCache::Pressure | ParameterCache::Composition); } catch (...) {};
    try { paramCache.updateAllPressures(fs); } catch (...) {};

    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        fs.restrictToPhase(static_cast<int>(phaseIdx));
        try { paramCache.updatePhase(fs, phaseIdx); } catch (...) {};
        try { paramCache.updatePhase(fs, phaseIdx, /*except=*/ParameterCache::None); } catch (...) {};
        try { paramCache.updatePhase(fs, phaseIdx, /*except=*/ParameterCache::Temperature | ParameterCache::Pressure | ParameterCache::Composition); } catch (...) {};
        try { paramCache.updateTemperature(fs, phaseIdx); } catch (...) {};
        try { paramCache.updatePressure(fs, phaseIdx); } catch (...) {};
        try { paramCache.updateComposition(fs, phaseIdx); } catch (...) {};
        try { paramCache.updateSingleMoleFraction(fs, phaseIdx, /*compIdx=*/0); } catch (...) {};
    }

    // some value to make sure the return values of the fluid system
    // are convertible to scalars
    LhsEval val = 0.0;
    Scalar scalarVal = 0.0;

    scalarVal = 2*scalarVal; // get rid of GCC warning (only occurs with paranoid warning flags)
    val = 2*val; // get rid of GCC warning (only occurs with paranoid warning flags)

    // actually check the fluid system API
    try { FluidSystem::init(); } catch (...) {};
    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
        fs.restrictToPhase(static_cast<int>(phaseIdx));
        fs.allowPressure(FluidSystem::isCompressible(phaseIdx));
        fs.allowComposition(true);
        fs.allowDensity(false);
        try { auto tmpVal EWOMS_UNUSED = FluidSystem::density(fs, paramCache, phaseIdx); static_assert(std::is_same<decltype(tmpVal), RhsEval>::value, "The default return value must be the scalar used by the fluid state!"); } catch (...) {};
        try { val = FluidSystem::template density<FluidState, LhsEval>(fs, paramCache, phaseIdx); } catch (...) {};
        try { scalarVal = FluidSystem::template density<FluidState, Scalar>(fs, paramCache, phaseIdx); } catch (...) {};

        fs.allowPressure(true);
        fs.allowDensity(true);
        try { auto tmpVal EWOMS_UNUSED = FluidSystem::viscosity(fs, paramCache, phaseIdx); static_assert(std::is_same<decltype(tmpVal), RhsEval>::value, "The default return value must be the scalar used by the fluid state!"); } catch (...) {};
        try { auto tmpVal EWOMS_UNUSED = FluidSystem::enthalpy(fs, paramCache, phaseIdx); static_assert(std::is_same<decltype(tmpVal), RhsEval>::value, "The default return value must be the scalar used by the fluid state!"); } catch (...) {};
        try { auto tmpVal EWOMS_UNUSED = FluidSystem::heatCapacity(fs, paramCache, phaseIdx); static_assert(std::is_same<decltype(tmpVal), RhsEval>::value, "The default return value must be the scalar used by the fluid state!"); } catch (...) {};
        try { auto tmpVal EWOMS_UNUSED= FluidSystem::thermalConductivity(fs, paramCache, phaseIdx); static_assert(std::is_same<decltype(tmpVal), RhsEval>::value, "The default return value must be the scalar used by the fluid state!"); } catch (...) {};
        try { val = FluidSystem::template viscosity<FluidState, LhsEval>(fs, paramCache, phaseIdx); } catch (...) {};
        try { val = FluidSystem::template enthalpy<FluidState, LhsEval>(fs, paramCache, phaseIdx); } catch (...) {};
        try { val = FluidSystem::template heatCapacity<FluidState, LhsEval>(fs, paramCache, phaseIdx); } catch (...) {};
        try { val = FluidSystem::template thermalConductivity<FluidState, LhsEval>(fs, paramCache, phaseIdx); } catch (...) {};
        try { scalarVal = FluidSystem::template viscosity<FluidState, Scalar>(fs, paramCache, phaseIdx); } catch (...) {};
        try { scalarVal = FluidSystem::template enthalpy<FluidState, Scalar>(fs, paramCache, phaseIdx); } catch (...) {};
        try { scalarVal = FluidSystem::template heatCapacity<FluidState, Scalar>(fs, paramCache, phaseIdx); } catch (...) {};
        try { scalarVal = FluidSystem::template thermalConductivity<FluidState, Scalar>(fs, paramCache, phaseIdx); } catch (...) {};

        for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx) {
            fs.allowComposition(!FluidSystem::isIdealMixture(phaseIdx));
            try { auto tmpVal EWOMS_UNUSED = FluidSystem::fugacityCoefficient(fs, paramCache, phaseIdx, compIdx); static_assert(std::is_same<decltype(tmpVal), RhsEval>::value, "The default return value must be the scalar used by the fluid state!"); } catch (...) {};
            try { val = FluidSystem::template fugacityCoefficient<FluidState, LhsEval>(fs, paramCache, phaseIdx, compIdx); } catch (...) {};
            try { scalarVal = FluidSystem::template fugacityCoefficient<FluidState, Scalar>(fs, paramCache, phaseIdx, compIdx); } catch (...) {};
            fs.allowComposition(true);
            try { auto tmpVal EWOMS_UNUSED = FluidSystem::diffusionCoefficient(fs, paramCache, phaseIdx, compIdx); static_assert(std::is_same<decltype(tmpVal), RhsEval>::value, "The default return value must be the scalar used by the fluid state!"); } catch (...) {};
            try { val = FluidSystem::template diffusionCoefficient<FluidState, LhsEval>(fs, paramCache, phaseIdx, compIdx); } catch (...) {};
            try { scalarVal = FluidSystem::template fugacityCoefficient<FluidState, Scalar>(fs, paramCache, phaseIdx, compIdx); } catch (...) {};
        }
    }

    // test for phaseName(), isLiquid() and isIdealGas()
    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
        std::string name EWOMS_UNUSED = FluidSystem::phaseName(phaseIdx);
        bool bVal = FluidSystem::isLiquid(phaseIdx);
        bVal = FluidSystem::isIdealGas(phaseIdx);
        bVal = !bVal; // get rid of GCC warning (only occurs with paranoid warning flags)
    }

    // test for molarMass() and componentName()
    for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx) {
        val = FluidSystem::molarMass(compIdx);
        std::string name = FluidSystem::componentName(compIdx);
    }

    std::cout << "----------------------------------\n";
}

#endif
