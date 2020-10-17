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
 * \brief This test makes sure that the API for fluid-matrix
 *        interactions is observed by all capillary pressure / relperm
 *        laws.
 */
#include "config.h"

// include the local AD framwork
#include <ewoms/common/densead/math.hh>

// include all capillary pressure laws
#include <ewoms/material/fluidmatrixinteractions/brookscorey.hh>
#include <ewoms/material/fluidmatrixinteractions/parkerlenhard.hh>
#include <ewoms/material/fluidmatrixinteractions/linearmaterial.hh>
#include <ewoms/material/fluidmatrixinteractions/vangenuchten.hh>
#include <ewoms/material/fluidmatrixinteractions/nullmaterial.hh>
#include <ewoms/material/fluidmatrixinteractions/regularizedbrookscorey.hh>
#include <ewoms/material/fluidmatrixinteractions/regularizedvangenuchten.hh>
#include <ewoms/material/fluidmatrixinteractions/efftoabslaw.hh>
#include <ewoms/material/fluidmatrixinteractions/piecewiselineartwophasematerial.hh>
#include <ewoms/material/fluidmatrixinteractions/splinetwophasematerial.hh>
#include <ewoms/material/fluidmatrixinteractions/threephaseparkervangenuchten.hh>
#include <ewoms/material/fluidmatrixinteractions/eclepstwophaselaw.hh>
#include <ewoms/material/fluidmatrixinteractions/eclhysteresistwophaselaw.hh>
#include <ewoms/material/fluidmatrixinteractions/ecldefaultmaterial.hh>
#include <ewoms/material/fluidmatrixinteractions/eclstone1material.hh>
#include <ewoms/material/fluidmatrixinteractions/eclstone2material.hh>
#include <ewoms/material/fluidmatrixinteractions/ecltwophasematerial.hh>
#include <ewoms/material/fluidmatrixinteractions/eclmultiplexermaterial.hh>

// include the helper classes to construct traits
#include <ewoms/material/fluidmatrixinteractions/materialtraits.hh>

// include some fluid states
#include <ewoms/material/fluidstates/compositionalfluidstate.hh>
#include <ewoms/material/fluidstates/immisciblefluidstate.hh>

// include some fluid systems
#include <ewoms/material/fluidsystems/twophaseimmisciblefluidsystem.hh>
#include <ewoms/material/fluidsystems/blackoilfluidsystem.hh>

// include some components
#include <ewoms/material/components/simpleh2o.hh>
#include <ewoms/material/components/n2.hh>

#include <ewoms/common/unused.hh>

#include <dune/common/parallel/mpihelper.hh>

// this function makes sure that a capillary pressure law adheres to
// the generic programming interface for such laws. This API _must_ be
// implemented by all capillary pressure laws. If there are no _very_
// strong reasons to do otherwise, numerical models should only use on
// this API.
template <class MaterialLaw, class FluidState>
void testGenericApi()
{
    while (0) {
        // ensure the presence of the required values
        static const int numPhases = MaterialLaw::numPhases;

        // check for the presence of the is*Dependent values
        static const bool EWOMS_UNUSED isSaturationDependent = MaterialLaw::isSaturationDependent;
        static const bool EWOMS_UNUSED isPressureDependent = MaterialLaw::isPressureDependent;
        static const bool EWOMS_UNUSED isTemperatureDependent = MaterialLaw::isTemperatureDependent;
        static const bool EWOMS_UNUSED isCompositionDependent = MaterialLaw::isCompositionDependent;

        // Make sure that the Traits, Params and Scalar typedefs are
        // exported by the material law
        typedef typename MaterialLaw::Params Params;
        typedef typename MaterialLaw::Traits Traits;
        typedef typename MaterialLaw::Scalar Scalar;
        typedef typename MaterialLaw::Traits::Scalar TraitsScalar;

        static_assert(std::is_same<Scalar, TraitsScalar>::value,
                      "The traits and the material law must use the same type as scalar value");
        static_assert(numPhases == Traits::numPhases,
                      "The traits and the material law must use the number of fluid phases");

        // check the API of the parameter class. setting the actual
        // parameter values is implementation specific. But all
        // parameters must be default and copy constructible as well
        // as exhibit the finalize() method!
        Params params;
        params.finalize();
        const Params paramsConst(params);

        // test the generic methods which need to be implemented by
        // all material laws
        const FluidState fs;

        {
            double destValues[numPhases];
            MaterialLaw::capillaryPressures(destValues, paramsConst, fs);
            MaterialLaw::saturations(destValues, paramsConst, fs);
            MaterialLaw::relativePermeabilities(destValues, paramsConst, fs);
        }

        {
            typename FluidState::Scalar destValuesEval[numPhases];
            MaterialLaw::capillaryPressures(destValuesEval, paramsConst, fs);
            MaterialLaw::saturations(destValuesEval, paramsConst, fs);
            MaterialLaw::relativePermeabilities(destValuesEval, paramsConst, fs);
        }
    }
}

// this function makes ensures that a pressure law adheres to the
// covenience programming interface for two-phase material laws. The
// main purpose of this interface is to simplify the implementation of
// nested material laws.
template <class MaterialLaw, class FluidState>
void testTwoPhaseApi()
{
    typedef typename MaterialLaw::Scalar Scalar;

    while (0) {
        static const int numPhases = MaterialLaw::numPhases;
        static_assert(numPhases == 2,
                      "The number of fluid phases for a twophase "
                      "capillary pressure law must be 2");
        static_assert(MaterialLaw::implementsTwoPhaseApi,
                      "This material law is expected to implement "
                      "the two-phase API!");

        static const int EWOMS_UNUSED wettingPhaseIdx = MaterialLaw::wettingPhaseIdx;
        static const int EWOMS_UNUSED nonWettingPhaseIdx = MaterialLaw::nonWettingPhaseIdx;

        // make sure the two-phase specific methods are present
        const FluidState fs;
        const typename MaterialLaw::Params params;

        Scalar v EWOMS_UNUSED;
        v = MaterialLaw::template pcnw<FluidState, Scalar>(params, fs);
        v = MaterialLaw::template Sw<FluidState, Scalar>(params, fs);
        v = MaterialLaw::template Sn<FluidState, Scalar>(params, fs);
        v = MaterialLaw::template krw<FluidState, Scalar>(params, fs);
        v = MaterialLaw::template krn<FluidState, Scalar>(params, fs);

        typename FluidState::Scalar vEval EWOMS_UNUSED;
        vEval = MaterialLaw::pcnw(params, fs);
        vEval = MaterialLaw::Sw(params, fs);
        vEval = MaterialLaw::Sn(params, fs);
        vEval = MaterialLaw::krw(params, fs);
        vEval = MaterialLaw::krn(params, fs);

    }
}

template <class MaterialLaw, class FluidState>
void testTwoPhaseSatApi()
{
    typedef typename MaterialLaw::Scalar Scalar;

    while (0) {
        static_assert(MaterialLaw::implementsTwoPhaseSatApi,
                      "This material law is expected to implement "
                      "the two-phase saturation only API!");
        static_assert(!MaterialLaw::isPressureDependent,
                      "Capillary pressure laws which implement the twophase saturation only "
                      "API cannot be dependent on the absolute phase pressures!");
        static_assert(!MaterialLaw::isTemperatureDependent,
                      "Capillary pressure laws which implement the twophase saturation only "
                      "API cannot be dependent on temperature!");
        static_assert(!MaterialLaw::isCompositionDependent,
                      "Capillary pressure laws which implement the twophase saturation only "
                      "API cannot be dependent on the phase compositions!");

        static const int EWOMS_UNUSED numPhases = MaterialLaw::numPhases;

        // make sure the two-phase specific methods are present
        const typename MaterialLaw::Params params;

        Scalar Sw = 0;
        Scalar v EWOMS_UNUSED;
        v = MaterialLaw::twoPhaseSatPcnw(params, Sw);
        v = MaterialLaw::twoPhaseSatSw(params, Sw);
        v = MaterialLaw::twoPhaseSatSn(params, Sw);
        v = MaterialLaw::twoPhaseSatKrw(params, Sw);
        v = MaterialLaw::twoPhaseSatKrn(params, Sw);

        typename FluidState::Scalar SwEval = 0;
        typename FluidState::Scalar vEval EWOMS_UNUSED;
        vEval = MaterialLaw::twoPhaseSatPcnw(params, SwEval);
        vEval = MaterialLaw::twoPhaseSatSw(params, SwEval);
        vEval = MaterialLaw::twoPhaseSatSn(params, SwEval);
        vEval = MaterialLaw::twoPhaseSatKrw(params, SwEval);
        vEval = MaterialLaw::twoPhaseSatKrn(params, SwEval);
    }
}

template <class MaterialLaw, class FluidState>
void testThreePhaseApi()
{
    typedef typename MaterialLaw::Scalar Scalar;

    while (0) {
        static const int numPhases = MaterialLaw::numPhases;
        static_assert(numPhases == 3,
                      "The number of fluid phases for a threephase "
                      "capillary pressure law must be 3");

        static const int EWOMS_UNUSED wettingPhaseIdx = MaterialLaw::wettingPhaseIdx;
        static const int EWOMS_UNUSED nonWettingPhaseIdx = MaterialLaw::nonWettingPhaseIdx;
        static const int EWOMS_UNUSED gasPhaseIdx = MaterialLaw::gasPhaseIdx;

        // make sure the two-phase specific methods are present
        const FluidState fs;
        const typename MaterialLaw::Params params;

        Scalar v EWOMS_UNUSED;
        v = MaterialLaw::template pcnw<FluidState, Scalar>(params, fs);
        v = MaterialLaw::template Sw<FluidState, Scalar>(params, fs);
        v = MaterialLaw::template Sn<FluidState, Scalar>(params, fs);
        v = MaterialLaw::template Sg<FluidState, Scalar>(params, fs);
        v = MaterialLaw::template krw<FluidState, Scalar>(params, fs);
        v = MaterialLaw::template krn<FluidState, Scalar>(params, fs);
        v = MaterialLaw::template krg<FluidState, Scalar>(params, fs);

        typename FluidState::Scalar vEval EWOMS_UNUSED;
        vEval = MaterialLaw::pcnw(params, fs);
        vEval = MaterialLaw::Sw(params, fs);
        vEval = MaterialLaw::Sn(params, fs);
        vEval = MaterialLaw::Sg(params, fs);
        vEval = MaterialLaw::krw(params, fs);
        vEval = MaterialLaw::krn(params, fs);
        vEval = MaterialLaw::krg(params, fs);
    }
}

template <class MaterialLaw>
void testThreePhaseSatApi()
{
}

template <class Scalar>
inline void testAll()
{
    typedef Ewoms::SimpleH2O<Scalar> H2O;
    typedef Ewoms::N2<Scalar> N2;

    typedef Ewoms::LiquidPhase<Scalar, H2O> Liquid;
    typedef Ewoms::GasPhase<Scalar, N2> Gas;

    typedef Ewoms::TwoPhaseImmiscibleFluidSystem<Scalar, Liquid, Gas> TwoPFluidSystem;
    typedef Ewoms::BlackOilFluidSystem<Scalar> ThreePFluidSystem;

    typedef Ewoms::TwoPhaseMaterialTraits<Scalar,
                                        TwoPFluidSystem::wettingPhaseIdx,
                                        TwoPFluidSystem::nonWettingPhaseIdx> TwoPhaseTraits;

    typedef Ewoms::ThreePhaseMaterialTraits<Scalar,
                                          ThreePFluidSystem::waterPhaseIdx,
                                          ThreePFluidSystem::oilPhaseIdx,
                                          ThreePFluidSystem::gasPhaseIdx> ThreePhaseTraits;

    typedef Ewoms::DenseAd::Evaluation<Scalar, 3> Evaluation;
    typedef Ewoms::ImmiscibleFluidState<Evaluation, TwoPFluidSystem> TwoPhaseFluidState;
    typedef Ewoms::ImmiscibleFluidState<Evaluation, ThreePFluidSystem> ThreePhaseFluidState;

    // test conformance to the capillary pressure APIs
    {
        typedef Ewoms::BrooksCorey<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }
    {
        typedef Ewoms::LinearMaterial<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();

        typedef Ewoms::EffToAbsLaw<MaterialLaw> TwoPAbsLaw;
        testGenericApi<TwoPAbsLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<TwoPAbsLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<TwoPAbsLaw, TwoPhaseFluidState>();

        typedef Ewoms::LinearMaterial<ThreePhaseTraits> ThreePMaterialLaw;
        testGenericApi<ThreePMaterialLaw, ThreePhaseFluidState>();
        testThreePhaseApi<ThreePMaterialLaw, ThreePhaseFluidState>();
        //testThreePhaseSatApi<ThreePMaterialLaw, ThreePhaseFluidState>();

        typedef Ewoms::EffToAbsLaw<ThreePMaterialLaw> ThreePAbsLaw;
        testGenericApi<ThreePAbsLaw, ThreePhaseFluidState>();
        testThreePhaseApi<ThreePAbsLaw, ThreePhaseFluidState>();
        //testThreePhaseSatApi<ThreePAbsLaw, ThreePhaseFluidState>();
    }
    {
        typedef Ewoms::BrooksCorey<TwoPhaseTraits> TwoPhaseMaterial;
        typedef Ewoms::EclDefaultMaterial<ThreePhaseTraits,
                                        /*GasOilMaterial=*/TwoPhaseMaterial,
                                        /*OilWaterMaterial=*/TwoPhaseMaterial> MaterialLaw;
        testGenericApi<MaterialLaw, ThreePhaseFluidState>();
        testThreePhaseApi<MaterialLaw, ThreePhaseFluidState>();
        //testThreePhaseSatApi<MaterialLaw, ThreePhaseFluidState>();
    }
    {
        typedef Ewoms::BrooksCorey<TwoPhaseTraits> TwoPhaseMaterial;
        typedef Ewoms::EclStone1Material<ThreePhaseTraits,
                                       /*GasOilMaterial=*/TwoPhaseMaterial,
                                       /*OilWaterMaterial=*/TwoPhaseMaterial> MaterialLaw;
        testGenericApi<MaterialLaw, ThreePhaseFluidState>();
        testThreePhaseApi<MaterialLaw, ThreePhaseFluidState>();
        //testThreePhaseSatApi<MaterialLaw, ThreePhaseFluidState>();
    }
    {
        typedef Ewoms::BrooksCorey<TwoPhaseTraits> TwoPhaseMaterial;
        typedef Ewoms::EclStone2Material<ThreePhaseTraits,
                                       /*GasOilMaterial=*/TwoPhaseMaterial,
                                       /*OilWaterMaterial=*/TwoPhaseMaterial> MaterialLaw;
        testGenericApi<MaterialLaw, ThreePhaseFluidState>();
        testThreePhaseApi<MaterialLaw, ThreePhaseFluidState>();
        //testThreePhaseSatApi<MaterialLaw, ThreePhaseFluidState>();
    }
    {
        typedef Ewoms::BrooksCorey<TwoPhaseTraits> TwoPhaseMaterial;
        typedef Ewoms::EclTwoPhaseMaterial<ThreePhaseTraits,
                                         /*GasOilMaterial=*/TwoPhaseMaterial,
                                         /*OilWaterMaterial=*/TwoPhaseMaterial> MaterialLaw;
        testGenericApi<MaterialLaw, ThreePhaseFluidState>();
        testThreePhaseApi<MaterialLaw, ThreePhaseFluidState>();
        //testThreePhaseSatApi<MaterialLaw, ThreePhaseFluidState>();
    }
    {
        typedef Ewoms::BrooksCorey<TwoPhaseTraits> TwoPhaseMaterial;
        typedef Ewoms::EclMultiplexerMaterial<ThreePhaseTraits,
                                            /*GasOilMaterial=*/TwoPhaseMaterial,
                                            /*OilWaterMaterial=*/TwoPhaseMaterial> MaterialLaw;
        testGenericApi<MaterialLaw, ThreePhaseFluidState>();
        testThreePhaseApi<MaterialLaw, ThreePhaseFluidState>();
        //testThreePhaseSatApi<MaterialLaw, ThreePhaseFluidState>();
    }
    {
        typedef Ewoms::ThreePhaseParkerVanGenuchten<ThreePhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, ThreePhaseFluidState>();
        testThreePhaseApi<MaterialLaw, ThreePhaseFluidState>();
        //testThreePhaseSatApi<MaterialLaw, ThreePhaseFluidState>();
    }
    {
        typedef Ewoms::NullMaterial<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }
    {
        typedef Ewoms::NullMaterial<ThreePhaseTraits> ThreePMaterialLaw;
        testGenericApi<ThreePMaterialLaw, ThreePhaseFluidState>();
        testThreePhaseApi<ThreePMaterialLaw, ThreePhaseFluidState>();
        //testThreePhaseSatApi<ThreePMaterialLaw, ThreePhaseFluidState>();
    }
    {
        typedef Ewoms::ParkerLenhard<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }
    {
        typedef Ewoms::PiecewiseLinearTwoPhaseMaterial<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }
    {
        typedef Ewoms::SplineTwoPhaseMaterial<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }
    {
        typedef Ewoms::VanGenuchten<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }
    {
        typedef Ewoms::RegularizedBrooksCorey<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }
    {
        typedef Ewoms::RegularizedVanGenuchten<TwoPhaseTraits> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }

    {
        typedef Ewoms::BrooksCorey<TwoPhaseTraits> RawMaterialLaw;
        typedef Ewoms::EclEpsTwoPhaseLaw<RawMaterialLaw> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }
    {
        typedef Ewoms::BrooksCorey<TwoPhaseTraits> RawMaterialLaw;
        typedef Ewoms::EclHysteresisTwoPhaseLaw<RawMaterialLaw> MaterialLaw;
        testGenericApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseApi<MaterialLaw, TwoPhaseFluidState>();
        testTwoPhaseSatApi<MaterialLaw, TwoPhaseFluidState>();
    }
}

int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);

    testAll<double>();
    testAll<float>();

    return 0;
}
