// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::EclStone1Material
 */
#ifndef OPM_ECL_STONE1_MATERIAL_HPP
#define OPM_ECL_STONE1_MATERIAL_HPP

#include "EclStone1MaterialParams.hpp"

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/MathToolbox.hpp>

#include <opm/common/Exceptions.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <algorithm>
#include <cmath>

namespace Opm {

/*!
 * \ingroup material
 *
 * \brief Implements the second phase capillary pressure/relperm law suggested by Stone
 *        as used by the ECLipse simulator.
 *
 * This material law is valid for three fluid phases and only depends
 * on the saturations.
 *
 * The required two-phase relations are supplied by means of template
 * arguments and can be an arbitrary other material laws. (Provided
 * that these only depend on saturation.)
 */
template <class TraitsT,
          class GasOilMaterialLawT,
          class OilWaterMaterialLawT,
          class ParamsT = EclStone1MaterialParams<TraitsT, GasOilMaterialLawT, OilWaterMaterialLawT> >
class EclStone1Material : public TraitsT
{
public:
    typedef GasOilMaterialLawT GasOilMaterialLaw;
    typedef OilWaterMaterialLawT OilWaterMaterialLaw;

    // some safety checks
    static_assert(TraitsT::numPhases == 3,
                  "The number of phases considered by this capillary pressure "
                  "law is always three!");
    static_assert(GasOilMaterialLaw::numPhases == 2,
                  "The number of phases considered by the gas-oil capillary "
                  "pressure law must be two!");
    static_assert(OilWaterMaterialLaw::numPhases == 2,
                  "The number of phases considered by the oil-water capillary "
                  "pressure law must be two!");
    static_assert(std::is_same<typename GasOilMaterialLaw::Scalar,
                               typename OilWaterMaterialLaw::Scalar>::value,
                  "The two two-phase capillary pressure laws must use the same "
                  "type of floating point values.");

    static_assert(GasOilMaterialLaw::implementsTwoPhaseSatApi,
                  "The gas-oil material law must implement the two-phase saturation "
                  "only API to for the default Ecl capillary pressure law!");
    static_assert(OilWaterMaterialLaw::implementsTwoPhaseSatApi,
                  "The oil-water material law must implement the two-phase saturation "
                  "only API to for the default Ecl capillary pressure law!");

    typedef TraitsT Traits;
    typedef ParamsT Params;
    typedef typename Traits::Scalar Scalar;

    static const int numPhases = 3;
    static const int waterPhaseIdx = Traits::wettingPhaseIdx;
    static const int oilPhaseIdx = Traits::nonWettingPhaseIdx;
    static const int gasPhaseIdx = Traits::gasPhaseIdx;

    //! Specify whether this material law implements the two-phase
    //! convenience API
    static const bool implementsTwoPhaseApi = false;

    //! Specify whether this material law implements the two-phase
    //! convenience API which only depends on the phase saturations
    static const bool implementsTwoPhaseSatApi = false;

    //! Specify whether the quantities defined by this material law
    //! are saturation dependent
    static const bool isSaturationDependent = true;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the absolute pressure
    static const bool isPressureDependent = false;

    //! Specify whether the quantities defined by this material law
    //! are temperature dependent
    static const bool isTemperatureDependent = false;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the phase composition
    static const bool isCompositionDependent = false;

    /*!
     * \brief Implements the default three phase capillary pressure law
     *        used by the ECLipse simulator.
     *
     * This material law is valid for three fluid phases and only
     * depends on the saturations.
     *
     * The required two-phase relations are supplied by means of template
     * arguments and can be an arbitrary other material laws.
     *
     * \param values Container for the return values
     * \param params Parameters
     * \param state The fluid state
     */
    template <class ContainerT, class FluidState>
    static void capillaryPressures(ContainerT &values,
                                   const Params &params,
                                   const FluidState &state)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;
        values[gasPhaseIdx] = pcgn<FluidState, Evaluation>(params, state);
        values[oilPhaseIdx] = 0;
        values[waterPhaseIdx] = - pcnw<FluidState, Evaluation>(params, state);
        Valgrind::CheckDefined(values[gasPhaseIdx]);
        Valgrind::CheckDefined(values[oilPhaseIdx]);
        Valgrind::CheckDefined(values[waterPhaseIdx]);
    }

    /*!
     * \brief Capillary pressure between the gas and the non-wetting
     *        liquid (i.e., oil) phase.
     *
     * This is defined as
     * \f[
     * p_{c,gn} = p_g - p_n
     * \f]
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation pcgn(const Params &params,
                           const FluidState &fs)
    {
        typedef MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& Sw = 1.0 - FsToolbox::template toLhs<Evaluation>(fs.saturation(gasPhaseIdx));
        return GasOilMaterialLaw::twoPhaseSatPcnw(params.gasOilParams(), Sw);
    }

    /*!
     * \brief Capillary pressure between the non-wetting liquid (i.e.,
     *        oil) and the wetting liquid (i.e., water) phase.
     *
     * This is defined as
     * \f[
     * p_{c,nw} = p_n - p_w
     * \f]
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation pcnw(const Params &params,
                           const FluidState &fs)
    {
        typedef MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& Sw = FsToolbox::template toLhs<Evaluation>(fs.saturation(waterPhaseIdx));
        Valgrind::CheckDefined(Sw);
        const auto& result = OilWaterMaterialLaw::twoPhaseSatPcnw(params.oilWaterParams(), Sw);
        Valgrind::CheckDefined(result);
        return result;
    }

    /*!
     * \brief The inverse of the capillary pressure
     */
    template <class ContainerT, class FluidState>
    static void saturations(ContainerT& /* values */,
                            const Params& /* params */,
                            const FluidState& /* fluidState */)
    {
        OPM_THROW(std::logic_error, "Not implemented: saturations()");
    }

    /*!
     * \brief The saturation of the gas phase.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sg(const Params& /* params */,
                         const FluidState& /* fluidState */)
    {
        OPM_THROW(std::logic_error, "Not implemented: Sg()");
    }

    /*!
     * \brief The saturation of the non-wetting (i.e., oil) phase.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sn(const Params& /* params */,
                         const FluidState& /* fluidState */)
    {
        OPM_THROW(std::logic_error, "Not implemented: Sn()");
    }

    /*!
     * \brief The saturation of the wetting (i.e., water) phase.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sw(const Params& /* params */,
                         const FluidState& /* fluidState */)
    {
        OPM_THROW(std::logic_error, "Not implemented: Sw()");
    }

    /*!
     * \brief The relative permeability of all phases.
     *
     * The relative permeability of the water phase it uses the same
     * value as the relative permeability for water in the water-oil
     * law with \f$S_o = 1 - S_w\f$. The gas relative permebility is
     * taken from the gas-oil material law, but with \f$S_o = 1 -
     * S_g\f$.  The relative permeability of the oil phase is
     * calculated using the relative permeabilities of the oil phase
     * in the two two-phase systems.
     *
     * A more detailed description can be found in the "Three phase
     * oil relative permeability models" section of the ECLipse
     * technical description.
     */
    template <class ContainerT, class FluidState>
    static void relativePermeabilities(ContainerT &values,
                                       const Params &params,
                                       const FluidState &fluidState)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[waterPhaseIdx] = krw<FluidState, Evaluation>(params, fluidState);
        values[oilPhaseIdx] = krn<FluidState, Evaluation>(params, fluidState);
        values[gasPhaseIdx] = krg<FluidState, Evaluation>(params, fluidState);
    }

    /*!
     * \brief The relative permeability of the gas phase.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krg(const Params &params,
                          const FluidState &fluidState)
    {
        typedef MathToolbox<typename FluidState::Scalar> FsToolbox;

        const Evaluation& Sw = 1 - FsToolbox::template toLhs<Evaluation>(fluidState.saturation(gasPhaseIdx));
        return GasOilMaterialLaw::twoPhaseSatKrn(params.gasOilParams(), Sw);
    }

    /*!
     * \brief The relative permeability of the wetting phase.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krw(const Params &params,
                          const FluidState &fluidState)
    {
        typedef MathToolbox<typename FluidState::Scalar> FsToolbox;

        const Evaluation& Sw = FsToolbox::template toLhs<Evaluation>(fluidState.saturation(waterPhaseIdx));
        return OilWaterMaterialLaw::twoPhaseSatKrw(params.oilWaterParams(), Sw);
    }

    /*!
     * \brief The relative permeability of the non-wetting (i.e., oil) phase.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krn(const Params &params,
                          const FluidState &fluidState)
    {
        typedef MathToolbox<Evaluation> Toolbox;
        typedef MathToolbox<typename FluidState::Scalar> FsToolbox;

        // the Eclipse docu is inconsistent with naming the variable of connate water: In
        // some places the connate water saturation is represented by "Swl", in others
        // "Swco" is used.
        Scalar Swco = params.Swl();

        // oil relperm at connate water saturations (with Sg=0)
        Scalar krocw = params.krocw();

        const Evaluation& Sw = FsToolbox::template toLhs<Evaluation>(fluidState.saturation(waterPhaseIdx));
        const Evaluation& Sg = FsToolbox::template toLhs<Evaluation>(fluidState.saturation(gasPhaseIdx));

        Evaluation kro_ow = OilWaterMaterialLaw::twoPhaseSatKrn(params.oilWaterParams(), Sw);
        Evaluation kro_go = GasOilMaterialLaw::twoPhaseSatKrw(params.gasOilParams(), 1 - Sg - Swco);

        Evaluation beta;
        if (Sw <= Swco)
            beta = 1.0;
        else {
            // there seems to be an error in the ECL documentation: using the approach to
            // the scaled saturations as described there leads to significant deviations
            // from the results produced by Eclipse 100.
            Evaluation SSw = (Sw - Swco)/(1.0 - Swco);
            Evaluation SSg = Sg/(1.0 - Swco);
            Evaluation SSo = 1.0 - SSw - SSg;

            if (SSw >= 1.0 || SSg >= 1.0)
                beta = 1.0;
            else
                beta = Toolbox::pow( SSo/((1 - SSw)*(1 - SSg)), params.eta());
        }

        return Toolbox::max(0.0, Toolbox::min(1.0, beta*kro_ow*kro_go/krocw));
    }

    /*!
     * \brief Update the hysteresis parameters after a time step.
     *
     * This assumes that the nested two-phase material laws are parameters for
     * EclHysteresisLaw. If they are not, calling this methid will cause a compiler
     * error. (But not calling it will still work.)
     */
    template <class FluidState>
    static void updateHysteresis(Params &params, const FluidState &fluidState)
    {
        typedef MathToolbox<typename FluidState::Scalar> FsToolbox;

        Scalar Sw = FsToolbox::value(fluidState.saturation(waterPhaseIdx));
        Scalar Sg = FsToolbox::value(fluidState.saturation(gasPhaseIdx));

        params.oilWaterParams().update(/*pcSw=*/Sw, /*krwSw=*/Sw, /*krnSw=*/Sw);
        params.gasOilParams().update(/*pcSw=*/1 - Sg, /*krwSw=*/1 - Sg, /*krnSw=*/1 - Sg);
    }
};
} // namespace Opm

#endif