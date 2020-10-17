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
 * \copydoc Ewoms::VanGenuchten
 */
#ifndef EWOMS_VAN_GENUCHTEN_HH
#define EWOMS_VAN_GENUCHTEN_HH

#include "vangenuchtenparams.hh"

#include <ewoms/common/mathtoolbox.hh>

#include <algorithm>
#include <cmath>
#include <cassert>

namespace Ewoms {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Implementation of the van Genuchten capillary pressure -
 *        saturation relation.
 *
 * This class only implements the "raw" van-Genuchten curves as static
 * members and doesn't concern itself converting absolute to effective
 * saturations and vice versa.
 *
 * The converion from and to effective saturations can be done using,
 * e.g. EffToAbsLaw.
 *
 * \see VanGenuchtenParams
 */
template <class TraitsT, class ParamsT = VanGenuchtenParams<TraitsT> >
class VanGenuchten : public TraitsT
{
public:
    //! The traits class for this material law
    typedef TraitsT Traits;

    //! The type of the parameter objects for this law
    typedef ParamsT Params;

    //! The type of the scalar values for this law
    typedef typename Traits::Scalar Scalar;

    //! The number of fluid phases
    static const int numPhases = Traits::numPhases;
    static_assert(numPhases == 2,
                  "The van Genuchten capillary pressure law only "
                  "applies to the case of two fluid phases");

    //! Specify whether this material law implements the two-phase
    //! convenience API
    static const bool implementsTwoPhaseApi = true;

    //! Specify whether this material law implements the two-phase
    //! convenience API which only depends on the phase saturations
    static const bool implementsTwoPhaseSatApi = true;

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
     * \brief The capillary pressure-saturation curves according to van Genuchten.
     *
     * Van Genuchten's empirical capillary pressure <-> saturation
     * function is given by
     * \f[
     * p_{c,wn} = p_n - p_w = ({S_w}^{-1/m} - 1)^{1/n}/\alpha
     * \f]
     *
     * \param values A random access container which stores the
     *               relative pressure of each fluid phase.
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state for which the capillary pressure
     *           ought to be calculated
     */
    template <class Container, class FluidState>
    static void capillaryPressures(Container& values, const Params& params, const FluidState& fs)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = 0.0; // reference phase
        values[Traits::nonWettingPhaseIdx] = pcnw<FluidState, Evaluation>(params, fs);
    }

    /*!
     * \brief Calculate the saturations of the phases starting from
     *        their pressure differences.
     */
    template <class Container, class FluidState>
    static void saturations(Container& values, const Params& params, const FluidState& fs)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = Sw<FluidState, Evaluation>(params, fs);
        values[Traits::nonWettingPhaseIdx] = 1 - values[Traits::wettingPhaseIdx];
    }

    /*!
     * \brief The relative permeability-saturation curves according to van Genuchten.
     *
     * \param values A random access container which stores the
     *               relative permeability of each fluid phase.
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state for which the relative permeabilities
     *           ought to be calculated
     */
    template <class Container, class FluidState>
    static void relativePermeabilities(Container& values, const Params& params, const FluidState& fs)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = krw<FluidState, Evaluation>(params, fs);
        values[Traits::nonWettingPhaseIdx] = krn<FluidState, Evaluation>(params, fs);
    }

    /*!
     * \brief The capillary pressure-saturation curve according to van Genuchten.
     *
     * Van Genuchten's empirical capillary pressure <-> saturation
     * function is given by
     * \f[
     * p_{c,wn} = p_n - p_w = ({S_w}^{-1/m} - 1)^{1/n}/\alpha
     * \f]
     *
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state for which the capillary pressure
     *           ought to be calculated
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation pcnw(const Params& params, const FluidState& fs)
    {
        const Evaluation& Sw =
            Ewoms::decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));

        assert(0 <= Sw && Sw <= 1);

        return twoPhaseSatPcnw(params, Sw);
    }

    /*!
     * \brief The saturation-capillary pressure curve according to van
     *        Genuchten using a material law specific API.
     *
     * The advantage of this model is that it is simpler to use
     * because the baggage of the fluid state API does not need to be
     * carried along. The disavantage of this is, that it is very
     * specific to the van Genuchten law (i.e., depends only on the
     * wetting phase saturation, assumes two fluid phases, etc)
     *
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param Sw The effective wetting phase saturation
     */
    template <class Evaluation>
    static Evaluation twoPhaseSatPcnw(const Params& params, const Evaluation& Sw)
    {
        return Ewoms::pow(Ewoms::pow(Sw, -1.0/params.vgM()) - 1, 1.0/params.vgN())/params.vgAlpha();
    }

    /*!
     * \brief The saturation-capillary pressure curve according to van Genuchten.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
     * S_w = {p_C}^{-1} = ((\alpha p_C)^n + 1)^{-m}
     * \f]
     *
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state containing valid phase pressures
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sw(const Params& params, const FluidState& fs)
    {
        Evaluation pC =
            Ewoms::decay<Evaluation>(fs.pressure(Traits::nonWettingPhaseIdx))
            - Ewoms::decay<Evaluation>(fs.pressure(Traits::wettingPhaseIdx));
        return twoPhaseSatSw(params, pC);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatSw(const Params& params, const Evaluation& pC)
    {
        assert(pC >= 0);

        return Ewoms::pow(Ewoms::pow(params.vgAlpha()*pC, params.vgN()) + 1, -params.vgM());
    }

    /*!
     * \brief Calculate the non-wetting phase saturations depending on
     *        the phase pressures.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sn(const Params& params, const FluidState& fs)
    { return 1 - Sw<FluidState, Evaluation>(params, fs); }

    template <class Evaluation>
    static Evaluation twoPhaseSatSn(const Params& params, const Evaluation& pC)
    { return 1 - twoPhaseSatSw(params, pC); }

    /*!
     * \brief The relative permeability for the wetting phase of the
     *        medium according to van Genuchten's curve with Mualem
     *        parameterization.
     *
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state for which the relative permeability
     *           ought to be calculated
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krw(const Params& params, const FluidState& fs)
    {
        const Evaluation& Sw =
            Ewoms::decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));

        return twoPhaseSatKrw(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrw(const Params& params, const Evaluation& Sw)
    {
        assert(0.0 <= Sw && Sw <= 1.0);

        Evaluation r = 1.0 - Ewoms::pow(1.0 - Ewoms::pow(Sw, 1/params.vgM()), params.vgM());
        return Ewoms::sqrt(Sw)*r*r;
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium according to van Genuchten.
     *
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state for which the derivative
     *           ought to be calculated
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krn(const Params& params, const FluidState& fs)
    {
        const Evaluation& Sw =
            1.0 - Ewoms::decay<Evaluation>(fs.saturation(Traits::nonWettingPhaseIdx));

        return twoPhaseSatKrn(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrn(const Params& params, Evaluation Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        return
            Ewoms::pow(1 - Sw, 1.0/3) *
            Ewoms::pow(1 - Ewoms::pow(Sw, 1/params.vgM()), 2*params.vgM());
    }
};
} // namespace Ewoms

#endif
