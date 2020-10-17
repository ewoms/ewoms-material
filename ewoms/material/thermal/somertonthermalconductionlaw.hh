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
 * \copydoc Ewoms::SomertonThermalConductionLaw
 */
#ifndef EWOMS_SOMERTON_THERMAL_CONDUCTION_LAW_HH
#define EWOMS_SOMERTON_THERMAL_CONDUCTION_LAW_HH

#include "somertonthermalconductionlawparams.hh"

#include <ewoms/common/spline.hh>

#include <ewoms/common/valgrind.hh>
#include <ewoms/common/mathtoolbox.hh>

#include <algorithm>

namespace Ewoms {

/*!
 * \ingroup material
 *
 * \brief Implements the Somerton law of thermal conductivity in a
 *        porous medium.
 *
 * See:
 *
 * W.H. Somerton, A.H. El-Shaarani and S.M. Mobarak: High
 * Temperature Behavior of Rocks Associated with Geothermal Type
 * Reservoirs, paper SPE-4897 presentet at SPE California Regional
 * Meeting 1974, 1974
 *
 * or
 *
 * H. Class: Theorie und numerische Modellierung nichtisothermer
 * Mehrphasenprozesse in NAPL kontaminierten poroesen Medien, PhD
 * thesis, Technical University of Braunschweig, 2000
 */
template <class FluidSystem,
          class ScalarT,
          class ParamsT = SomertonThermalConductionLawParams<FluidSystem::numPhases, ScalarT> >
class SomertonThermalConductionLaw
{
    enum { numPhases = FluidSystem::numPhases };

public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief Given a fluid state, return the effective thermal conductivity [W/m^2 / (K/m)] of the porous
     *        medium.
     *
     * For two phases, the Somerton law is given by:
     * \f[
     \lambda_{pm} =
     \lambda_{ful,g} +
     \sqrt{S_w}(\lambda_{ful,w} - \lambda_{vac}) +
     \sqrt{S_n}(\lambda_{ful,n} - \lambda_{vac})
     \f]
     *
     * where \f$\lambda_{vac}\f$ is the thermal conductivity of the
     * porous medium at vacuum, \f$\lambda_{ful,\alpha}\f$ is the thermal
     * conductivty of the porous medium if it is fully saturated by
     * phase \f$\alpha\f$ and \f$S_\alpha\f$ is the saturation of
     * phase \f$\alpha\f$.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation thermalConductivity(const Params& params,
                                       const FluidState& fluidState)
    {
        Valgrind::CheckDefined(params.vacuumLambda());

        Evaluation lambda = 0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Valgrind::CheckDefined(params.fullySaturatedLambda(phaseIdx));

            if (FluidSystem::isLiquid(phaseIdx)) {
                const auto& sat = Ewoms::decay<Evaluation>(fluidState.saturation(phaseIdx));
                lambda +=
                    regularizedSqrt_(Ewoms::max(0.0, Ewoms::min(1.0, sat)))
                    * (params.fullySaturatedLambda(phaseIdx) - params.vacuumLambda());
            }
            else { // gas phase
                lambda += params.fullySaturatedLambda(phaseIdx) - params.vacuumLambda();
            }
        };

        lambda += params.vacuumLambda();
        assert(lambda >= 0);
        return lambda;
    }

protected:
    template <class Evaluation>
    static Evaluation regularizedSqrt_(const Evaluation& x)
    {
        typedef Ewoms::Spline<Scalar> Spline;

        static const Scalar xMin = 1e-2;
        static const Scalar sqrtXMin = std::sqrt(xMin);
        static const Scalar fPrimeXMin = 1.0/(2*std::sqrt(xMin));
        static const Scalar fPrime0 = 2*fPrimeXMin;
        static const Spline sqrtRegSpline(0, xMin, // x0, x1
                                          0, sqrtXMin, // y0, y1
                                          fPrime0, fPrimeXMin); // m0, m1

        if (x > xMin)
            return Ewoms::sqrt(x);
        else if (x <= 0)
            return fPrime0 * x;
        else
            return sqrtRegSpline.eval(x);
    }
};
} // namespace Ewoms

#endif
