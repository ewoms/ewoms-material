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
 * \copydoc Ewoms::BinaryCoeff::H2O_CO2
 */
#ifndef EWOMS_BINARY_COEFF_H2O_CO2_HH
#define EWOMS_BINARY_COEFF_H2O_CO2_HH

#include <ewoms/material/binarycoefficients/henryiapws.hh>
#include <ewoms/material/binarycoefficients/fullermethod.hh>

#include <ewoms/material/components/h2o.hh>
#include <ewoms/material/components/simpleco2.hh>

namespace Ewoms {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and CO2.
 */
class H2O_CO2
{
public:
    /*!
     * \brief Henry coefficent \f$[N/m^2]\f$  for molecular CO2 in liquid water.
     *
     * See:
     *
     * IAPWS: "Guideline on the Henry's Constant and Vapor-Liquid
     * Distribution Constant for Gases in H2O and D2O at High
     * Temperatures"
     * http://www.iapws.org/relguide/HenGuide.pdf
     */
    template <class Scalar, class Evaluation = Scalar>
    static Evaluation henry(const Evaluation& temperature)
    {
        const Scalar E = 1672.9376;
        const Scalar F = 28.1751;
        const Scalar G = -112.4619;
        const Scalar H = 85.3807;

        return henryIAPWS(E, F, G, H, temperature);
    }

    /*!
     * \brief Binary diffusion coefficent [m^2/s] for molecular water and CO2.
     *
     * To calculate the values, the \ref fullerMethod is used.
     */
    template <class Scalar, class Evaluation = Scalar>
    static Evaluation gasDiffCoeff(const Evaluation& temperature, const Evaluation& pressure)
    {
        typedef Ewoms::H2O<Scalar> H2O;
        typedef Ewoms::SimpleCO2<Scalar> CO2;

        // atomic diffusion volumes
        const Scalar SigmaNu[2] = { 13.1 /* H2O */,  26.9 /* CO2 */ };
        // molar masses [g/mol]
        const Scalar M[2] = { H2O::molarMass()*1e3, CO2::molarMass()*1e3 };

        return fullerMethod(M, SigmaNu, temperature, pressure);
    }

    /*!
     * \brief Diffusion coefficent [m^2/s] for molecular CO2 in liquid water.
     */
    template <class Scalar, class Evaluation = Scalar>
    static Evaluation liquidDiffCoeff(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    { throw std::runtime_error("Not implemented: Binary liquid diffusion coefficients of CO2 and CH4"); }
};

} // namespace BinaryCoeff
} // namespace Ewoms

#endif
