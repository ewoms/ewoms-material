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
 * \copydoc Ewoms::BinaryCoeff::fullerMethod
 */
#ifndef EWOMS_FULLERMETHOD_HH
#define EWOMS_FULLERMETHOD_HH

#include <ewoms/common/means.hh>
#include <ewoms/common/mathtoolbox.hh>

#include <cmath>

namespace Ewoms {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Estimate binary diffusion coefficents \f$\mathrm{[m^2/s]}\f$ in gases according to
 *        the method by Fuller.
 *
 * \param M molar masses \f$\mathrm{[g/mol]}\f$
 * \param SigmaNu atomic diffusion volume
 * \param temperature The temperature \f$\mathrm{[K]}\f$
 * \param pressure phase pressure \f$\mathrm{[Pa]}\f$
 *
 * This function estimates the diffusion coefficents in binary gases
 * using to the method proposed by Fuller. This method and is only
 * valid at "low" pressures.
 *
 * See: R. Reid, et al.: The Properties of Gases and Liquids, 4th
 * edition, McGraw-Hill, 1987, pp. 587-588
 */
template <class Scalar, class Evaluation = Scalar>
inline Evaluation fullerMethod(const Scalar* M, // molar masses [g/mol]
                               const Scalar* SigmaNu, // atomic diffusion volume
                               const Evaluation& temperature, // [K]
                               const Evaluation& pressure) // [Pa]
{
    // "effective" molar mass in [g/m^3]
    Scalar Mab = Ewoms::harmonicMean(M[0], M[1]);

    // Fuller's method
    const Evaluation& tmp = std::pow(SigmaNu[0], 1./3) + std::pow(SigmaNu[1], 1./3);
    return 1e-4 * (143.0*Ewoms::pow(temperature, 1.75))/(pressure*std::sqrt(Mab)*tmp*tmp);
}

} // namespace BinaryCoeff
} // namespace Ewoms

#endif
