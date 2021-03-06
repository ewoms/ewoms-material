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
 * \brief The IAPWS formulation of Henry coefficients in water.
 */
#ifndef EWOMS_HENRY_IAPWS_HH
#define EWOMS_HENRY_IAPWS_HH

#include <ewoms/material/components/h2o.hh>

namespace Ewoms
{
/*!
 * \ingroup Binarycoefficients
 * \brief The Henry constants in liquid water using the IAPWS 2004 formulation.
 *
 * This function calculates \f$K_D\f$, see:
 *
 * IAPWS: "Guideline on the Henry's Constant and Vapor-Liquid Distribution Constant for
 * Gases in H2O and D2O at High Temperatures" http://www.iapws.org/relguide/HenGuide.pdf
 */
template <class Scalar, class Evaluation>
inline Evaluation henryIAPWS(Scalar E,
                             Scalar F,
                             Scalar G,
                             Scalar H,
                             const Evaluation& temperature)
{
    typedef Ewoms::H2O<Evaluation> H2O;

    Evaluation Tr = temperature/H2O::criticalTemperature();
    Evaluation tau = 1 - Tr;

    static const Scalar c[6] = {
        1.99274064, 1.09965342, -0.510839303,
        -1.75493479,-45.5170352, -6.7469445e5
    };
    static const Scalar d[6] = {
        1/3.0, 2/3.0, 5/3.0,
        16/3.0, 43/3.0, 110/3.0
    };
    static const Scalar q = -0.023767;

    Evaluation f = 0;
    for (int i = 0; i < 6; ++i) {
        f += c[i]*Ewoms::pow(tau, d[i]);
    }

    const Evaluation& exponent =
        q*F +
        E/temperature*f +
        (F +
         G*Ewoms::pow(tau, 2.0/3) +
         H*tau)*
        Ewoms::exp((H2O::tripleTemperature() - temperature)/100);
    // CAUTION: K_D is formulated in mole fractions. We have to
    // multiply it with the vapor pressure of water in order to get
    // derivative of the partial pressure.
    return Ewoms::exp(exponent)*H2O::vaporPressure(temperature);
}
} // namespace Ewoms

#endif // EWOMS_HENRY_IAPWS_HH
