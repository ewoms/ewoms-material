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
 * \copydoc Ewoms::BinaryCoeff::H2O_Xylene
 */
#ifndef EWOMS_BINARY_COEFF_H2O_XYLENE_HH
#define EWOMS_BINARY_COEFF_H2O_XYLENE_HH

#include <ewoms/material/components/h2o.hh>
#include <ewoms/material/components/xylene.hh>

namespace Ewoms {
namespace BinaryCoeff {

/*!
 * \brief Binary coefficients for water and xylene.
 */
class H2O_Xylene
{
public:
    /*!
     * \brief Henry coefficent \f$[N/m^2]\f$  for xylene in liquid water.
     *
     * See:
     *
     *  Sanders1999 Henry collection
     */

    template <class Evaluation>
    static Evaluation henry(const Evaluation& /*temperature*/)
    {
        // after Sanders
        double sanderH = 1.5e-1;    //[M/atm]
        //conversion to our Henry definition
        double ewomsH = sanderH / 101.325; // has now [(mol/m^3)/Pa]
        ewomsH *= 18.02e-6;  //multiplied by molar volume of reference phase = water
        return 1.0/ewomsH; // [Pa]
    }

    /*!
     * \brief Binary diffusion coefficent [m^2/s] for molecular water and xylene.
     *
     */
    template <class Evaluation>
    static Evaluation gasDiffCoeff(Evaluation temperature, Evaluation pressure)
    {
        typedef Ewoms::H2O<double> H2O;
        typedef Ewoms::Xylene<double> Xylene;

        temperature = Ewoms::max(temperature, 1e-9); // regularization
        temperature = Ewoms::min(temperature, 500.0); // regularization
        pressure = Ewoms::max(pressure, 0.0); // regularization
        pressure = Ewoms::min(pressure, 1e8); // regularization

        const double M_x = 1e3*Xylene::molarMass(); // [g/mol] molecular weight of xylene
        const double M_w = 1e3*H2O::molarMass(); // [g/mol] molecular weight of water
        const double Tb_x = 412.9;        // [K] boiling temperature of xylene
        const double Tb_w = 373.15;       // [K] boiling temperature of water (at p_atm)
        const double V_B_w = 18.0;                // [cm^3/mol] LeBas molal volume of water
        const double sigma_w = 1.18*std::pow(V_B_w, 0.333);     // charact. length of air
        const double T_scal_w = 1.15*Tb_w;     // [K] (molec. energy of attraction/Boltzmann constant)
        const double V_B_x = 140.4;       // [cm^3/mol] LeBas molal volume of xylene
        const double sigma_x = 1.18*std::pow(V_B_x, 0.333);     // charact. length of xylene
        const double sigma_wx = 0.5*(sigma_w + sigma_x);
        const double T_scal_x = 1.15*Tb_x;
        const double T_scal_wx = std::sqrt(T_scal_w*T_scal_x);

        const Evaluation& T_star = Ewoms::max(temperature/T_scal_wx, 1e-5);

        const Evaluation& Omega = 1.06036/Ewoms::pow(T_star, 0.1561) + 0.193/Ewoms::exp(T_star*0.47635)
            + 1.03587/Ewoms::exp(T_star*1.52996) + 1.76474/Ewoms::exp(T_star*3.89411);
        const double  B_ = 0.00217 - 0.0005*std::sqrt(1.0/M_w + 1.0/M_x);
        const double Mr = (M_w + M_x)/(M_w*M_x);
        return 1e-4
            *(B_*Ewoms::pow(temperature,1.6)*std::sqrt(Mr))
            /(1e-5*pressure*std::pow(sigma_wx, 2.0)*Omega);
    }

    /*!
     * \brief Diffusion coefficent [m^2/s] for xylene in liquid water.
     *
     * \todo
     */
    template <class Evaluation>
    static Evaluation liquidDiffCoeff(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    {
        return 1.e-9;  // This is just an order of magnitude. Please improve it!
    }
};

} // namespace BinaryCoeff
} // namespace Ewoms

#endif
