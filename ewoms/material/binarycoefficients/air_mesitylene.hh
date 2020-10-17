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
 * \copydoc Ewoms::BinaryCoeff::Air_Mesitylene
 */
#ifndef EWOMS_BINARY_COEFF_AIR_MESITYLENE_HH
#define EWOMS_BINARY_COEFF_AIR_MESITYLENE_HH

#include <ewoms/material/components/air.hh>
#include <ewoms/material/components/mesitylene.hh>

namespace Ewoms
{
namespace BinaryCoeff
{

/*!
 * \brief Binary coefficients for water and mesitylene.
 */
class Air_Mesitylene
{
public:
    /*!
     *
     */
    template <class Evaluation>
    static Evaluation henry(const Evaluation& /*temperature*/)
    { throw std::runtime_error("Not implemented: Henry coefficient of air in mesitylene"); }

    /*!
     * \brief Binary diffusion coefficent [m^2/s] for air and mesitylene.
     * I used the method according to Wilke and Lee
     * see Handbook of chem. property's Estimation Methods
     * W.J. Lyman, W.F. Reehl, D.H. Rosenblatt
     *
     */
    template <class Evaluation>
    static Evaluation gasDiffCoeff(Evaluation temperature, Evaluation pressure)
    {
        typedef Ewoms::Air<double> Air;
        typedef Ewoms::Mesitylene<double> Mesitylene;

        temperature = Ewoms::max(temperature, 1e-9); // regularization
        temperature = Ewoms::min(temperature, 500.0); // regularization
        pressure = Ewoms::max(pressure, 0.0); // regularization
        pressure = Ewoms::min(pressure, 1e8); // regularization

        const double M_m = 1e3*Mesitylene::molarMass(); // [g/mol] molecular weight of mesitylene
        const double M_a = 1e3*Air::molarMass(); // [g/mol] molecular weight of air
        const double Tb_m = 437.9;        // [K] boiling temperature of mesitylene
        const double sigma_a = 3.711;     // charact. length of air
        const double T_scal_a = 78.6;     // [K] (molec. energy of attraction/Boltzmann constant)
        const double V_B_m = 162.6;       // [cm^3/mol] LeBas molal volume of mesitylene
        const double sigma_m = 1.18*std::pow(V_B_m, 0.333);     // charact. length of mesitylene
        const double sigma_am = 0.5*(sigma_a + sigma_m);
        const double T_scal_m = 1.15*Tb_m;
        const double T_scal_am = std::sqrt(T_scal_a*T_scal_m);

        Evaluation T_star = temperature/T_scal_am;
        T_star = Ewoms::max(T_star, 1e-5); // regularization

        const Evaluation Omega = 1.06036/Ewoms::pow(T_star, 0.1561) + 0.193/Ewoms::exp(T_star*0.47635)
            + 1.03587/Ewoms::exp(T_star*1.52996) + 1.76474/Ewoms::exp(T_star*3.89411);
        const double B_ = 0.00217 - 0.0005*std::sqrt(1.0/M_a + 1.0/M_m);
        const double Mr = (M_a + M_m)/(M_a*M_m);
        const Evaluation D_am = (B_*Ewoms::pow(temperature, 1.5) * std::sqrt(Mr))
                           /(1e-5*pressure*std::pow(sigma_am, 2.0) * Omega); // [cm^2/s]

        return 1e-4*D_am; // [m^2/s]
    }

    /*!
     * \brief Diffusion coefficent [m^2/s] for molecular mesitylene in liquid water.
     */
    template <class Evaluation>
    static Evaluation liquidDiffCoeff(const Evaluation& /*temperature*/, const Evaluation& /*pressure*/)
    { throw std::runtime_error("Not implemented: Binary liquid diffusion coefficients of air and mesitylene"); }
};

} // namespace BinaryCoeff
} // namespace Ewoms

#endif
