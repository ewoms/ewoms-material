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
 * \copydoc Ewoms::CO2
 */
#ifndef EWOMS_CO2_HH
#define EWOMS_CO2_HH

#include <ewoms/material/constants.hh>
#include <ewoms/material/idealgas.hh>
#include <ewoms/material/components/component.hh>
#include <ewoms/common/mathtoolbox.hh>

#include <cmath>
#include <iostream>

namespace Ewoms {

/*!
 * \brief A class for the CO2 fluid properties
 *
 * Under reservoir conditions, CO2 is typically in supercritical
 * state. These properties can be provided in tabulated form, which is
 * necessary for this component. The template is used by the
 * fluidsystem \c BrineCO2FluidSystem. If thermodynamic precision
 * is not a top priority, the much simpler component \c Ewoms::SimpleCO2 can be
 * used instead
 */
template <class Scalar, class CO2Tables>
class CO2 : public Component<Scalar, CO2<Scalar, CO2Tables> >
{
    static const Scalar R;
    static bool warningPrinted;

public:
    /*!
     * \brief A human readable name for the CO2.
     */
    static const char* name()
    { return "CO2"; }

    /*!
     * \brief The mass in [kg] of one mole of CO2.
     */
    static Scalar molarMass()
    { return 44e-3; }

    /*!
     * \brief Returns the critical temperature [K] of CO2
     */
    static Scalar criticalTemperature()
    { return 273.15 + 30.95; /* [K] */ }

    /*!
     * \brief Returns the critical pressure [Pa] of CO2
     */
    static Scalar criticalPressure()
    { return 73.8e5; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature [K]at CO2's triple point.
     */
    static Scalar tripleTemperature()
    { return 273.15 - 56.35; /* [K] */ }

    /*!
     * \brief Returns the pressure [Pa] at CO2's triple point.
     */
    static Scalar triplePressure()
    { return 5.11e5; /* [N/m^2] */ }

    /*!
     * \brief Returns the pressure [Pa] at CO2's triple point.
     */
    static Scalar minTabulatedPressure()
    { return CO2Tables::tabulatedEnthalpy.minPress(); /* [N/m^2] */ }

    /*!
     * \brief Returns the pressure [Pa] at CO2's triple point.
     */
    static Scalar maxTabulatedPressure()
    { return CO2Tables::tabulatedEnthalpy.maxPress(); /* [N/m^2] */ }

    /*!
     * \brief Returns the pressure [Pa] at CO2's triple point.
     */
    static Scalar minTabulatedTemperature()
    { return CO2Tables::tabulatedEnthalpy.minTemp(); /* [N/m^2] */ }

    /*!
     * \brief Returns the pressure [Pa] at CO2's triple point.
     */
    static Scalar maxTabulatedTemperature()
    { return CO2Tables::tabulatedEnthalpy.maxTemp(); /* [N/m^2] */ }

    /*!
     * \brief The vapor pressure in [N/m^2] of pure CO2
     *        at a given temperature.
     *
     * See:
     *
     * R. Span and W. Wagner: A New Equation of State for Carbon
     * Dioxide Covering the Fluid Region from the Triple‐Point
     * Temperature to 1100 K at Pressures up to 800 MPa. Journal of
     * Physical and Chemical Reference Data, 25 (6), pp. 1509-1596,
     * 1996
     */
    template <class Evaluation>
    static Evaluation vaporPressure(const Evaluation& T)
    {
        static const Scalar a[4] =
            { -7.0602087, 1.9391218, -1.6463597, -3.2995634 };
        static const Scalar t[4] =
            { 1.0, 1.5, 2.0, 4.0 };

        // this is on page 1524 of the reference
        Evaluation exponent = 0;
        Evaluation Tred = T/criticalTemperature();
        for (int i = 0; i < 5; ++i)
            exponent += a[i]*Ewoms::pow(1 - Tred, t[i]);
        exponent *= 1.0/Tred;

        return Ewoms::exp(exponent)*criticalPressure();
    }

    /*!
     * \brief Returns true iff the gas phase is assumed to be compressible
     */
    static bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static bool gasIsIdeal()
    { return false; }

    /*!
     * \brief Specific enthalpy of gaseous CO2 [J/kg].
     */
    template <class Evaluation>
    static Evaluation gasEnthalpy(const Evaluation& temperature,
                                  const Evaluation& pressure)
    {
        return CO2Tables::tabulatedEnthalpy.eval(temperature, pressure);
    }

    /*!
     * \brief Specific internal energy of CO2 [J/kg].
     */
    template <class Evaluation>
    static Evaluation gasInternalEnergy(const Evaluation& temperature,
                                        const Evaluation& pressure)
    {
        const Evaluation& h = gasEnthalpy(temperature, pressure);
        const Evaluation& rho = gasDensity(temperature, pressure);

        return h - (pressure / rho);
    }

    /*!
     * \brief The density of CO2 at a given pressure and temperature [kg/m^3].
     */
    template <class Evaluation>
    static Evaluation gasDensity(const Evaluation& temperature, const Evaluation& pressure)
    {
        return CO2Tables::tabulatedDensity.eval(temperature, pressure);
    }

    /*!
     * \brief The dynamic viscosity [Pa s] of CO2.
     *
     * Equations given in: - Vesovic et al., 1990
     *                        - Fenhour etl al., 1998
     */
    template <class Evaluation>
    static Evaluation gasViscosity(Evaluation temperature, const Evaluation& pressure)
    {
        const Scalar a0 = 0.235156;
        const Scalar a1 = -0.491266;
        const Scalar a2 = 5.211155e-2;
        const Scalar a3 = 5.347906e-2;
        const Scalar a4 = -1.537102e-2;

        const Scalar d11 = 0.4071119e-2;
        const Scalar d21 = 0.7198037e-4;
        const Scalar d64 = 0.2411697e-16;
        const Scalar d81 = 0.2971072e-22;
        const Scalar d82 = -0.1627888e-22;

        const Scalar ESP = 251.196;

        if(temperature < 275.) // regularization
            temperature = 275.0;
        Evaluation TStar = temperature/ESP;

        // mu0: viscosity in zero-density limit
        const Evaluation& logTStar = Ewoms::log(TStar);
        Evaluation SigmaStar = Ewoms::exp(a0 + logTStar*(a1 + logTStar*(a2 + logTStar*(a3 + logTStar*a4))));

        Evaluation mu0 = 1.00697*Ewoms::sqrt(temperature) / SigmaStar;

        const Evaluation& rho = gasDensity(temperature, pressure); // CO2 mass density [kg/m^3]

        // dmu : excess viscosity at elevated density
        Evaluation dmu =
            d11*rho
            + d21*rho*rho
            + d64*Ewoms::pow(rho, 6.0)/(TStar*TStar*TStar)
            + d81*Ewoms::pow(rho, 8.0)
            + d82*Ewoms::pow(rho, 8.0)/TStar;

        return (mu0 + dmu)/1.0e6; // conversion to [Pa s]
    }

    /*!
     * \brief Specific isobaric heat capacity of the component [J/kg]
     *        as a liquid.
     *
     * This function uses the fact that heat capacity is the partial
     * derivative of enthalpy function with respect to temperature.
     *
     * \param temperature Temperature of component \f$\mathrm{[K]}\f$
     * \param pressure Pressure of component \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasHeatCapacity(const Evaluation& temperature, const Evaluation& pressure)
    {
        Scalar eps = 1e-6;

        // use central differences here because one-sided methods do
        // not come with a performance improvement. (central ones are
        // more accurate, though...)
        const Evaluation& h1 = gasEnthalpy(temperature - eps, pressure);
        const Evaluation& h2 = gasEnthalpy(temperature + eps, pressure);

        return (h2 - h1) / (2*eps) ;
    }
};

template <class Scalar, class CO2Tables>
bool CO2<Scalar, CO2Tables>::warningPrinted = false;

template <class Scalar, class CO2Tables>
const Scalar CO2<Scalar, CO2Tables>::R = Constants<Scalar>::R;

} // namespace Ewoms

#endif
