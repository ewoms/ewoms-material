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
 * \copydoc Ewoms::SimpleHuDuanH2O
 */
#ifndef EWOMS_SIMPLE_HU_DUAN_H2O_HH
#define EWOMS_SIMPLE_HU_DUAN_H2O_HH

#include "component.hh"
#include "iapws/common.hh"

#include <ewoms/material/idealgas.hh>

#include <ewoms/common/mathtoolbox.hh>

#include <ewoms/common/unused.hh>

#include <cmath>

namespace Ewoms {

 /*!
 * \ingroup Components
 *
 * \brief A simple version of  pure water with density from Hu et al.
 *
 * Compared to the water formulation of IAPWS'97, this class provides
 * a much simpler component that represents the thermodynamic
 * properties of pure water. This implies that the likelyhood for
 * bugs in this class is reduced and the numerical performance is
 * increased. (At the cost of accuracy for the representation of the
 * physical quantities, of course.)
 *
 * Density from Hu, Duan, Zhu and Chou: PVTx properties of the CO2-H2O and CO2-H2O-NaCl
 * systems below 647 K: Assessment of experimental data and
 * thermodynamics models, Chemical Geology, 2007.
 *
 * \tparam Scalar The type used for representing scalar values
 */
template <class Scalar>
class SimpleHuDuanH2O : public Component<Scalar, SimpleHuDuanH2O<Scalar> >
{
    typedef Ewoms::IdealGas<Scalar> IdealGas;
    typedef IAPWS::Common<Scalar> Common;

    static const Scalar R;  // specific gas constant of water

public:
    /*!
     * \brief A human readable name for the water.
     */
    static const char* name()
    { return "H2O"; }

    /*!
     * \brief Returns true iff the gas phase is assumed to be compressible
     */
    static bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true iff the liquid phase is assumed to be compressible
     */
    static bool liquidIsCompressible()
    { return false; }

    /*!
     * \brief Returns true iff the gas phase is assumed to be ideal
     */
    static bool gasIsIdeal()
    { return true; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of water.
     */
    static Scalar molarMass()
    { return 18e-3; }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of water.
     */
    static Scalar criticalTemperature()
    { return 647.096; /* [K] */ }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of water.
     */
    static Scalar criticalPressure()
    { return 22.064e6; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at water's triple point.
     */
    static Scalar tripleTemperature()
    { return 273.16; /* [K] */ }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at water's triple point.
     */
    static Scalar triplePressure()
    { return 611.657; /* [N/m^2] */ }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure water
     *        at a given temperature.
     *
     *\param T temperature of component in \f$\mathrm{[K]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf
     */
    template <class Evaluation>
    static Evaluation vaporPressure(const Evaluation& T)
    {
        if (T > criticalTemperature())
            return criticalPressure();
        if (T < tripleTemperature())
            return 0; // water is solid: We don't take sublimation into account

        static const Scalar n[10] = {
            0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
            0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
            -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
            0.65017534844798e3
        };

        Evaluation sigma = T + n[8]/(T - n[9]);

        Evaluation A = (sigma + n[0])*sigma + n[1];
        Evaluation B = (n[2]*sigma + n[3])*sigma + n[4];
        Evaluation C = (n[5]*sigma + n[6])*sigma + n[7];

        Evaluation tmp = 2.0*C/(Ewoms::sqrt(B*B - 4.0*A*C) - B);
        tmp *= tmp;
        tmp *= tmp;

        return 1e6*tmp;
    }

    /*!
     * \brief Specific enthalpy of water steam \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasEnthalpy(const Evaluation& temperature,
                                  const Evaluation& /*pressure*/)
    { return 1.976e3*temperature + 40.65e3/molarMass(); }

    /*!
     * \copydoc Component::gasHeatCapacity
     */
    template <class Evaluation>
    static Evaluation gasHeatCapacity(const Evaluation& temperature EWOMS_UNUSED,
                                      const Evaluation& pressure EWOMS_UNUSED)
    { return 1.976e3; }

    /*!
     * \brief Specific enthalpy of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidEnthalpy(const Evaluation& temperature,
                                     const Evaluation& /*pressure*/)
    { return 4180*temperature; }

    /*!
     * \copydoc Component::liquidHeatCapacity
     */
    template <class Evaluation>
    static Evaluation liquidHeatCapacity(const Evaluation& temperature EWOMS_UNUSED,
                                         const Evaluation& pressure EWOMS_UNUSED)
    { return 4.184e3; }

    /*!
     * \brief Specific internal energy of steam \f$\mathrm{[J/kg]}\f$.
     *
     *        Definition of enthalpy: \f$h= u + pv = u + p / \rho\f$.
     *
     *        Rearranging for internal energy yields: \f$u = h - pv\f$.
     *
     *        Exploiting the Ideal Gas assumption (\f$pv = R_{\textnormal{specific}} T\f$)gives: \f$u = h - R / M T \f$.
     *
     *        The universal gas constant can only be used in the case of molar formulations.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasInternalEnergy(const Evaluation& temperature,
                                        const Evaluation& pressure)
    {
        return
            gasEnthalpy(temperature, pressure) -
            1/molarMass()* // conversion from [J/(mol K)] to [J/(kg K)]
            IdealGas::R*temperature; // = pressure *spec. volume for an ideal gas
    }

    /*!
     * \brief Specific internal energy of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidInternalEnergy(const Evaluation& temperature,
                                           const Evaluation& pressure)
    {
        return
            liquidEnthalpy(temperature, pressure) -
            pressure/liquidDensity(temperature, pressure);
    }

    /*!
     * \brief Specific heat conductivity of liquid water \f$\mathrm{[W/(m K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidThermalConductivity(const Evaluation& /*temperature*/,
                                                const Evaluation& /*pressure*/)
    {
        return 0.578078; // conductivity of liquid water [W / (m K ) ] IAPWS evaluated at p=.1 MPa, T=8°C
    }

    /*!
     * \brief Specific heat conductivity of steam \f$\mathrm{[W/(m K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasThermalConductivity(const Evaluation& /*temperature*/,
                                             const Evaluation& /*pressure*/)
    {
        return 0.028224; // conductivity of steam [W / (m K ) ] IAPWS evaluated at p=.1 MPa, T=8°C
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of steam at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation gasDensity(const Evaluation& temperature, const Evaluation& pressure)
    {
        // Assume an ideal gas
        return molarMass()*IdealGas::molarDensity(temperature, pressure);
    }

    /*!
     * \brief The pressure of steam in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    template <class Evaluation>
    static Evaluation gasPressure(const Evaluation& temperature, const Evaluation& density)
    {
        // Assume an ideal gas
        return IdealGas::pressure(temperature, density/molarMass());
    }

    /*!
     * \brief The density of pure water at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidDensity(const Evaluation& temperature, const Evaluation& pressure)
    {
        return liquidDensity_(temperature, pressure);
    }

    /*!
     * \brief The pressure of water in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidPressure(const Evaluation& /*temperature*/, const Evaluation& /*density*/)
    {
        throw std::logic_error("The liquid pressure is undefined for incompressible fluids");
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of steam.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     * \param regularize defines, if the functions is regularized or not, set to true by default
     */
    template <class Evaluation>
    static Evaluation gasViscosity(const Evaluation& /*temperature*/,
                                   const Evaluation& /*pressure*/)
    {
        return 1e-05;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure water.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidViscosity(const Evaluation& temperature, const Evaluation& pressure)
    {
        if (temperature > 570) {
            std::ostringstream oss;
            oss << "Viscosity of water based on Hu et al is too different from IAPWS for T above 570K and "
                << "(T = " << temperature << ")";
            throw NumericalIssue(oss.str());
        }

        const Evaluation& rho = liquidDensity(temperature, pressure);
        return Common::viscosity(temperature, rho);
    }

private:

    /*!
     * \brief The density of pure water in \f$\mathrm{[kg/m3]}\f$
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template <class Evaluation>
    static Evaluation liquidDensity_(const Evaluation& T, const Evaluation& pressure) {
        // Hu, Duan, Zhu and Chou: PVTx properties of the CO2-H2O and CO2-H2O-NaCl
        // systems below 647 K: Assessment of experimental data and
        // thermodynamics models, Chemical Geology, 2007.
        if (T > 647 || pressure > 100e6) {
            std::ostringstream oss;
            oss << "Density of water is only implemented for temperatures below 647K and "
                << "pressures below 100MPa. (T = " << T << ", p=" << pressure;
            throw NumericalIssue(oss.str());
        }

        Evaluation p = pressure / 1e6; // to MPa
        Scalar Mw = molarMass() * 1e3; //kg/kmol

        static const Scalar k0[5] = { 3.27225e-07, -4.20950e-04, 2.32594e-01, -4.16920e+01, 5.71292e+03 };
        static const Scalar k1[5] = { -2.32306e-10, 2.91138e-07, -1.49662e-04, 3.59860e-02, -3.55071 };
        static const Scalar k2[3] = { 2.57241e-14, -1.24336e-11, 5.42707e-07 };
        static const Scalar k3[3] = { -4.42028e-18, 2.10007e-15, -8.11491e-11 };
        Evaluation k0_eval = 1e-3 * (((k0[0]*T + k0[1])*T + k0[2])*T + k0[3] + k0[4]/T);
        Evaluation k1_eval = 1e-2 * (((k1[0]*T + k1[1])*T + k1[2])*T + k1[3] + k1[4]/T);
        Evaluation k2_eval = 1e-1 * ((k2[0]*T + k2[1])*T*T + k2[2]);
        Evaluation k3_eval = (k3[0]*T + k3[1])*T*T + k3[2];

        // molar volum (m³/kmol):
        Evaluation vw = ((k3_eval*p + k2_eval)*p + k1_eval)*p + k0_eval;

        // density kg/m3
        return Mw / vw;

    }

};

template <class Scalar>
const Scalar SimpleHuDuanH2O<Scalar>::R = Ewoms::Constants<Scalar>::R / 18e-3;

} // namespace Ewoms

#endif
