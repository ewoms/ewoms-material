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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Ewoms::Co2GasPvt
 */
#ifndef EWOMS_CO2_GAS_PVT_HH
#define EWOMS_CO2_GAS_PVT_HH

#include <ewoms/material/constants.hh>

#include <ewoms/common/tabulated1dfunction.hh>
#include <ewoms/material/components/co2.hh>
#include <ewoms/common/uniformtabulated2dfunction.hh>

#if HAVE_ECL_INPUT
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tablemanager.hh>
#endif

#include <vector>

namespace Ewoms {

namespace CO2GasPvtDefaultTables {
#include <ewoms/material/components/co2tables.inc>
}

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the gas phase
 * for CO2
 */
template <class Scalar, class CO2 = Ewoms::CO2<Scalar, Ewoms::CO2GasPvtDefaultTables::CO2Tables> >
class Co2GasPvt
{
    typedef std::vector<std::pair<Scalar, Scalar> > SamplingPoints;

public:
    typedef Ewoms::Tabulated1DFunction<Scalar> TabulatedOneDFunction;

    explicit Co2GasPvt() = default;
    Co2GasPvt(const std::vector<Scalar>& gasReferenceDensity)
        : gasReferenceDensity_(gasReferenceDensity)
    {}

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize the parameters for co2 gas using an ECL deck.
     */
    void initFromEclState(const EclipseState& eclState, const Schedule&)
    {
        if( !eclState.getTableManager().getDensityTable().empty()) {
            std::cerr << "WARNING: CO2STOR is enabled but DENSITY is in the deck. \n" <<
                         "The surface density is computed based on CO2-BRINE PVT at standard conditions (STCOND) and DENSITY is ignored " << std::endl;
        }

        if( eclState.getTableManager().hasTables("PVDG") || !eclState.getTableManager().getPvtgTables().empty()) {
            std::cerr << "WARNING: CO2STOR is enabled but PVDG or PVTG is in the deck. \n" <<
                         "CO2 PVT properties are computed based on the Span-Wagner pvt model and PVDG/PVTG input is ignored. " << std::endl;
        }

        // We only supported single pvt region for the co2-brine module
        size_t numRegions = 1;
        setNumRegions(numRegions);
        size_t regionIdx = 0;
        Scalar T_ref = eclState.getTableManager().stCond().temperature;
        Scalar P_ref = eclState.getTableManager().stCond().pressure;
        gasReferenceDensity_[regionIdx] = CO2::gasDensity(T_ref, P_ref);
        initEnd();
    }
#endif

    void setNumRegions(size_t numRegions)
    {
        gasReferenceDensity_.resize(numRegions);
    }

    /*!
     * \brief Initialize the reference densities of all fluids for a given PVT region
     */
    void setReferenceDensities(unsigned regionIdx,
                               Scalar /*rhoRefOil*/,
                               Scalar rhoRefGas,
                               Scalar /*rhoRefWater*/)
    {
        gasReferenceDensity_[regionIdx] = rhoRefGas;
    }

    /*!
     * \brief Finish initializing the oil phase PVT properties.
     */
    void initEnd()
    {

    }

    /*!
     * \brief Return the number of PVT regions which are considered by this PVT-object.
     */
    unsigned numRegions() const
    { return gasReferenceDensity_.size(); }

    /*!
     * \brief Returns the specific enthalpy [J/kg] of gas given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx EWOMS_UNUSED,
                        const Evaluation& temperature EWOMS_UNUSED,
                        const Evaluation& pressure EWOMS_UNUSED,
                        const Evaluation& Rv EWOMS_UNUSED) const
    {
        throw std::runtime_error("Requested the enthalpy of gas but the thermal option is not enabled");
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& /*Rv*/) const
    { return saturatedViscosity(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of oil saturated gas at given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned /*regionIdx*/,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const
    {
        return CO2::gasViscosity(temperature, pressure);
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& /*Rv*/) const
    { return saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the formation volume factor [-] of oil saturated gas at given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& temperature,
                                                     const Evaluation& pressure) const
    {
        return CO2::gasDensity(temperature, pressure)/gasReferenceDensity_[regionIdx];
    }

    /*!
     * \brief Returns the saturation pressure of the gas phase [Pa]
     *        depending on its mass fraction of the oil component
     *
     * \param Rv The surface volume of oil component dissolved in what will yield one cubic meter of gas at the surface [-]
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned /*regionIdx*/,
                                  const Evaluation& /*temperature*/,
                                  const Evaluation& /*Rv*/) const
    { return 0.0; /* this is dry gas! */ }

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of the oil phase.
     */
    template <class Evaluation>
    Evaluation saturatedOilVaporizationFactor(unsigned /*regionIdx*/,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& /*pressure*/,
                                              const Evaluation& /*oilSaturation*/,
                                              const Evaluation& /*maxOilSaturation*/) const
    { return 0.0; /* this is dry gas! */ }

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of the oil phase.
     */
    template <class Evaluation>
    Evaluation saturatedOilVaporizationFactor(unsigned /*regionIdx*/,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& /*pressure*/) const
    { return 0.0; /* this is dry gas! */ }

    const Scalar gasReferenceDensity(unsigned regionIdx) const
    { return gasReferenceDensity_[regionIdx]; }

    bool operator==(const Co2GasPvt<Scalar>& data) const
    {
        return gasReferenceDensity_ == data.gasReferenceDensity_;
    }

private:
    std::vector<Scalar> gasReferenceDensity_;
};

} // namespace Ewoms

#endif
