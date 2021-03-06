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
 * \copydoc Ewoms::GasPvtMultiplexer
 */
#ifndef EWOMS_GAS_PVT_MULTIPLEXER_HH
#define EWOMS_GAS_PVT_MULTIPLEXER_HH

#include "drygaspvt.hh"
#include "wetgaspvt.hh"
#include "gaspvtthermal.hh"
#include "co2gaspvt.hh"

#if HAVE_ECL_INPUT
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#endif

namespace Ewoms {
#define EWOMS_GAS_PVT_MULTIPLEXER_CALL(codeToCall)                        \
    switch (gasPvtApproach_) {                                          \
    case DryGasPvt: {                                                   \
        auto& pvtImpl = getRealPvt<DryGasPvt>();                        \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case WetGasPvt: {                                                   \
        auto& pvtImpl = getRealPvt<WetGasPvt>();                        \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case ThermalGasPvt: {                                               \
        auto& pvtImpl = getRealPvt<ThermalGasPvt>();                    \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case Co2GasPvt: {                                               \
        auto& pvtImpl = getRealPvt<Co2GasPvt>();                    \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case NoGasPvt:                                                      \
        throw std::logic_error("Not implemented: Gas PVT of this deck!"); \
    } \

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the gas
 *        phase in the black-oil model.
 *
 * This is a multiplexer class which forwards all calls to the real implementation.
 *
 * Note that, since the main application for this class is the black oil fluid system,
 * the API exposed by this class is pretty specific to the assumptions made by the black
 * oil model.
 */
template <class Scalar, bool enableThermal = true>
class GasPvtMultiplexer
{
public:
    typedef Ewoms::GasPvtThermal<Scalar> GasPvtThermal;

    enum GasPvtApproach {
        NoGasPvt,
        DryGasPvt,
        WetGasPvt,
        ThermalGasPvt,
        Co2GasPvt
    };

    GasPvtMultiplexer()
    {
        gasPvtApproach_ = NoGasPvt;
        realGasPvt_ = nullptr;
    }

    GasPvtMultiplexer(GasPvtApproach approach, void* realGasPvt)
        : gasPvtApproach_(approach)
        , realGasPvt_(realGasPvt)
    { }

    GasPvtMultiplexer(const GasPvtMultiplexer<Scalar,enableThermal>& data)
    {
        *this = data;
    }

    ~GasPvtMultiplexer()
    {
        switch (gasPvtApproach_) {
        case DryGasPvt: {
            delete &getRealPvt<DryGasPvt>();
            break;
        }
        case WetGasPvt: {
            delete &getRealPvt<WetGasPvt>();
            break;
        }
        case ThermalGasPvt: {
            delete &getRealPvt<ThermalGasPvt>();
            break;
        }
        case Co2GasPvt: {
            delete &getRealPvt<Co2GasPvt>();
            break;
        }
        case NoGasPvt:
            break;
        }
    }

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize the parameters for gas using an ECL deck.
     *
     * This method assumes that the deck features valid DENSITY and PVDG keywords.
     */
    void initFromEclState(const EclipseState& eclState, const Schedule& schedule)
    {
        if (!eclState.runspec().phases().active(Phase::GAS))
            return;
        if (eclState.runspec().co2Storage())
            setApproach(Co2GasPvt);
        else if (enableThermal && eclState.getSimulationConfig().isThermal())
            setApproach(ThermalGasPvt);
        else if (!eclState.getTableManager().getPvtgTables().empty())
            setApproach(WetGasPvt);
        else if (eclState.getTableManager().hasTables("PVDG"))
            setApproach(DryGasPvt);

        EWOMS_GAS_PVT_MULTIPLEXER_CALL(pvtImpl.initFromEclState(eclState, schedule));
    }
#endif // HAVE_ECL_INPUT

    void setApproach(GasPvtApproach gasPvtAppr)
    {
        switch (gasPvtAppr) {
        case DryGasPvt:
            realGasPvt_ = new Ewoms::DryGasPvt<Scalar>;
            break;

        case WetGasPvt:
            realGasPvt_ = new Ewoms::WetGasPvt<Scalar>;
            break;

        case ThermalGasPvt:
            realGasPvt_ = new Ewoms::GasPvtThermal<Scalar>;
            break;

        case Co2GasPvt:
            realGasPvt_ = new Ewoms::Co2GasPvt<Scalar>;
            break;

        case NoGasPvt:
            throw std::logic_error("Not implemented: Gas PVT of this deck!");
        }

        gasPvtApproach_ = gasPvtAppr;
    }

    void initEnd()
    { EWOMS_GAS_PVT_MULTIPLEXER_CALL(pvtImpl.initEnd()); }

    /*!
     * \brief Return the number of PVT regions which are considered by this PVT-object.
     */
    unsigned numRegions() const
    { EWOMS_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.numRegions()); return 1; }

    /*!
     * \brief Return the reference density which are considered by this PVT-object.
     */
    const Scalar gasReferenceDensity(unsigned regionIdx)
    { EWOMS_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.gasReferenceDensity(regionIdx)); return 2.; }

    /*!
     * \brief Returns the specific enthalpy [J/kg] of gas given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                        const Evaluation& temperature,
                        const Evaluation& pressure,
                        const Evaluation& Rv) const
    { EWOMS_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.internalEnergy(regionIdx, temperature, pressure, Rv)); return 0; }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation = Scalar>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& Rv) const
    { EWOMS_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.viscosity(regionIdx, temperature, pressure, Rv)); return 0; }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of oil saturated gas given a set of parameters.
     */
    template <class Evaluation = Scalar>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const
    { EWOMS_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedViscosity(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation = Scalar>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& Rv) const
    { EWOMS_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rv)); return 0; }

    /*!
     * \brief Returns the formation volume factor [-] of oil saturated gas given a set of parameters.
     */
    template <class Evaluation = Scalar>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& temperature,
                                                     const Evaluation& pressure) const
    { EWOMS_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of oil saturated gas.
     */
    template <class Evaluation = Scalar>
    Evaluation saturatedOilVaporizationFactor(unsigned regionIdx,
                                              const Evaluation& temperature,
                                              const Evaluation& pressure) const
    { EWOMS_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedOilVaporizationFactor(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of oil saturated gas.
     */
    template <class Evaluation = Scalar>
    Evaluation saturatedOilVaporizationFactor(unsigned regionIdx,
                                              const Evaluation& temperature,
                                              const Evaluation& pressure,
                                              const Evaluation& oilSaturation,
                                              const Evaluation& maxOilSaturation) const
    { EWOMS_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedOilVaporizationFactor(regionIdx, temperature, pressure, oilSaturation, maxOilSaturation)); return 0; }

    /*!
     * \brief Returns the saturation pressure of the gas phase [Pa]
     *        depending on its mass fraction of the oil component
     *
     * \param Rv The surface volume of oil component dissolved in what will yield one cubic meter of gas at the surface [-]
     */
    template <class Evaluation = Scalar>
    Evaluation saturationPressure(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& Rv) const
    { EWOMS_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturationPressure(regionIdx, temperature, Rv)); return 0; }

    /*!
     * \brief Returns the concrete approach for calculating the PVT relations.
     *
     * (This is only determined at runtime.)
     */
    GasPvtApproach gasPvtApproach() const
    { return gasPvtApproach_; }

    // get the parameter object for the dry gas case
    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == DryGasPvt, Ewoms::DryGasPvt<Scalar> >::type& getRealPvt()
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<Ewoms::DryGasPvt<Scalar>* >(realGasPvt_);
    }

    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == DryGasPvt, const Ewoms::DryGasPvt<Scalar> >::type& getRealPvt() const
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<const Ewoms::DryGasPvt<Scalar>* >(realGasPvt_);
    }

    // get the parameter object for the wet gas case
    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == WetGasPvt, Ewoms::WetGasPvt<Scalar> >::type& getRealPvt()
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<Ewoms::WetGasPvt<Scalar>* >(realGasPvt_);
    }

    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == WetGasPvt, const Ewoms::WetGasPvt<Scalar> >::type& getRealPvt() const
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<const Ewoms::WetGasPvt<Scalar>* >(realGasPvt_);
    }

    // get the parameter object for the thermal gas case
    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == ThermalGasPvt, Ewoms::GasPvtThermal<Scalar> >::type& getRealPvt()
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<Ewoms::GasPvtThermal<Scalar>* >(realGasPvt_);
    }

    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == ThermalGasPvt, const Ewoms::GasPvtThermal<Scalar> >::type& getRealPvt() const
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<const Ewoms::GasPvtThermal<Scalar>* >(realGasPvt_);
    }

    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == Co2GasPvt, Ewoms::Co2GasPvt<Scalar> >::type& getRealPvt()
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<Ewoms::Co2GasPvt<Scalar>* >(realGasPvt_);
    }

    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == Co2GasPvt, const Ewoms::Co2GasPvt<Scalar> >::type& getRealPvt() const
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<const Ewoms::Co2GasPvt<Scalar>* >(realGasPvt_);
    }

    const void* realGasPvt() const { return realGasPvt_; }

    bool operator==(const GasPvtMultiplexer<Scalar,enableThermal>& data) const
    {
        if (this->gasPvtApproach() != data.gasPvtApproach())
            return false;

        switch (gasPvtApproach_) {
        case DryGasPvt:
            return *static_cast<const Ewoms::DryGasPvt<Scalar>*>(realGasPvt_) ==
                   *static_cast<const Ewoms::DryGasPvt<Scalar>*>(data.realGasPvt_);
        case WetGasPvt:
            return *static_cast<const Ewoms::WetGasPvt<Scalar>*>(realGasPvt_) ==
                   *static_cast<const Ewoms::WetGasPvt<Scalar>*>(data.realGasPvt_);
        case ThermalGasPvt:
            return *static_cast<const Ewoms::GasPvtThermal<Scalar>*>(realGasPvt_) ==
                   *static_cast<const Ewoms::GasPvtThermal<Scalar>*>(data.realGasPvt_);
        case Co2GasPvt:
            return *static_cast<const Ewoms::Co2GasPvt<Scalar>*>(realGasPvt_) ==
                    *static_cast<const Ewoms::Co2GasPvt<Scalar>*>(data.realGasPvt_);
        default:
            return true;
        }
    }

    GasPvtMultiplexer<Scalar,enableThermal>& operator=(const GasPvtMultiplexer<Scalar,enableThermal>& data)
    {
        gasPvtApproach_ = data.gasPvtApproach_;
        switch (gasPvtApproach_) {
        case DryGasPvt:
            realGasPvt_ = new Ewoms::DryGasPvt<Scalar>(*static_cast<const Ewoms::DryGasPvt<Scalar>*>(data.realGasPvt_));
            break;
        case WetGasPvt:
            realGasPvt_ = new Ewoms::WetGasPvt<Scalar>(*static_cast<const Ewoms::WetGasPvt<Scalar>*>(data.realGasPvt_));
            break;
        case ThermalGasPvt:
            realGasPvt_ = new Ewoms::GasPvtThermal<Scalar>(*static_cast<const Ewoms::GasPvtThermal<Scalar>*>(data.realGasPvt_));
            break;
        case Co2GasPvt:
            realGasPvt_ = new Ewoms::Co2GasPvt<Scalar>(*static_cast<const Ewoms::Co2GasPvt<Scalar>*>(data.realGasPvt_));
            break;
        default:
            break;
        }

        return *this;
    }

private:
    GasPvtApproach gasPvtApproach_;
    void* realGasPvt_;
};

#undef EWOMS_GAS_PVT_MULTIPLEXER_CALL

} // namespace Ewoms

#endif
