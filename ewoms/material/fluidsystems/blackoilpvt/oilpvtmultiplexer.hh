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
 * \copydoc Ewoms::OilPvtMultiplexer
 */
#ifndef EWOMS_OIL_PVT_MULTIPLEXER_HH
#define EWOMS_OIL_PVT_MULTIPLEXER_HH

#include "constantcompressibilityoilpvt.hh"
#include "deadoilpvt.hh"
#include "liveoilpvt.hh"
#include "oilpvtthermal.hh"
#include "brineco2pvt.hh"

#if HAVE_ECL_INPUT
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/runspec.hh>
#endif

namespace Ewoms {
#define EWOMS_OIL_PVT_MULTIPLEXER_CALL(codeToCall)                        \
    switch (approach_) {                                                \
    case ConstantCompressibilityOilPvt: {                               \
        auto& pvtImpl = getRealPvt<ConstantCompressibilityOilPvt>();    \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case DeadOilPvt: {                                                  \
        auto& pvtImpl = getRealPvt<DeadOilPvt>();                       \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case LiveOilPvt: {                                                  \
        auto& pvtImpl = getRealPvt<LiveOilPvt>();                       \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case ThermalOilPvt: {                                               \
        auto& pvtImpl = getRealPvt<ThermalOilPvt>();                    \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case BrineCo2Pvt: {                                                 \
        auto& pvtImpl = getRealPvt<BrineCo2Pvt>();                      \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case NoOilPvt:                                                      \
        throw std::logic_error("Not implemented: Oil PVT of this deck!"); \
    }                                                                     \

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the oil
 *        phase in the black-oil model.
 *
 * This is the base class which which provides an API for the actual PVT implementation
 * classes which based on dynamic polymorphism. The rationale to use dynamic polymorphism
 * here is that this enables the fluid system to easily switch the used PVT relations for
 * the individual fluid phases.
 *
 * Note that, since the application for this class is the black-oil fluid system, the API
 * exposed by this class is pretty specific to the black-oil model.
 */
template <class Scalar, bool enableThermal = true>
class OilPvtMultiplexer
{
public:
    typedef Ewoms::OilPvtThermal<Scalar> OilPvtThermal;

    enum OilPvtApproach {
        NoOilPvt,
        LiveOilPvt,
        DeadOilPvt,
        ConstantCompressibilityOilPvt,
        ThermalOilPvt,
        BrineCo2Pvt
    };

    OilPvtMultiplexer()
    {
        approach_ = NoOilPvt;
        realOilPvt_ = nullptr;
    }

    OilPvtMultiplexer(OilPvtApproach approach, void* realOilPvt)
        : approach_(approach)
        , realOilPvt_(realOilPvt)
    { }

    OilPvtMultiplexer(const OilPvtMultiplexer<Scalar,enableThermal>& data)
    {
        *this = data;
    }

    ~OilPvtMultiplexer()
    {
        switch (approach_) {
        case LiveOilPvt: {
            delete &getRealPvt<LiveOilPvt>();
            break;
        }
        case DeadOilPvt: {
            delete &getRealPvt<DeadOilPvt>();
            break;
        }
        case ConstantCompressibilityOilPvt: {
            delete &getRealPvt<ConstantCompressibilityOilPvt>();
            break;
        }
        case ThermalOilPvt: {
            delete &getRealPvt<ThermalOilPvt>();
            break;
        }
        case BrineCo2Pvt: {
            delete &getRealPvt<BrineCo2Pvt>();
            break;
        }

        case NoOilPvt:
            break;
        }
    }

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize the parameters for water using an ECL state.
     *
     * This method assumes that the deck features valid DENSITY and PVTO/PVDO/PVCDO keywords.
     */
    void initFromEclState(const EclipseState& eclState, const Schedule& schedule)
    {
        if (!eclState.runspec().phases().active(Phase::OIL))
            return;
        // TODO move the BrineCo2 approach to the waterPvtMultiplexer
        // when a proper gas-water simulator is supported
        if (eclState.runspec().co2Storage())
            setApproach(BrineCo2Pvt);
        else if (enableThermal && eclState.getSimulationConfig().isThermal())
            setApproach(ThermalOilPvt);
        else if (!eclState.getTableManager().getPvcdoTable().empty())
            setApproach(ConstantCompressibilityOilPvt);
        else if (eclState.getTableManager().hasTables("PVDO"))
            setApproach(DeadOilPvt);
        else if (!eclState.getTableManager().getPvtoTables().empty())
            setApproach(LiveOilPvt);

        EWOMS_OIL_PVT_MULTIPLEXER_CALL(pvtImpl.initFromEclState(eclState, schedule));
    }
#endif // HAVE_ECL_INPUT

    void initEnd()
    { EWOMS_OIL_PVT_MULTIPLEXER_CALL(pvtImpl.initEnd()); }

    /*!
     * \brief Return the number of PVT regions which are considered by this PVT-object.
     */
    unsigned numRegions() const
    { EWOMS_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.numRegions()); return 1; }

    /*!
     * \brief Return the reference density which are considered by this PVT-object.
     */
    const Scalar oilReferenceDensity(unsigned regionIdx)
    { EWOMS_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.oilReferenceDensity(regionIdx)); return 700.; }

    /*!
     * \brief Returns the specific enthalpy [J/kg] oil given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                        const Evaluation& temperature,
                        const Evaluation& pressure,
                        const Evaluation& Rs) const
    { EWOMS_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.internalEnergy(regionIdx, temperature, pressure, Rs)); return 0; }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& Rs) const
    { EWOMS_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.viscosity(regionIdx, temperature, pressure, Rs)); return 0; }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const
    { EWOMS_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedViscosity(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& Rs) const
    { EWOMS_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rs)); return 0; }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& temperature,
                                                     const Evaluation& pressure) const
    { EWOMS_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of saturated oil.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure) const
    { EWOMS_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedGasDissolutionFactor(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of saturated oil.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure,
                                             const Evaluation& oilSaturation,
                                             const Evaluation& maxOilSaturation) const
    { EWOMS_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedGasDissolutionFactor(regionIdx, temperature, pressure, oilSaturation, maxOilSaturation)); return 0; }

    /*!
     * \brief Returns the saturation pressure [Pa] of oil given the mass fraction of the
     *        gas component in the oil phase.
     *
     * Calling this method only makes sense for live oil. All other implementations of
     * the black-oil PVT interface will just throw an exception...
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& Rs) const
    { EWOMS_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturationPressure(regionIdx, temperature, Rs)); return 0; }

    void setApproach(OilPvtApproach appr)
    {
        switch (appr) {
        case LiveOilPvt:
            realOilPvt_ = new Ewoms::LiveOilPvt<Scalar>;
            break;

        case DeadOilPvt:
            realOilPvt_ = new Ewoms::DeadOilPvt<Scalar>;
            break;

        case ConstantCompressibilityOilPvt:
            realOilPvt_ = new Ewoms::ConstantCompressibilityOilPvt<Scalar>;
            break;

        case ThermalOilPvt:
            realOilPvt_ = new Ewoms::OilPvtThermal<Scalar>;
            break;

        case BrineCo2Pvt:
            realOilPvt_ = new Ewoms::BrineCo2Pvt<Scalar>;
            break;

        case NoOilPvt:
            throw std::logic_error("Not implemented: Oil PVT of this deck!");
        }

        approach_ = appr;
    }

    /*!
     * \brief Returns the concrete approach for calculating the PVT relations.
     *
     * (This is only determined at runtime.)
     */
    OilPvtApproach approach() const
    { return approach_; }

    // get the concrete parameter object for the oil phase
    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == LiveOilPvt, Ewoms::LiveOilPvt<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<Ewoms::LiveOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == LiveOilPvt, const Ewoms::LiveOilPvt<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<Ewoms::LiveOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == DeadOilPvt, Ewoms::DeadOilPvt<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<Ewoms::DeadOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == DeadOilPvt, const Ewoms::DeadOilPvt<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<Ewoms::DeadOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == ConstantCompressibilityOilPvt, Ewoms::ConstantCompressibilityOilPvt<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<Ewoms::ConstantCompressibilityOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == ConstantCompressibilityOilPvt, const Ewoms::ConstantCompressibilityOilPvt<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<Ewoms::ConstantCompressibilityOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == ThermalOilPvt, Ewoms::OilPvtThermal<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<Ewoms::OilPvtThermal<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == ThermalOilPvt, const Ewoms::OilPvtThermal<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<const Ewoms::OilPvtThermal<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == BrineCo2Pvt, Ewoms::BrineCo2Pvt<Scalar>>::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<Ewoms::BrineCo2Pvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == BrineCo2Pvt, const Ewoms::BrineCo2Pvt<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<const Ewoms::BrineCo2Pvt<Scalar>* >(realOilPvt_);
    }

    const void* realOilPvt() const { return realOilPvt_; }

    bool operator==(const OilPvtMultiplexer<Scalar,enableThermal>& data) const
    {
        if (this->approach() != data.approach())
            return false;

        switch (approach_) {
        case ConstantCompressibilityOilPvt:
            return *static_cast<const Ewoms::ConstantCompressibilityOilPvt<Scalar>*>(realOilPvt_) ==
                   *static_cast<const Ewoms::ConstantCompressibilityOilPvt<Scalar>*>(data.realOilPvt_);
        case DeadOilPvt:
            return *static_cast<const Ewoms::DeadOilPvt<Scalar>*>(realOilPvt_) ==
                   *static_cast<const Ewoms::DeadOilPvt<Scalar>*>(data.realOilPvt_);
        case LiveOilPvt:
            return *static_cast<const Ewoms::LiveOilPvt<Scalar>*>(realOilPvt_) ==
                   *static_cast<const Ewoms::LiveOilPvt<Scalar>*>(data.realOilPvt_);
        case ThermalOilPvt:
            return *static_cast<const Ewoms::OilPvtThermal<Scalar>*>(realOilPvt_) ==
                   *static_cast<const Ewoms::OilPvtThermal<Scalar>*>(data.realOilPvt_);
        case BrineCo2Pvt:
            return *static_cast<const Ewoms::BrineCo2Pvt<Scalar>*>(realOilPvt_) ==
                    *static_cast<const Ewoms::BrineCo2Pvt<Scalar>*>(data.realOilPvt_);
        default:
            return true;
        }
    }

    OilPvtMultiplexer<Scalar,enableThermal>& operator=(const OilPvtMultiplexer<Scalar,enableThermal>& data)
    {
        approach_ = data.approach_;
        switch (approach_) {
        case ConstantCompressibilityOilPvt:
            realOilPvt_ = new Ewoms::ConstantCompressibilityOilPvt<Scalar>(*static_cast<const Ewoms::ConstantCompressibilityOilPvt<Scalar>*>(data.realOilPvt_));
            break;
        case DeadOilPvt:
            realOilPvt_ = new Ewoms::DeadOilPvt<Scalar>(*static_cast<const Ewoms::DeadOilPvt<Scalar>*>(data.realOilPvt_));
            break;
        case LiveOilPvt:
            realOilPvt_ = new Ewoms::LiveOilPvt<Scalar>(*static_cast<const Ewoms::LiveOilPvt<Scalar>*>(data.realOilPvt_));
            break;
        case ThermalOilPvt:
            realOilPvt_ = new Ewoms::OilPvtThermal<Scalar>(*static_cast<const Ewoms::OilPvtThermal<Scalar>*>(data.realOilPvt_));
            break;
        case BrineCo2Pvt:
            realOilPvt_ = new Ewoms::BrineCo2Pvt<Scalar>(*static_cast<const Ewoms::BrineCo2Pvt<Scalar>*>(data.realOilPvt_));
            break;
        default:
            break;
        }

        return *this;
    }

private:
    OilPvtApproach approach_;
    void* realOilPvt_;
};

#undef EWOMS_OIL_PVT_MULTIPLEXER_CALL

} // namespace Ewoms

#endif
