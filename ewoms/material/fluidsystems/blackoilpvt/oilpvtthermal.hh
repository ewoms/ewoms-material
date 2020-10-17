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
 * \copydoc Ewoms::OilPvtThermal
 */
#ifndef EWOMS_OIL_PVT_THERMAL_HH
#define EWOMS_OIL_PVT_THERMAL_HH

#include <ewoms/material/constants.hh>

#include <ewoms/common/final.hh>
#include <ewoms/common/uniformxtabulated2dfunction.hh>
#include <ewoms/common/tabulated1dfunction.hh>
#include <ewoms/common/spline.hh>

#if HAVE_ECL_INPUT
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/simpletable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tablemanager.hh>
#endif

namespace Ewoms {
template <class Scalar, bool enableThermal>
class OilPvtMultiplexer;

/*!
 * \brief This class implements temperature dependence of the PVT properties of oil
 *
 * Note that this _only_ implements the temperature part, i.e., it requires the
 * isothermal properties as input.
 */
template <class Scalar>
class OilPvtThermal
{
public:
    typedef Ewoms::Tabulated1DFunction<Scalar> TabulatedOneDFunction;
    typedef OilPvtMultiplexer<Scalar, /*enableThermal=*/false> IsothermalPvt;

    OilPvtThermal()
    {
        enableThermalDensity_ = false;
        enableThermalViscosity_ = false;
        enableInternalEnergy_ = false;
        isothermalPvt_ = nullptr;
    }

    OilPvtThermal(IsothermalPvt* isothermalPvt,
                  const std::vector<TabulatedOneDFunction>& oilvisctCurves,
                  const std::vector<Scalar>& viscrefPress,
                  const std::vector<Scalar>& viscrefRs,
                  const std::vector<Scalar>& viscRef,
                  const std::vector<Scalar>& oildentRefTemp,
                  const std::vector<Scalar>& oildentCT1,
                  const std::vector<Scalar>& oildentCT2,
                  const std::vector<TabulatedOneDFunction>& internalEnergyCurves,
                  bool enableThermalDensity,
                  bool enableThermalViscosity,
                  bool enableInternalEnergy)
        : isothermalPvt_(isothermalPvt)
        , oilvisctCurves_(oilvisctCurves)
        , viscrefPress_(viscrefPress)
        , viscrefRs_(viscrefRs)
        , viscRef_(viscRef)
        , oildentRefTemp_(oildentRefTemp)
        , oildentCT1_(oildentCT1)
        , oildentCT2_(oildentCT2)
        , internalEnergyCurves_(internalEnergyCurves)
        , enableThermalDensity_(enableThermalDensity)
        , enableThermalViscosity_(enableThermalViscosity)
        , enableInternalEnergy_(enableInternalEnergy)
    { }

    OilPvtThermal(const OilPvtThermal& data)
    { *this = data; }

    ~OilPvtThermal()
    { delete isothermalPvt_; }

#if HAVE_ECL_INPUT
    /*!
     * \brief Implement the temperature part of the oil PVT properties.
     */
    void initFromEclState(const EclipseState& eclState, const Schedule& schedule)
    {
        //////
        // initialize the isothermal part
        //////
        isothermalPvt_ = new IsothermalPvt;
        isothermalPvt_->initFromEclState(eclState, schedule);

        //////
        // initialize the thermal part
        //////
        const auto& tables = eclState.getTableManager();

        enableThermalDensity_ = tables.OilDenT().size() > 0;
        enableThermalViscosity_ = tables.hasTables("OILVISCT");
        enableInternalEnergy_ = tables.hasTables("SPECHEAT");

        unsigned numRegions = isothermalPvt_->numRegions();
        setNumRegions(numRegions);

        // viscosity
        if (enableThermalViscosity_) {
            if (tables.getViscrefTable().empty())
                throw std::runtime_error("VISCREF is required when OILVISCT is present");

            const auto& oilvisctTables = tables.getOilvisctTables();
            const auto& viscrefTable = tables.getViscrefTable();

            assert(oilvisctTables.size() == numRegions);
            assert(viscrefTable.size() == numRegions);

            for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const auto& TCol = oilvisctTables[regionIdx].getColumn("Temperature").vectorCopy();
                const auto& muCol = oilvisctTables[regionIdx].getColumn("Viscosity").vectorCopy();
                oilvisctCurves_[regionIdx].setXYContainers(TCol, muCol);

                viscrefPress_[regionIdx] = viscrefTable[regionIdx].reference_pressure;
                viscrefRs_[regionIdx] = viscrefTable[regionIdx].reference_rs;

                // temperature used to calculate the reference viscosity [K]. the
                // value does not really matter if the underlying PVT object really
                // is isothermal...
                Scalar Tref = 273.15 + 20;

                // compute the reference viscosity using the isothermal PVT object.
                viscRef_[regionIdx] =
                    isothermalPvt_->viscosity(regionIdx,
                                              Tref,
                                              viscrefPress_[regionIdx],
                                              viscrefRs_[regionIdx]);
            }
        }

        // temperature dependence of oil density
        const auto& oilDenT = tables.OilDenT();
        if (oilDenT.size() > 0) {
            assert(oilDenT.size() == numRegions);
            for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const auto& record = oilDenT[regionIdx];

                oildentRefTemp_[regionIdx] = record.T0;
                oildentCT1_[regionIdx] = record.C1;
                oildentCT2_[regionIdx] = record.C2;
            }
        }

        if (enableInternalEnergy_) {
            // the specific internal energy of liquid oil. be aware that ecl only specifies the
            // heat capacity (via the SPECHEAT keyword) and we need to integrate it
            // ourselfs to get the internal energy
            for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const auto& specheatTable = tables.getSpecheatTables()[regionIdx];
                const auto& temperatureColumn = specheatTable.getColumn("TEMPERATURE");
                const auto& cvOilColumn = specheatTable.getColumn("CV_OIL");

                std::vector<double> uSamples(temperatureColumn.size());

                Scalar u = temperatureColumn[0]*cvOilColumn[0];
                for (size_t i = 0;; ++i) {
                    uSamples[i] = u;

                    if (i >= temperatureColumn.size() - 1)
                        break;

                    // integrate to the heat capacity from the current sampling point to the next
                    // one. this leads to a quadratic polynomial.
                    Scalar c_v0 = cvOilColumn[i];
                    Scalar c_v1 = cvOilColumn[i + 1];
                    Scalar T0 = temperatureColumn[i];
                    Scalar T1 = temperatureColumn[i + 1];
                    u += 0.5*(c_v0 + c_v1)*(T1 - T0);
                }

                internalEnergyCurves_[regionIdx].setXYContainers(temperatureColumn.vectorCopy(), uSamples);
            }
        }
    }
#endif // HAVE_ECL_INPUT

    /*!
     * \brief Set the number of PVT-regions considered by this object.
     */
    void setNumRegions(size_t numRegions)
    {
        oilvisctCurves_.resize(numRegions);
        viscrefPress_.resize(numRegions);
        viscrefRs_.resize(numRegions);
        viscRef_.resize(numRegions);
        internalEnergyCurves_.resize(numRegions);
    }

    /*!
     * \brief Finish initializing the thermal part of the oil phase PVT properties.
     */
    void initEnd()
    { }

    /*!
     * \brief Returns true iff the density of the oil phase is temperature dependent.
     */
    bool enableThermalDensity() const
    { return enableThermalDensity_; }

    /*!
     * \brief Returns true iff the viscosity of the oil phase is temperature dependent.
     */
    bool enableThermalViscosity() const
    { return enableThermalViscosity_; }

    size_t numRegions() const
    { return viscrefRs_.size(); }

    /*!
     * \brief Returns the specific internal energy [J/kg] of oil given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                              const Evaluation& temperature,
                              const Evaluation& pressure EWOMS_UNUSED,
                              const Evaluation& Rs EWOMS_UNUSED) const
    {
        if (!enableInternalEnergy_)
            throw std::runtime_error("Requested the internal energy of oil but it is disabled");

        // compute the specific internal energy for the specified tempature. We use linear
        // interpolation here despite the fact that the underlying heat capacities are
        // piecewise linear (which leads to a quadratic function)
        return internalEnergyCurves_[regionIdx].eval(temperature, /*extrapolate=*/true);
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& Rs) const
    {
        const auto& isothermalMu = isothermalPvt_->viscosity(regionIdx, temperature, pressure, Rs);
        if (!enableThermalViscosity())
            return isothermalMu;

        // compute the viscosity deviation due to temperature
        const auto& muOilvisct = oilvisctCurves_[regionIdx].eval(temperature);
        return muOilvisct/viscRef_[regionIdx]*isothermalMu;
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const
    {
        const auto& isothermalMu = isothermalPvt_->saturatedViscosity(regionIdx, temperature, pressure);
        if (!enableThermalViscosity())
            return isothermalMu;

        // compute the viscosity deviation due to temperature
        const auto& muOilvisct = oilvisctCurves_[regionIdx].eval(temperature, true);
        return muOilvisct/viscRef_[regionIdx]*isothermalMu;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& Rs) const
    {
        const auto& b =
            isothermalPvt_->inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rs);

        if (!enableThermalDensity())
            return b;

        // we use the same approach as for the for water here, but with the eWoms-specific
        // OILDENT keyword.
        Scalar TRef = oildentRefTemp_[regionIdx];
        Scalar cT1 = oildentCT1_[regionIdx];
        Scalar cT2 = oildentCT2_[regionIdx];
        const Evaluation& Y = temperature - TRef;

        return b/(1 + (cT1 + cT2*Y)*Y);
    }

    /*!
     * \brief Returns the formation volume factor [-] of gas-saturated oil phase.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& temperature,
                                                     const Evaluation& pressure) const
    {
        const auto& b =
            isothermalPvt_->saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure);

        if (!enableThermalDensity())
            return b;

        // we use the same approach as for the for water here, but with the eWoms-specific
        // OILDENT keyword.
        Scalar TRef = oildentRefTemp_[regionIdx];
        Scalar cT1 = oildentCT1_[regionIdx];
        Scalar cT2 = oildentCT2_[regionIdx];
        const Evaluation& Y = temperature - TRef;

        return b/(1 + (cT1 + cT2*Y)*Y);
    }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     *
     * This method implements temperature dependence and requires the isothermal gas
     * dissolution factor for gas saturated oil and temperature as inputs. Currently it
     * is just a dummy method which passes through the isothermal gas dissolution factor.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure) const
    { return isothermalPvt_->saturatedGasDissolutionFactor(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     *
     * This method implements temperature dependence and requires the isothermal gas
     * dissolution factor for gas saturated oil and temperature as inputs. Currently it
     * is just a dummy method which passes through the isothermal gas dissolution factor.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure,
                                             const Evaluation& oilSaturation,
                                             const Evaluation& maxOilSaturation) const
    { return isothermalPvt_->saturatedGasDissolutionFactor(regionIdx, temperature, pressure, oilSaturation, maxOilSaturation); }

    /*!
     * \brief Returns the saturation pressure of the oil phase [Pa]
     *
     * This method implements temperature dependence and requires isothermal satuation
     * pressure and temperature as inputs. Currently it is just a dummy method which
     * passes through the isothermal saturation pressure.
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const
    { return isothermalPvt_->saturationPressure(regionIdx, temperature, pressure); }

    const IsothermalPvt* isoThermalPvt() const
    { return isothermalPvt_; }

    const Scalar oilReferenceDensity(unsigned regionIdx) const
    { return isothermalPvt_->oilReferenceDensity(regionIdx); }

    const std::vector<TabulatedOneDFunction>& oilvisctCurves() const
    { return oilvisctCurves_; }

    const std::vector<Scalar>& viscrefPress() const
    { return viscrefPress_; }

    const std::vector<Scalar>& viscrefRs() const
    { return viscrefRs_; }

    const std::vector<Scalar>& viscRef() const
    { return viscRef_; }

    const std::vector<Scalar>& oildentRefTemp() const
    { return oildentRefTemp_; }

    const std::vector<Scalar>& oildentCT1() const
    { return oildentCT1_; }

    const std::vector<Scalar>& oildentCT2() const
    { return oildentCT2_; }

    const std::vector<TabulatedOneDFunction> internalEnergyCurves() const
    { return internalEnergyCurves_; }

    bool enableInternalEnergy() const
    { return enableInternalEnergy_; }

    bool operator==(const OilPvtThermal<Scalar>& data) const
    {
        if (isothermalPvt_ && !data.isothermalPvt_)
            return false;
        if (!isothermalPvt_ && data.isothermalPvt_)
            return false;

        return (!this->isoThermalPvt() ||
                (*this->isoThermalPvt() == *data.isoThermalPvt())) &&
                this->oilvisctCurves() == data.oilvisctCurves() &&
                this->viscrefPress() == data.viscrefPress() &&
                this->viscrefRs() == data.viscrefRs() &&
                this->viscRef() == data.viscRef() &&
                this->oildentRefTemp() == data.oildentRefTemp() &&
                this->oildentCT1() == data.oildentCT1() &&
                this->oildentCT2() == data.oildentCT2() &&
                this->internalEnergyCurves() == data.internalEnergyCurves() &&
                this->enableThermalDensity() == data.enableThermalDensity() &&
                this->enableThermalViscosity() == data.enableThermalViscosity() &&
                this->enableInternalEnergy() == data.enableInternalEnergy();
    }

    OilPvtThermal<Scalar>& operator=(const OilPvtThermal<Scalar>& data)
    {
        if (data.isothermalPvt_)
            isothermalPvt_ = new IsothermalPvt(*data.isothermalPvt_);
        else
            isothermalPvt_ = nullptr;
        oilvisctCurves_ = data.oilvisctCurves_;
        viscrefPress_ = data.viscrefPress_;
        viscrefRs_ = data.viscrefRs_;
        viscRef_ = data.viscRef_;
        oildentRefTemp_ = data.oildentRefTemp_;
        oildentCT1_ = data.oildentCT1_;
        oildentCT2_ = data.oildentCT2_;
        internalEnergyCurves_ = data.internalEnergyCurves_;
        enableThermalDensity_ = data.enableThermalDensity_;
        enableThermalViscosity_ = data.enableThermalViscosity_;
        enableInternalEnergy_ = data.enableInternalEnergy_;

        return *this;
    }

private:
    IsothermalPvt* isothermalPvt_;

    // The PVT properties needed for temperature dependence of the viscosity. We need
    // to store one value per PVT region.
    std::vector<TabulatedOneDFunction> oilvisctCurves_;
    std::vector<Scalar> viscrefPress_;
    std::vector<Scalar> viscrefRs_;
    std::vector<Scalar> viscRef_;

    // The PVT properties needed for temperature dependence of the density.
    std::vector<Scalar> oildentRefTemp_;
    std::vector<Scalar> oildentCT1_;
    std::vector<Scalar> oildentCT2_;

    // piecewise linear curve representing the internal energy of oil
    std::vector<TabulatedOneDFunction> internalEnergyCurves_;

    bool enableThermalDensity_;
    bool enableThermalViscosity_;
    bool enableInternalEnergy_;
};

} // namespace Ewoms

#endif
