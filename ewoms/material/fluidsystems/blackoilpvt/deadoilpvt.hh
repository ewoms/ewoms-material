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
 * \copydoc Ewoms::DeadOilPvt
 */
#ifndef EWOMS_DEAD_OIL_PVT_HH
#define EWOMS_DEAD_OIL_PVT_HH

#include <ewoms/common/uniformxtabulated2dfunction.hh>
#include <ewoms/common/tabulated1dfunction.hh>
#include <ewoms/common/spline.hh>

#if HAVE_ECL_INPUT
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/pvdotable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tablemanager.hh>
#endif

namespace Ewoms {
/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the oil phase
 *        without dissolved gas.
 */
template <class Scalar>
class DeadOilPvt
{

public:
    typedef Ewoms::Tabulated1DFunction<Scalar> TabulatedOneDFunction;

    DeadOilPvt() = default;
    DeadOilPvt(const std::vector<Scalar>& oilReferenceDensity,
               const std::vector<TabulatedOneDFunction>& inverseOilB,
               const std::vector<TabulatedOneDFunction>& oilMu,
               const std::vector<TabulatedOneDFunction>& inverseOilBMu)
        : oilReferenceDensity_(oilReferenceDensity)
        , inverseOilB_(inverseOilB)
        , oilMu_(oilMu)
        , inverseOilBMu_(inverseOilBMu)
    { }

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize the oil parameters via the data specified by the PVDO ECL keyword.
     */
    void initFromEclState(const EclipseState& eclState, const Schedule&)
    {
        const auto& pvdoTables = eclState.getTableManager().getPvdoTables();
        const auto& densityTable = eclState.getTableManager().getDensityTable();

        assert(pvdoTables.size() == densityTable.size());

        size_t numRegions = pvdoTables.size();
        setNumRegions(numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            Scalar rhoRefO = densityTable[regionIdx].oil;
            Scalar rhoRefG = densityTable[regionIdx].gas;
            Scalar rhoRefW = densityTable[regionIdx].water;

            setReferenceDensities(regionIdx, rhoRefO, rhoRefG, rhoRefW);

            const auto& pvdoTable = pvdoTables.getTable<PvdoTable>(regionIdx);

            const auto& BColumn(pvdoTable.getFormationFactorColumn());
            std::vector<Scalar> invBColumn(BColumn.size());
            for (unsigned i = 0; i < invBColumn.size(); ++i)
                invBColumn[i] = 1/BColumn[i];

            inverseOilB_[regionIdx].setXYArrays(pvdoTable.numRows(),
                                                pvdoTable.getPressureColumn(),
                                                invBColumn);
            oilMu_[regionIdx].setXYArrays(pvdoTable.numRows(),
                                          pvdoTable.getPressureColumn(),
                                          pvdoTable.getViscosityColumn());
        }

        initEnd();
    }
#endif // HAVE_ECL_INPUT

    void setNumRegions(size_t numRegions)
    {
        oilReferenceDensity_.resize(numRegions);
        inverseOilB_.resize(numRegions);
        inverseOilBMu_.resize(numRegions);
        oilMu_.resize(numRegions);
    }

    /*!
     * \brief Initialize the reference densities of all fluids for a given PVT region
     */
    void setReferenceDensities(unsigned regionIdx,
                               Scalar rhoRefOil,
                               Scalar /*rhoRefGas*/,
                               Scalar /*rhoRefWater*/)
    {
        oilReferenceDensity_[regionIdx] = rhoRefOil;
    }

    /*!
     * \brief Initialize the function for the oil formation volume factor
     *
     * The oil formation volume factor \f$B_o\f$ is a function of \f$(p_o, X_o^G)\f$ and
     * represents the partial density of the oil component in the oil phase at a given
     * pressure.
     *
     * This method sets \f$1/B_o(p_o)\f$. Note that the mass fraction of the gas
     * component in the oil phase is missing when assuming dead oil.
     */
    void setInverseOilFormationVolumeFactor(unsigned regionIdx, const TabulatedOneDFunction& invBo)
    { inverseOilB_[regionIdx] = invBo; }

    /*!
     * \brief Initialize the viscosity of the oil phase.
     *
     * This is a function of \f$(R_s, p_o)\f$...
     */
    void setOilViscosity(unsigned regionIdx, const TabulatedOneDFunction& muo)
    { oilMu_[regionIdx] = muo; }

    /*!
     * \brief Finish initializing the oil phase PVT properties.
     */
    void initEnd()
    {
        // calculate the final 2D functions which are used for interpolation.
        size_t numRegions = oilMu_.size();
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            // calculate the table which stores the inverse of the product of the oil
            // formation volume factor and the oil viscosity
            const auto& oilMu = oilMu_[regionIdx];
            const auto& invOilB = inverseOilB_[regionIdx];
            assert(oilMu.numSamples() == invOilB.numSamples());

            std::vector<Scalar> invBMuColumn;
            std::vector<Scalar> pressureColumn;
            invBMuColumn.resize(oilMu.numSamples());
            pressureColumn.resize(oilMu.numSamples());

            for (unsigned pIdx = 0; pIdx < oilMu.numSamples(); ++pIdx) {
                pressureColumn[pIdx] = invOilB.xAt(pIdx);
                invBMuColumn[pIdx] = invOilB.valueAt(pIdx)*1/oilMu.valueAt(pIdx);
            }

            inverseOilBMu_[regionIdx].setXYArrays(pressureColumn.size(),
                                                  pressureColumn,
                                                  invBMuColumn);
        }
    }

    /*!
     * \brief Return the number of PVT regions which are considered by this PVT-object.
     */
    unsigned numRegions() const
    { return inverseOilBMu_.size(); }

    /*!
     * \brief Returns the specific enthalpy [J/kg] of oil given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx EWOMS_UNUSED,
                        const Evaluation& temperature EWOMS_UNUSED,
                        const Evaluation& pressure EWOMS_UNUSED,
                        const Evaluation& Rs EWOMS_UNUSED) const
    {
        throw std::runtime_error("Requested the enthalpy of oil but the thermal option is not enabled");
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& /*Rs*/) const
    { return saturatedViscosity(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of gas saturated oil given a pressure.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& /*temperature*/,
                                  const Evaluation& pressure) const
    {
        const Evaluation& invBo = inverseOilB_[regionIdx].eval(pressure, /*extrapolate=*/true);
        const Evaluation& invMuoBo = inverseOilBMu_[regionIdx].eval(pressure, /*extrapolate=*/true);

        return invBo/invMuoBo;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& /*temperature*/,
                                            const Evaluation& pressure,
                                            const Evaluation& /*Rs*/) const
    { return inverseOilB_[regionIdx].eval(pressure, /*extrapolate=*/true); }

    /*!
     * \brief Returns the formation volume factor [-] of saturated oil.
     *
     * Note that by definition, dead oil is always gas saturated.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& pressure) const
    { return inverseOilB_[regionIdx].eval(pressure, /*extrapolate=*/true); }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned /*regionIdx*/,
                                             const Evaluation& /*temperature*/,
                                             const Evaluation& /*pressure*/) const
    { return 0.0; /* this is dead oil! */ }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned /*regionIdx*/,
                                             const Evaluation& /*temperature*/,
                                             const Evaluation& /*pressure*/,
                                             const Evaluation& /*oilSaturation*/,
                                             const Evaluation& /*maxOilSaturation*/) const
    { return 0.0; /* this is dead oil! */ }

    /*!
     * \brief Returns the saturation pressure of the oil phase [Pa]
     *        depending on its mass fraction of the gas component
     *
     * \param Rs The surface volume of gas component dissolved in what will yield one cubic meter of oil at the surface [-]
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned /*regionIdx*/,
                                  const Evaluation& /*temperature*/,
                                  const Evaluation& /*Rs*/) const
    { return 0.0; /* this is dead oil, so there isn't any meaningful saturation pressure! */ }

    template <class Evaluation>
    Evaluation saturatedGasMassFraction(unsigned /*regionIdx*/,
                                        const Evaluation& /*temperature*/,
                                        const Evaluation& /*pressure*/) const
    { return 0.0; /* this is dead oil! */ }

    template <class Evaluation>
    Evaluation saturatedGasMoleFraction(unsigned /*regionIdx*/,
                                        const Evaluation& /*temperature*/,
                                        const Evaluation& /*pressure*/) const
    { return 0.0; /* this is dead oil! */ }

    const Scalar oilReferenceDensity(unsigned regionIdx) const
    { return oilReferenceDensity_[regionIdx]; }

    const std::vector<TabulatedOneDFunction>& inverseOilB() const
    { return inverseOilB_; }

    const std::vector<TabulatedOneDFunction>& oilMu() const
    { return oilMu_; }

    const std::vector<TabulatedOneDFunction>& inverseOilBMu() const
    { return inverseOilBMu_; }

    bool operator==(const DeadOilPvt<Scalar>& data) const
    {
        return this->oilReferenceDensity() == data.oilReferenceDensity() &&
               this->inverseOilB() == data.inverseOilB() &&
               this->oilMu() == data.oilMu() &&
               this->inverseOilBMu() == data.inverseOilBMu();
    }

private:
    std::vector<Scalar> oilReferenceDensity_;
    std::vector<TabulatedOneDFunction> inverseOilB_;
    std::vector<TabulatedOneDFunction> oilMu_;
    std::vector<TabulatedOneDFunction> inverseOilBMu_;
};

} // namespace Ewoms

#endif
