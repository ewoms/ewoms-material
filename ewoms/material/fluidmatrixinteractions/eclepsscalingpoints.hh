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
 * \copydoc Ewoms::EclEpsTwoPhaseLawPoints
 */
#ifndef EWOMS_ECL_EPS_SCALING_POINTS_HH
#define EWOMS_ECL_EPS_SCALING_POINTS_HH

#include "eclepsconfig.hh"
#include "eclepsgridproperties.hh"

#if HAVE_ECL_INPUT
#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/deck/deckrecord.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/sgfntable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/sgoftable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/slgoftable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/sof2table.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/sof3table.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/swfntable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/swoftable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tablemanager.hh>
#endif

#include <ewoms/common/means.hh>

#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <algorithm>

namespace Ewoms {

/*!
 * \brief This structure represents all values which can be possibly used as scaling
 *        points in the endpoint scaling code.
 *
 * Depending on the exact configuration, some of these quantities are not used as actual
 * scaling points. It is easier to extract all of them at once, though.
 */
template <class Scalar>
struct EclEpsScalingPointsInfo
{
    // connate saturations
    Scalar Swl; // oil
    Scalar Sgl; // gas
    Scalar Sowl; // oil for the oil-water system
    Scalar Sogl; // oil for the gas-oil system

    // critical water and gas saturations
    Scalar krCriticalEps; // relative permeability below which a saturation is considered
                          // to be critical
    Scalar Swcr; // oil
    Scalar Sgcr; // gas
    Scalar Sowcr; // oil for the oil-water system
    Scalar Sogcr; // oil for the gas-oil system

    // maximum saturations
    Scalar Swu; // oil
    Scalar Sgu; // gas
    Scalar Sowu; // oil for the oil-water system
    Scalar Sogu; // oil for the gas-oil system

    // maximum capillary pressures
    Scalar maxPcow; // maximum capillary pressure of the oil-water system
    Scalar maxPcgo; // maximum capillary pressure of the gas-oil system

    // the Leverett capillary pressure scaling factors. (those only make sense for the
    // scaled points, for the unscaled ones they are 1.0.)
    Scalar pcowLeverettFactor;
    Scalar pcgoLeverettFactor;

    // maximum relative permabilities
    Scalar maxKrw; // maximum relative permability of water
    Scalar maxKrow; // maximum relative permability of oil in the oil-water system
    Scalar maxKrog; // maximum relative permability of oil in the gas-oil system
    Scalar maxKrg; // maximum relative permability of gas

    bool operator==(const EclEpsScalingPointsInfo<Scalar>& data) const
    {
        return Swl == data.Swl &&
               Sgl == data.Sgl &&
               Sowl == data.Sowl &&
               Sogl == data.Sogl &&
               krCriticalEps == data.krCriticalEps &&
               Swcr == data.Swcr &&
               Sgcr == data.Sgcr &&
               Sowcr == data.Sowcr &&
               Sogcr == data.Sogcr &&
               Swu == data.Swu &&
               Sgu == data.Sgu &&
               Sowu == data.Sowu &&
               Sogu == data.Sogu &&
               maxPcow == data.maxPcow &&
               maxPcgo == data.maxPcgo &&
               pcowLeverettFactor == data.pcowLeverettFactor &&
               pcgoLeverettFactor == data.pcgoLeverettFactor &&
               maxKrw == data.maxKrw &&
               maxKrow == data.maxKrow &&
               maxKrog == data.maxKrog &&
               maxKrg == data.maxKrg;
    }

    void print() const
    {
        std::cout << "    Swl: " << Swl << "\n"
                  << "    Sgl: " << Sgl << "\n"
                  << "    Sowl: " << Sowl << "\n"
                  << "    Sogl: " << Sogl << "\n"
                  << "    Swcr: " << Swcr << "\n"
                  << "    Sgcr: " << Sgcr << "\n"
                  << "    Sowcr: " << Sowcr << "\n"
                  << "    Sogcr: " << Sogcr << "\n"
                  << "    Swu: " << Swu << "\n"
                  << "    Sgu: " << Sgu << "\n"
                  << "    Sowu: " << Sowu << "\n"
                  << "    Sogu: " << Sogu << "\n"
                  << "    maxPcow: " << maxPcow << "\n"
                  << "    maxPcgo: " << maxPcgo << "\n"
                  << "    pcowLeverettFactor: " << pcowLeverettFactor << "\n"
                  << "    pcgoLeverettFactor: " << pcgoLeverettFactor << "\n"
                  << "    maxKrw: " << maxKrw << "\n"
                  << "    maxKrg: " << maxKrg << "\n"
                  << "    maxKrow: " << maxKrow << "\n"
                  << "    maxKrog: " << maxKrog << "\n";
    }

#if HAVE_ECL_INPUT
    /*!
     * \brief Extract the values of the unscaled scaling parameters.
     *
     * I.e., the values which are used for the nested Fluid-Matrix interactions and which
     * are produced by them.
     */
    void extractUnscaled(const Ewoms::EclipseState& eclState,
                         unsigned satRegionIdx)
    {
        // determine the value of the relative permeability below which the corresponding
        // saturation is considered to be critical
        const auto& satFuncCtrls = eclState.runspec().saturationFunctionControls();
        krCriticalEps = satFuncCtrls.minimumRelpermMobilityThreshold();

        const auto& tables = eclState.getTableManager();
        const TableContainer&  swofTables = tables.getSwofTables();
        const TableContainer&  sgofTables = tables.getSgofTables();
        const TableContainer& slgofTables = tables.getSlgofTables();
        const TableContainer&  swfnTables = tables.getSwfnTables();
        const TableContainer&  sgfnTables = tables.getSgfnTables();
        const TableContainer&  sof3Tables = tables.getSof3Tables();
        const TableContainer&  sof2Tables = tables.getSof2Tables();

        bool hasWater = eclState.runspec().phases().active(Phase::WATER);
        bool hasGas = eclState.runspec().phases().active(Phase::GAS);
        bool hasOil = eclState.runspec().phases().active(Phase::OIL);

        if (int(hasWater) + int(hasGas) + int(hasOil) == 1) {
            return;
        }
        else if (!hasWater) {
            Swl = 0.0;
            Swu = 0.0;
            Swcr = 0.0;
            bool family1 = (!sgofTables.empty() || !slgofTables.empty());
            bool family2 = !sgfnTables.empty() && !sof2Tables.empty();
            if (family1) {
                if (!sgofTables.empty())
                    extractUnscaledSgof_(sgofTables.getTable<SgofTable>(satRegionIdx));
                else {
                    assert(!slgofTables.empty());
                    extractUnscaledSlgof_(slgofTables.getTable<SlgofTable>(satRegionIdx));
                }
            }
            else if (family2) {
                extractUnscaledSgfn_(sgfnTables.getTable<SgfnTable>(satRegionIdx));
                extractUnscaledSof2_(sof2Tables.getTable<Sof2Table>(satRegionIdx));
            }
            else {
                throw std::domain_error("No valid saturation keyword family specified");
            }
            return;
        }
        else if (!hasGas) {
            Sgl = 0.0;
            Sgu = 0.0;
            Sgcr = 0.0;
            bool family1 = !swofTables.empty();
            bool family2 = !swfnTables.empty() && !sof2Tables.empty();
            if (family1) {
                extractUnscaledSwof_(swofTables.getTable<SwofTable>(satRegionIdx));
            }
            else if (family2) {
                extractUnscaledSwfn_(swfnTables.getTable<SwfnTable>(satRegionIdx));
                extractUnscaledSof2_(sof2Tables.getTable<Sof2Table>(satRegionIdx));
            }
            else {
                throw std::domain_error("No valid saturation keyword family specified");
            }
            return;
        }

        bool family1 = (!sgofTables.empty() || !slgofTables.empty()) && !swofTables.empty();
        bool family2 = !swfnTables.empty() && !sgfnTables.empty() && !sof3Tables.empty();

        // so far, only water-oil and oil-gas simulations are supported, i.e.,
        // there's no gas-water yet.
        if (!hasWater || !hasGas || !hasOil)
            throw std::domain_error("The specified phase configuration is not suppored");

        if (family1) {
            extractUnscaledSwof_(swofTables.getTable<SwofTable>(satRegionIdx));

            if (!sgofTables.empty()) {
                // gas-oil parameters are specified using the SGOF keyword
                extractUnscaledSgof_(sgofTables.getTable<SgofTable>(satRegionIdx));
            }
            else {
                // gas-oil parameters are specified using the SLGOF keyword
                assert(!slgofTables.empty());

                extractUnscaledSlgof_(slgofTables.getTable<SlgofTable>(satRegionIdx));
            }
        }
        else if (family2) {
            extractUnscaledSwfn_(swfnTables.getTable<SwfnTable>(satRegionIdx));
            extractUnscaledSgfn_(sgfnTables.getTable<SgfnTable>(satRegionIdx));
            extractUnscaledSof3_(sof3Tables.getTable<Sof3Table>(satRegionIdx));
        }
        else {
            throw std::domain_error("No valid saturation keyword family specified");
        }

        // there are no "unscaled" Leverett factors, so we just set them to 1.0
        pcowLeverettFactor = 1.0;
        pcgoLeverettFactor = 1.0;
    }

    /*!
     * \brief Extract the values of the scaled scaling parameters.
     *
     * I.e., the values which are "seen" by the physical model.
     */
    void extractScaled(const Ewoms::EclipseState& eclState,
                       const EclEpsGridProperties& epsProperties,
                       unsigned activeIndex)
    {
        // overwrite the unscaled values with the values for the cell if it is
        // explicitly specified by the corresponding keyword.
        extractGridPropertyValue_(Swl, epsProperties.compressedSwl, activeIndex);
        extractGridPropertyValue_(Sgl, epsProperties.compressedSgl, activeIndex);
        extractGridPropertyValue_(Swcr, epsProperties.compressedSwcr, activeIndex);
        extractGridPropertyValue_(Sgcr, epsProperties.compressedSgcr, activeIndex);

        extractGridPropertyValue_(Sowcr, epsProperties.compressedSowcr, activeIndex);
        extractGridPropertyValue_(Sogcr, epsProperties.compressedSogcr, activeIndex);
        extractGridPropertyValue_(Swu, epsProperties.compressedSwu, activeIndex);
        extractGridPropertyValue_(Sgu, epsProperties.compressedSgu, activeIndex);
        extractGridPropertyValue_(maxPcow, epsProperties.compressedPcw, activeIndex);
        extractGridPropertyValue_(maxPcgo, epsProperties.compressedPcg, activeIndex);
        extractGridPropertyValue_(maxKrw, epsProperties.compressedKrw, activeIndex);
        extractGridPropertyValue_(maxKrg, epsProperties.compressedKrg, activeIndex);
        // quite likely that's wrong!
        extractGridPropertyValue_(maxKrow, epsProperties.compressedKro, activeIndex);
        extractGridPropertyValue_(maxKrog, epsProperties.compressedKro, activeIndex);

        // compute the Leverett capillary pressure scaling factors if applicable.  note
        // that this needs to be done using non-SI units to make it correspond to the
        // documentation.
        pcowLeverettFactor = 1.0;
        pcgoLeverettFactor = 1.0;
        if (eclState.getTableManager().useJFunc()) {
            const auto& jfunc = eclState.getTableManager().getJFunc();
            const auto& jfuncDir = jfunc.direction();

            Scalar perm;
            if (jfuncDir == Ewoms::JFunc::Direction::X)
                perm = epsProperties.compressedPermx[activeIndex];
            else if (jfuncDir == Ewoms::JFunc::Direction::Y)
                perm = epsProperties.compressedPermy[activeIndex];
            else if (jfuncDir == Ewoms::JFunc::Direction::Z)
                perm = epsProperties.compressedPermz[activeIndex];
            else if (jfuncDir == Ewoms::JFunc::Direction::XY)
            {
                // TODO: verify that this really is the arithmetic mean. (the
                // documentation just says that the "average" should be used, IMO the
                // harmonic mean would be more appropriate because that's what's usually
                // applied when calculating the fluxes.)
                double permx = epsProperties.compressedPermx[activeIndex];
                double permy = epsProperties.compressedPermy[activeIndex];
                perm = Ewoms::arithmeticMean(permx, permy);
            }
            else
                throw std::runtime_error("Illegal direction indicator for the JFUNC "
                                         "keyword ("+std::to_string(int(jfuncDir))+")");

            // convert permeability from m^2 to mD
            perm *= 1.01325e15;

            Scalar poro = epsProperties.compressedPoro[activeIndex];
            Scalar alpha = jfunc.alphaFactor();
            Scalar beta = jfunc.betaFactor();

            // the part of the Leverett capillary pressure which does not depend on
            // surface tension.
            Scalar commonFactor = std::pow(poro, alpha)/std::pow(perm, beta);

            // multiply the documented constant by 10^5 because we want the pressures
            // in [Pa], not in [bar]
            const Scalar Uconst = 0.318316 * 1e5;

            // compute the oil-water Leverett factor.
            const auto& jfuncFlag = jfunc.flag();
            if (jfuncFlag == Ewoms::JFunc::Flag::WATER || jfuncFlag == Ewoms::JFunc::Flag::BOTH) {
                // note that we use the surface tension in terms of [dyn/cm]
                Scalar gamma =
                    jfunc.owSurfaceTension();
                pcowLeverettFactor = commonFactor*gamma*Uconst;
            }

            // compute the gas-oil Leverett factor.
            if (jfuncFlag == Ewoms::JFunc::Flag::GAS || jfuncFlag == Ewoms::JFunc::Flag::BOTH) {
                // note that we use the surface tension in terms of [dyn/cm]
                Scalar gamma =
                    jfunc.goSurfaceTension();
                pcgoLeverettFactor = commonFactor*gamma*Uconst;
            }
        }
    }
#endif

private:
#if HAVE_ECL_INPUT
    void extractUnscaledSgof_(const Ewoms::SgofTable& sgofTable)
    {
        // minimum gas and oil-in-gas-oil saturation
        Sgl = sgofTable.getSgColumn().front();
        Sogl = 1.0 - sgofTable.getSgColumn().back();

        // maximum gas and oil-in-gas-oil saturation
        Sgu = sgofTable.getSgColumn().back();
        Sogu = 1.0 - sgofTable.getSgColumn().front();

        // critical gas saturation
        Sgcr = 0.0;
        for (size_t rowIdx = 0; rowIdx < sgofTable.numRows(); ++ rowIdx) {
            if (sgofTable.getKrgColumn()[rowIdx] > krCriticalEps)
                break;

            Sgcr = sgofTable.getSgColumn()[rowIdx];
        }

        // critical oil saturation of gas-oil system
        Sogcr = 0.0;
        for (int rowIdx = static_cast<int>(sgofTable.numRows() - 1);
             rowIdx >= 0;
             -- rowIdx)
        {
            if (sgofTable.getKrogColumn()[static_cast<size_t>(rowIdx)] > krCriticalEps)
                break;

            Sogcr = 1.0 - sgofTable.getSgColumn()[static_cast<size_t>(rowIdx)];
        }

        // maximum gas-oil capillary pressure
        maxPcgo = sgofTable.getPcogColumn().back();

        // maximum gas-* relperms
        maxKrg = sgofTable.getKrgColumn().back();
        maxKrog = sgofTable.getKrogColumn().front();
    }

    void extractUnscaledSlgof_(const Ewoms::SlgofTable& slgofTable)
    {
        // minimum gas and oil-in-gas-oil saturation
        Sgl = 1.0 - slgofTable.getSlColumn().back();
        Sogl = slgofTable.getSlColumn().front();

        // maximum gas and oil-in-gas-oil saturation
        Sgu = 1.0 - slgofTable.getSlColumn().front();
        Sogu = slgofTable.getSlColumn().back();

        // critical gas saturation
        Sgcr = 0.0;
        for (int rowIdx = static_cast<int>(slgofTable.numRows()) - 1;
             rowIdx >= 0;
             -- rowIdx)
        {
            if (slgofTable.getKrgColumn()[static_cast<size_t>(rowIdx)] > krCriticalEps)
                break;

            Sgcr = 1 - slgofTable.getSlColumn()[static_cast<size_t>(rowIdx)];
        }

        // critical oil saturation of gas-oil system
        Sogcr = 0.0;
        for (size_t rowIdx = 0; rowIdx < slgofTable.numRows(); ++ rowIdx) {
            if (slgofTable.getKrogColumn()[rowIdx] > krCriticalEps)
                break;

            Sogcr = slgofTable.getSlColumn()[rowIdx];
        }

        // maximum gas-oil capillary pressure
        maxPcgo = slgofTable.getPcogColumn().front();

        // maximum gas-* relperms
        maxKrg = slgofTable.getKrgColumn().front();
        maxKrog = slgofTable.getKrogColumn().back();
    }

    void extractUnscaledSwof_(const Ewoms::SwofTable& swofTable)
    {
        // connate saturations
        Swl = swofTable.getSwColumn().front();
        Sowl = 1.0 - swofTable.getSwColumn().back();

        // maximum water and oil-in-oil-water saturations
        Swu = swofTable.getSwColumn().back();
        Sowu = 1.0 - swofTable.getSwColumn().front();

        // critical water saturation
        Swcr = 0.0;
        for (size_t rowIdx = 0; rowIdx < swofTable.numRows(); ++ rowIdx) {
            if (swofTable.getKrwColumn()[rowIdx] > krCriticalEps)
                break;

            Swcr = swofTable.getSwColumn()[rowIdx];
        }

        // critical oil saturation of oil-water system
        Sowcr = 0.0;
        for (int rowIdx = static_cast<int>(swofTable.numRows()) - 1;
             rowIdx >= 0;
             -- rowIdx)
        {
            if (swofTable.getKrowColumn()[static_cast<size_t>(rowIdx)] > krCriticalEps)
                break;

            Sowcr = 1.0 - swofTable.getSwColumn()[static_cast<size_t>(rowIdx)];
        }

        // maximum oil-water capillary pressures
        maxPcow = swofTable.getPcowColumn().front();

        // maximum water-* relative permeabilities
        maxKrw = swofTable.getKrwColumn().back();
        maxKrow = swofTable.getKrowColumn().front();
    }

    void extractUnscaledSwfn_(const Ewoms::SwfnTable& swfnTable)
    {
        // connate water saturation
        Swl = swfnTable.getSwColumn().front();

        // maximum water saturation
        Swu = swfnTable.getSwColumn().back();

        // critical water saturation
        Swcr = 0.0;
        for (size_t rowIdx = 0; rowIdx < swfnTable.numRows(); ++ rowIdx) {
            if (swfnTable.getKrwColumn()[rowIdx] > krCriticalEps)
                break;

            Swcr = swfnTable.getSwColumn()[rowIdx];
        }

        // maximum oil-water capillary pressure
        maxPcow = swfnTable.getPcowColumn().front();

        // maximum water relative permeability
        maxKrw = swfnTable.getKrwColumn().back();
    }

    void extractUnscaledSgfn_(const Ewoms::SgfnTable& sgfnTable)
    {
        // connate gas saturation
        Sgl = sgfnTable.getSgColumn().front();

        // maximum gas saturations
        Sgu = sgfnTable.getSgColumn().back();
        Sogu = 1 - sgfnTable.getSgColumn().front();

        // critical gas saturation
        Sgcr = 0.0;
        for (size_t rowIdx = 0; rowIdx < sgfnTable.numRows(); ++ rowIdx) {
            if (sgfnTable.getKrgColumn()[rowIdx] > krCriticalEps)
                break;

            Sgcr = sgfnTable.getSgColumn()[rowIdx];
        }

        // maximum capillary pressure
        maxPcgo = sgfnTable.getPcogColumn().back();

        // maximum relative gas permeability
        maxKrg = sgfnTable.getKrgColumn().back();
    }

    void extractUnscaledSof3_(const Ewoms::Sof3Table& sof3Table)
    {
        // connate oil saturations
        Sowl = sof3Table.getSoColumn().front() + Sgl;
        Sogl = sof3Table.getSoColumn().front() + Swl;

        // maximum oil saturations
        Sowu = sof3Table.getSoColumn().back();

        // critical oil saturation of oil-water system
        Sowcr = 0.0;
        for (size_t rowIdx = 0 ; rowIdx < sof3Table.numRows(); ++ rowIdx) {
            if (sof3Table.getKrowColumn()[rowIdx] > krCriticalEps) {
                break;
            }

            Sowcr = sof3Table.getSoColumn()[rowIdx];
        }

        // critical oil saturation of gas-oil system
        Sogcr = 0.0;
        for (size_t rowIdx = 0 ; rowIdx < sof3Table.numRows(); ++ rowIdx) {
            if (sof3Table.getKrogColumn()[rowIdx] > krCriticalEps)
                break;

            Sogcr = sof3Table.getSoColumn()[rowIdx];
        }

        // maximum relative oil permeabilities
        maxKrow = sof3Table.getKrowColumn().back();
        maxKrog = sof3Table.getKrogColumn().back();
    }

    void extractUnscaledSof2_(const Ewoms::Sof2Table& sof2Table)
    {
        // connate oil saturations
        Sowl = sof2Table.getSoColumn().front() + Sgl;
        Sogl = sof2Table.getSoColumn().front() + Swl;

        // maximum oil saturations
        Sowu = sof2Table.getSoColumn().back();

        // critical oil saturation of oil-water system or critical oil saturation of
        // gas-oil system
        Sowcr = 0.0;
        for (size_t rowIdx = 0 ; rowIdx < sof2Table.numRows(); ++ rowIdx) {
            if (sof2Table.getKroColumn()[rowIdx] > krCriticalEps) {
                break;
            }

            Sowcr = sof2Table.getSoColumn()[rowIdx];
        }
        Sogcr = Sowcr;

        // maximum relative oil permeabilities
        maxKrow = sof2Table.getKroColumn().back();
        maxKrog = maxKrow;
    }
#endif // HAVE_ECL_INPUT

    void extractGridPropertyValue_(Scalar& targetValue,
                                   const std::vector<double>& propData,
                                   unsigned activeCellIdx)
    {
        if (propData.empty())
            return;

        targetValue = propData[activeCellIdx];
    }
};

/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Represents the points on the X and Y axis to be scaled if endpoint scaling is
 *        used.
 */
template <class Scalar>
class EclEpsScalingPoints
{
public:
    /*!
     * \brief Assigns the scaling points which actually ought to be used.
     */
    void init(const EclEpsScalingPointsInfo<Scalar>& epsInfo,
              const EclEpsConfig& config,
              EclTwoPhaseSystemType epsSystemType)
    {
        if (epsSystemType == EclOilWaterSystem) {
            // saturation scaling for capillary pressure
            saturationPcPoints_[0] = epsInfo.Swl;
            saturationPcPoints_[1] = epsInfo.Swu;

            // krw saturation scaling endpoints
            if (config.enableThreePointKrSatScaling()) {
                saturationKrwPoints_[0] = epsInfo.Swcr;
                saturationKrwPoints_[1] = 1.0 - epsInfo.Sowcr - epsInfo.Sgl;
                saturationKrwPoints_[2] = epsInfo.Swu;
            }
            else {
                saturationKrwPoints_[0] = epsInfo.Swcr;
                saturationKrwPoints_[1] = epsInfo.Swu;
            }

            // krn saturation scaling endpoints (with the non-wetting phase being oil).
            // because ewoms-material specifies non-wetting phase relperms in terms of the
            // wetting phase saturations, the code here uses 1 minus the values specified
            // by the Eclipse TD and the order of the scaling points is reversed
            if (config.enableThreePointKrSatScaling()) {
                saturationKrnPoints_[2] = 1.0 - epsInfo.Sowcr;
                saturationKrnPoints_[1] = epsInfo.Swcr + epsInfo.Sgl;
                saturationKrnPoints_[0] = epsInfo.Swl + epsInfo.Sgl;
            }
            else {
                saturationKrnPoints_[1] = 1 - epsInfo.Sowcr;
                saturationKrnPoints_[0] = epsInfo.Swl + epsInfo.Sgl;
            }

            if (config.enableLeverettScaling())
                maxPcnwOrLeverettFactor_ = epsInfo.pcowLeverettFactor;
            else
                maxPcnwOrLeverettFactor_ = epsInfo.maxPcow;
            maxKrw_ = epsInfo.maxKrw;
            maxKrn_ = epsInfo.maxKrow;
        }
        else {
            assert(epsSystemType == EclGasOilSystem);

            // saturation scaling for capillary pressure
            saturationPcPoints_[0] = 1.0 - epsInfo.Sgu;
            saturationPcPoints_[1] = 1.0 - epsInfo.Sgl;

            // krw saturation scaling endpoints
            if (config.enableThreePointKrSatScaling()) {
                saturationKrwPoints_[0] = epsInfo.Sogcr;
                saturationKrwPoints_[1] = 1 - epsInfo.Sgcr - epsInfo.Swl;
                saturationKrwPoints_[2] = 1 - epsInfo.Swl - epsInfo.Sgl;
            }
            else {
                saturationKrwPoints_[0] = epsInfo.Sogcr;
                saturationKrwPoints_[1] = 1 - epsInfo.Swl - epsInfo.Sgl;
            }

            // krn saturation scaling endpoints (with the non-wetting phase being gas).
            // because ewoms-material specifies non-wetting phase relperms in terms of the
            // wetting phase saturations, the code here uses 1 minus the values specified
            // by the Eclipse TD and the order of the scaling points is reversed
            if (config.enableThreePointKrSatScaling()) {
                saturationKrnPoints_[2] = 1.0 - epsInfo.Sgcr;
                saturationKrnPoints_[1] = epsInfo.Sogcr + epsInfo.Swl;
                saturationKrnPoints_[0] = 1.0 - epsInfo.Sgu;
            }
            else {
                saturationKrnPoints_[1] = 1.0 - epsInfo.Sgcr;
                saturationKrnPoints_[0] = 1.0 - epsInfo.Sgu;
            }

            if (config.enableLeverettScaling())
                maxPcnwOrLeverettFactor_ = epsInfo.pcgoLeverettFactor;
            else
                maxPcnwOrLeverettFactor_ = epsInfo.maxPcgo;

            maxKrw_ = epsInfo.maxKrog;
            maxKrn_ = epsInfo.maxKrg;
        }
    }

    /*!
     * \brief Sets an saturation value for capillary pressure saturation scaling
     */
    void setSaturationPcPoint(unsigned pointIdx, Scalar value)
    { saturationPcPoints_[pointIdx] = value; }

    /*!
     * \brief Returns the points used for capillary pressure saturation scaling
     */
    const std::array<Scalar, 2>& saturationPcPoints() const
    { return saturationPcPoints_; }

    /*!
     * \brief Sets an saturation value for wetting-phase relperm saturation scaling
     */
    void setSaturationKrwPoint(unsigned pointIdx, Scalar value)
    { saturationKrwPoints_[pointIdx] = value; }

    /*!
     * \brief Returns the points used for wetting phase relperm saturation scaling
     */
    const std::array<Scalar, 3>& saturationKrwPoints() const
    { return saturationKrwPoints_; }

    /*!
     * \brief Sets an saturation value for non-wetting phase relperm saturation scaling
     */
    void setSaturationKrnPoint(unsigned pointIdx, Scalar value)
    { saturationKrnPoints_[pointIdx] = value; }

    /*!
     * \brief Returns the points used for non-wetting phase relperm saturation scaling
     */
    const std::array<Scalar, 3>& saturationKrnPoints() const
    { return saturationKrnPoints_; }

    /*!
     * \brief Sets the maximum capillary pressure
     */
    void setMaxPcnw(Scalar value)
    { maxPcnwOrLeverettFactor_ = value; }

    /*!
     * \brief Returns the maximum capillary pressure
     */
    Scalar maxPcnw() const
    { return maxPcnwOrLeverettFactor_; }

    /*!
     * \brief Sets the Leverett scaling factor for capillary pressure
     */
    void setLeverettFactor(Scalar value)
    { maxPcnwOrLeverettFactor_ = value; }

    /*!
     * \brief Returns the Leverett scaling factor for capillary pressure
     */
    Scalar leverettFactor() const
    { return maxPcnwOrLeverettFactor_; }

    /*!
     * \brief Sets the maximum wetting phase relative permeability
     */
    void setMaxKrw(Scalar value)
    { maxKrw_ = value; }

    /*!
     * \brief Returns the maximum wetting phase relative permeability
     */
    Scalar maxKrw() const
    { return maxKrw_; }

    /*!
     * \brief Sets the maximum wetting phase relative permeability
     */
    void setMaxKrn(Scalar value)
    { maxKrn_ = value; }

    /*!
     * \brief Returns the maximum wetting phase relative permeability
     */
    Scalar maxKrn() const
    { return maxKrn_; }

    void print() const
    {
        std::cout << "    saturationKrnPoints_[0]: " << saturationKrnPoints_[0] << "\n"
                  << "    saturationKrnPoints_[1]: " << saturationKrnPoints_[1] << "\n"
                  << "    saturationKrnPoints_[2]: " << saturationKrnPoints_[2] << "\n";
    }

private:
    // The the points used for the "y-axis" scaling of capillary pressure
    Scalar maxPcnwOrLeverettFactor_;

    // The the points used for the "y-axis" scaling of wetting phase relative permability
    Scalar maxKrw_;

    // The the points used for the "y-axis" scaling of non-wetting phase relative permability
    Scalar maxKrn_;

    // The the points used for saturation ("x-axis") scaling of capillary pressure
    std::array<Scalar, 2> saturationPcPoints_;

    // The the points used for saturation ("x-axis") scaling of wetting phase relative permeability
    std::array<Scalar, 3> saturationKrwPoints_;

    // The the points used for saturation ("x-axis") scaling of non-wetting phase relative permeability
    std::array<Scalar, 3> saturationKrnPoints_;
};

} // namespace Ewoms

#endif
