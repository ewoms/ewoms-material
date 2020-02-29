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
 * \copydoc Ewoms::EclThermalLawManager
 */
#if ! HAVE_ECL_INPUT
#error "Eclipse input support in ewoms-eclio is required to use the ECL thermal law manager!"
#endif

#ifndef EWOMS_ECL_THERMAL_LAW_MANAGER_HH
#define EWOMS_ECL_THERMAL_LAW_MANAGER_HH

#include "eclsolidenergylawmultiplexer.hh"
#include "eclsolidenergylawmultiplexerparams.hh"

#include "eclthermalconductionlawmultiplexer.hh"
#include "eclthermalconductionlawmultiplexerparams.hh"

#include <ewoms/common/exceptions.hh>

#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tablemanager.hh>
#include <ewoms/eclio/parser/deck/deck.hh>

namespace Ewoms {

/*!
 * \ingroup fluidmatrixinteractions
 *
 * \brief Provides an simple way to create and manage the thermal law objects
 *        for a complete ECL deck.
 */
template <class Scalar, class FluidSystem>
class EclThermalLawManager
{
public:
    typedef EclSolidEnergyLawMultiplexer<Scalar, FluidSystem> SolidEnergyLaw;
    typedef typename SolidEnergyLaw::Params SolidEnergyLawParams;
    typedef typename SolidEnergyLawParams::HeatcrLawParams HeatcrLawParams;
    typedef typename SolidEnergyLawParams::SpecrockLawParams SpecrockLawParams;

    typedef EclThermalConductionLawMultiplexer<Scalar, FluidSystem> ThermalConductionLaw;
    typedef typename ThermalConductionLaw::Params ThermalConductionLawParams;

    EclThermalLawManager()
    {
        solidEnergyApproach_ = SolidEnergyLawParams::undefinedApproach;
        thermalConductivityApproach_ = ThermalConductionLawParams::undefinedApproach;
    }

    void initParamsForElements(const Ewoms::EclipseState& eclState,
                               const std::vector<int>& compressedToCartesianElemIdx)
    {
        const auto& fieldProps = eclState.fieldProps();
        const auto& tableManager = eclState.getTableManager();
        bool hasHeatcr = fieldProps.has_double("HEATCR");
        bool hasThconr = fieldProps.has_double("THCONR");
        bool hasThc =
            fieldProps.has_double("THCROCK")
            || fieldProps.has_double("THCOIL")
            || fieldProps.has_double("THCGAS")
            || fieldProps.has_double("THCWATER");

        if (hasHeatcr)
            initHeatcr_(eclState, compressedToCartesianElemIdx);
        else if (tableManager.hasTables("SPECROCK"))
            initSpecrock_(eclState, compressedToCartesianElemIdx);
        else
            initNullRockEnergy_();

        if (hasThconr)
            initThconr_(eclState, compressedToCartesianElemIdx);
        else if (hasThc)
            initThc_(eclState, compressedToCartesianElemIdx);
        else
            initNullCond_();
    }

    const SolidEnergyLawParams& solidEnergyLawParams(unsigned elemIdx) const
    {
        switch (solidEnergyApproach_) {
        case SolidEnergyLawParams::heatcrApproach:
            assert(0 <= elemIdx && elemIdx <  solidEnergyLawParams_.size());
            return solidEnergyLawParams_[elemIdx];

        case SolidEnergyLawParams::specrockApproach:
        {
            assert(0 <= elemIdx && elemIdx <  elemToSatnumIdx_.size());
            unsigned satnumIdx = elemToSatnumIdx_[elemIdx];
            assert(0 <= satnumIdx && satnumIdx <  solidEnergyLawParams_.size());
            return solidEnergyLawParams_[satnumIdx];
        }

        case SolidEnergyLawParams::nullApproach:
            return solidEnergyLawParams_[0];

        default:
            throw std::runtime_error("Attempting to retrieve solid energy storage parameters "
                                     "without a known approach being defined by the deck.");
        }
    }

    const ThermalConductionLawParams& thermalConductionLawParams(unsigned elemIdx) const
    {
        switch (thermalConductivityApproach_) {
        case ThermalConductionLawParams::thconrApproach:
        case ThermalConductionLawParams::thcApproach:
            assert(0 <= elemIdx && elemIdx <  thermalConductionLawParams_.size());
            return thermalConductionLawParams_[elemIdx];

        case ThermalConductionLawParams::nullApproach:
            return thermalConductionLawParams_[0];

        default:
            throw std::runtime_error("Attempting to retrieve thermal conduction parameters without "
                                     "a known approach being defined by the deck.");
        }
    }

private:
    /*!
     * \brief Initialize the parameters for the solid energy law using using HEATCR and friends.
     */
    void initHeatcr_(const Ewoms::EclipseState& eclState,
                     const std::vector<int>& compressedToCartesianElemIdx)
    {
        solidEnergyApproach_ = SolidEnergyLawParams::heatcrApproach;
        // actually the value of the reference temperature does not matter for energy
        // conservation. We set it anyway to faciliate comparisons with ECL
        HeatcrLawParams::setReferenceTemperature(FluidSystem::surfaceTemperature);

        const auto& fieldProps = eclState.fieldProps();
        const std::vector<double>& heatcrData  = fieldProps.get_global_double("HEATCR");
        const std::vector<double>& heatcrtData = fieldProps.get_global_double("HEATCRT");
        unsigned numElems = compressedToCartesianElemIdx.size();
        solidEnergyLawParams_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemParam = solidEnergyLawParams_[elemIdx];
            elemParam.setSolidEnergyApproach(SolidEnergyLawParams::heatcrApproach);
            auto& heatcrElemParams = elemParam.template getRealParams<SolidEnergyLawParams::heatcrApproach>();
            unsigned cartesianElemIdx = compressedToCartesianElemIdx[elemIdx];

            heatcrElemParams.setReferenceRockHeatCapacity(heatcrData[cartesianElemIdx]);
            heatcrElemParams.setDRockHeatCapacity_dT(heatcrtData[cartesianElemIdx]);
            heatcrElemParams.finalize();
            elemParam.finalize();
        }
    }

    /*!
     * \brief Initialize the parameters for the solid energy law using using SPECROCK and friends.
     */
    void initSpecrock_(const Ewoms::EclipseState& eclState,
                       const std::vector<int>& compressedToCartesianElemIdx)
    {
        solidEnergyApproach_ = SolidEnergyLawParams::specrockApproach;

        // initialize the element index -> SATNUM index mapping
        const auto& fieldProps = eclState.fieldProps();
        const std::vector<int>& satnumData = fieldProps.get_global_int("SATNUM");
        elemToSatnumIdx_.resize(compressedToCartesianElemIdx.size());
        for (unsigned elemIdx = 0; elemIdx < compressedToCartesianElemIdx.size(); ++ elemIdx) {
            unsigned cartesianElemIdx = compressedToCartesianElemIdx[elemIdx];

            // satnumData contains Fortran-style indices, i.e., they start with 1 instead
            // of 0!
            elemToSatnumIdx_[elemIdx] = satnumData[cartesianElemIdx] - 1;
        }
        // internalize the SPECROCK table
        unsigned numSatRegions = eclState.runspec().tabdims().getNumSatTables();
        const auto& tableManager = eclState.getTableManager();
        solidEnergyLawParams_.resize(numSatRegions);
        for (unsigned satnumIdx = 0; satnumIdx < numSatRegions; ++satnumIdx) {
            const auto& specrockTable = tableManager.getSpecrockTables()[satnumIdx];

            auto& multiplexerParams = solidEnergyLawParams_[satnumIdx];

            multiplexerParams.setSolidEnergyApproach(SolidEnergyLawParams::specrockApproach);

            auto& specrockParams = multiplexerParams.template getRealParams<SolidEnergyLawParams::specrockApproach>();
            const auto& temperatureColumn = specrockTable.getColumn("TEMPERATURE");
            const auto& cvRockColumn = specrockTable.getColumn("CV_ROCK");
            specrockParams.setHeatCapacities(temperatureColumn, cvRockColumn);
            specrockParams.finalize();

            multiplexerParams.finalize();
        }
    }

    /*!
     * \brief Specify the solid energy law by setting heat capacity of rock to 0
     */
    void initNullRockEnergy_()
    {
        solidEnergyApproach_ = SolidEnergyLawParams::nullApproach;

        solidEnergyLawParams_.resize(1);
        solidEnergyLawParams_[0].finalize();
    }

    /*!
     * \brief Initialize the parameters for the thermal conduction law using THCONR and friends.
     */
    void initThconr_(const Ewoms::EclipseState& eclState,
                     const std::vector<int>& compressedToCartesianElemIdx)
    {
        thermalConductivityApproach_ = ThermalConductionLawParams::thconrApproach;

        const auto& fieldProps = eclState.fieldProps();
        std::vector<double> thconrData;
        std::vector<double> thconsfData;
        if (fieldProps.has_double("THCONR"))
            thconrData  = fieldProps.get_global_double("THCONR");

        if (fieldProps.has_double("THCONSF"))
            thconsfData = fieldProps.get_global_double("THCONSF");

        unsigned numElems = compressedToCartesianElemIdx.size();
        thermalConductionLawParams_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemParams = thermalConductionLawParams_[elemIdx];
            elemParams.setThermalConductionApproach(ThermalConductionLawParams::thconrApproach);
            auto& thconrElemParams = elemParams.template getRealParams<ThermalConductionLawParams::thconrApproach>();

            int cartElemIdx = compressedToCartesianElemIdx[elemIdx];
            double thconr = thconrData.empty()   ? 0.0 : thconrData[cartElemIdx];
            double thconsf = thconsfData.empty() ? 0.0 : thconsfData[cartElemIdx];
            thconrElemParams.setReferenceTotalThermalConductivity(thconr);
            thconrElemParams.setDTotalThermalConductivity_dSg(thconsf);

            thconrElemParams.finalize();
            elemParams.finalize();
        }
    }

    /*!
     * \brief Initialize the parameters for the thermal conduction law using THCROCK and friends.
     */
    void initThc_(const Ewoms::EclipseState& eclState,
                  const std::vector<int>& compressedToCartesianElemIdx)
    {
        thermalConductivityApproach_ = ThermalConductionLawParams::thcApproach;

        const auto& fieldProps = eclState.fieldProps();
        std::vector<double> thcrockData;
        std::vector<double> thcoilData;
        std::vector<double> thcgasData;
        std::vector<double> thcwaterData = fieldProps.get_global_double("THCWATER");

        if (fieldProps.has_double("THCROCK"))
            thcrockData = fieldProps.get_global_double("THCROCK");

        if (fieldProps.has_double("THCOIL"))
            thcoilData = fieldProps.get_global_double("THCOIL");

        if (fieldProps.has_double("THCGAS"))
            thcgasData = fieldProps.get_global_double("THCGAS");

        if (fieldProps.has_double("THCWATER"))
            thcwaterData = fieldProps.get_global_double("THCWATER");

        const std::vector<double>& poroData = fieldProps.get_global_double("PORO");

        unsigned numElems = compressedToCartesianElemIdx.size();
        thermalConductionLawParams_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemParams = thermalConductionLawParams_[elemIdx];
            elemParams.setThermalConductionApproach(ThermalConductionLawParams::thcApproach);
            auto& thcElemParams = elemParams.template getRealParams<ThermalConductionLawParams::thcApproach>();

            int cartElemIdx = compressedToCartesianElemIdx[elemIdx];
            thcElemParams.setPorosity(poroData[cartElemIdx]);
            double thcrock = thcrockData.empty()    ? 0.0 : thcrockData[cartElemIdx];
            double thcoil = thcoilData.empty()      ? 0.0 : thcoilData[cartElemIdx];
            double thcgas = thcgasData.empty()      ? 0.0 : thcgasData[cartElemIdx];
            double thcwater = thcwaterData.empty()  ? 0.0 : thcwaterData[cartElemIdx];
            thcElemParams.setThcrock(thcrock);
            thcElemParams.setThcoil(thcoil);
            thcElemParams.setThcgas(thcgas);
            thcElemParams.setThcwater(thcwater);

            thcElemParams.finalize();
            elemParams.finalize();
        }
    }

    /*!
     * \brief Disable thermal conductivity
     */
    void initNullCond_()
    {
        thermalConductivityApproach_ = ThermalConductionLawParams::nullApproach;

        thermalConductionLawParams_.resize(1);
        thermalConductionLawParams_[0].finalize();
    }

private:
    typename ThermalConductionLawParams::ThermalConductionApproach thermalConductivityApproach_;
    typename SolidEnergyLawParams::SolidEnergyApproach solidEnergyApproach_;

    std::vector<unsigned> elemToSatnumIdx_;

    std::vector<SolidEnergyLawParams> solidEnergyLawParams_;
    std::vector<ThermalConductionLawParams> thermalConductionLawParams_;
};
} // namespace Ewoms

#endif
