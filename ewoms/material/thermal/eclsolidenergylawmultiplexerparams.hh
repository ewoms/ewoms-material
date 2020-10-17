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
 * \copydoc Ewoms::EclSolidEnergyLawMultiplexerParams
 */
#ifndef EWOMS_ECL_SOLID_ENERGY_LAW_MULTIPLEXER_PARAMS_HH
#define EWOMS_ECL_SOLID_ENERGY_LAW_MULTIPLEXER_PARAMS_HH

#include "eclheatcrlawparams.hh"
#include "eclspecrocklawparams.hh"

#include <ewoms/material/common/ensurefinalized.hh>

#include <memory>

namespace Ewoms {

/*!
 * \brief The default implementation of a parameter object for the
 *        ECL thermal law.
 */
template <class ScalarT>
class EclSolidEnergyLawMultiplexerParams : public EnsureFinalized
{
    typedef void* ParamPointerType;

public:
    typedef ScalarT Scalar;

    enum SolidEnergyApproach {
        undefinedApproach,
        heatcrApproach, // keywords: HEATCR, HEATCRT, STCOND
        specrockApproach, // keyword: SPECROCK
        nullApproach, // (no keywords)
    };

    typedef Ewoms::EclHeatcrLawParams<ScalarT> HeatcrLawParams;
    typedef Ewoms::EclSpecrockLawParams<ScalarT> SpecrockLawParams;

    EclSolidEnergyLawMultiplexerParams(const EclSolidEnergyLawMultiplexerParams&) = default;

    EclSolidEnergyLawMultiplexerParams()
    { solidEnergyApproach_ = undefinedApproach; }

    ~EclSolidEnergyLawMultiplexerParams()
    { destroy_(); }

    void setSolidEnergyApproach(SolidEnergyApproach newApproach)
    {
        destroy_();

        solidEnergyApproach_ = newApproach;
        switch (solidEnergyApproach()) {
        case undefinedApproach:
            throw std::logic_error("Cannot set the approach for solid energy storage to 'undefined'!");

        case heatcrApproach:
            realParams_ = new HeatcrLawParams;
            break;

        case specrockApproach:
            realParams_ = new SpecrockLawParams;
            break;

        case nullApproach:
            realParams_ = nullptr;
            break;
        }
    }

    SolidEnergyApproach solidEnergyApproach() const
    { return solidEnergyApproach_; }

    // get the parameter object for the HEATCR case
    template <SolidEnergyApproach approachV>
    typename std::enable_if<approachV == heatcrApproach, HeatcrLawParams>::type&
    getRealParams()
    {
        assert(solidEnergyApproach() == approachV);
        return *static_cast<HeatcrLawParams*>(realParams_);
    }

    template <SolidEnergyApproach approachV>
    typename std::enable_if<approachV == heatcrApproach, const HeatcrLawParams>::type&
    getRealParams() const
    {
        assert(solidEnergyApproach() == approachV);
        return *static_cast<const HeatcrLawParams*>(realParams_);
    }

    // get the parameter object for the SPECROCK case
    template <SolidEnergyApproach approachV>
    typename std::enable_if<approachV == specrockApproach, SpecrockLawParams>::type&
    getRealParams()
    {
        assert(solidEnergyApproach() == approachV);
        return *static_cast<SpecrockLawParams*>(realParams_);
    }

    template <SolidEnergyApproach approachV>
    typename std::enable_if<approachV == specrockApproach, const SpecrockLawParams>::type&
    getRealParams() const
    {
        assert(solidEnergyApproach() == approachV);
        return *static_cast<const SpecrockLawParams*>(realParams_);
    }

private:
    void destroy_()
    {
        switch (solidEnergyApproach()) {
        case undefinedApproach:
            break;

        case heatcrApproach:
            delete static_cast<HeatcrLawParams*>(realParams_);
            break;

        case specrockApproach:
            delete static_cast<SpecrockLawParams*>(realParams_);
            break;

        case nullApproach:
            break;
        }

        solidEnergyApproach_ = undefinedApproach;
    }

    SolidEnergyApproach solidEnergyApproach_;
    ParamPointerType realParams_;
};

} // namespace Ewoms

#endif
