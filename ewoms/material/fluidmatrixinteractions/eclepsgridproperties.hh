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
 * \copydoc Ewoms::EclEpsTwoPhaseLawPoints
 */
#ifndef EWOMS_ECL_EPS_GRID_PROPERTIES_HH
#define EWOMS_ECL_EPS_GRID_PROPERTIES_HH

#include "eclepsconfig.hh"

#include <ewoms/common/means.hh>

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

#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <memory>

namespace Ewoms {
/*!
 * \brief Collects all grid properties which are relevant for end point scaling.
 *
 * This class is used for both, the drainage and the imbibition variants of the ECL
 * keywords.
 */
class EclEpsGridProperties
{
public:
    // extract the field properties for required to initialize the ECL-style fluid-matrix
    // interactions. note that for this, the EclipseState's fieldProps must exhibit the
    // correct set of active set (confer FieldPropsManager::reset_actnum())
    EclEpsGridProperties(const Ewoms::EclipseState& eclState,
                         bool useImbibition,
                         const std::vector<int>& compressedToCartesianElemIdx)
    {
        std::string kwPrefix = useImbibition?"I":"";

        const auto& fieldProps = eclState.fieldProps();

        if (useImbibition)
            compressedSatnum = fieldProps.get_copy<int>("IMBNUM");
        else
            compressedSatnum = fieldProps.get_copy<int>("SATNUM");

        compressedSwl = fieldProps.get_copy<double>(kwPrefix+"SWL");
        compressedSgl = fieldProps.get_copy<double>(kwPrefix+"SGL");
        compressedSwcr = fieldProps.get_copy<double>(kwPrefix+"SWCR");
        compressedSgcr = fieldProps.get_copy<double>(kwPrefix+"SGCR");
        compressedSowcr = fieldProps.get_copy<double>(kwPrefix+"SOWCR");
        compressedSogcr = fieldProps.get_copy<double>(kwPrefix+"SOGCR");
        compressedSwu = fieldProps.get_copy<double>(kwPrefix+"SWU");
        compressedSgu = fieldProps.get_copy<double>(kwPrefix+"SGU");
        compressedPcw = fieldProps.get_copy<double>(kwPrefix+"PCW");
        compressedPcg = fieldProps.get_copy<double>(kwPrefix+"PCG");
        compressedKrw = fieldProps.get_copy<double>(kwPrefix+"KRW");
        compressedKro = fieldProps.get_copy<double>(kwPrefix+"KRO");
        compressedKrg = fieldProps.get_copy<double>(kwPrefix+"KRG");

        // _may_ be needed to calculate the Leverett capillary pressure scaling factor
        compressedPoro = fieldProps.get_copy<double>("PORO");

        if (fieldProps.has_double("PERMX"))
            compressedPermx = fieldProps.get_copy<double>("PERMX");
        else
            compressedPermx.clear();

        if (fieldProps.has_double("PERMY"))
            compressedPermy = fieldProps.get_copy<double>("PERMY");
        else
            compressedPermy = compressedPermx;

        if (fieldProps.has_double("PERMZ"))
            compressedPermz = fieldProps.get_copy<double>("PERMZ");
        else
            compressedPermz = compressedPermx;
    }

    std::vector<int> compressedSatnum;
    std::vector<double> compressedSwl;
    std::vector<double> compressedSgl;
    std::vector<double> compressedSwcr;
    std::vector<double> compressedSgcr;
    std::vector<double> compressedSowcr;
    std::vector<double> compressedSogcr;
    std::vector<double> compressedSwu;
    std::vector<double> compressedSgu;
    std::vector<double> compressedPcw;
    std::vector<double> compressedPcg;
    std::vector<double> compressedKrw;
    std::vector<double> compressedKro;
    std::vector<double> compressedKrg;

    std::vector<double> compressedPermx;
    std::vector<double> compressedPermy;
    std::vector<double> compressedPermz;
    std::vector<double> compressedPoro;
};

} // namespace Ewoms

#endif

