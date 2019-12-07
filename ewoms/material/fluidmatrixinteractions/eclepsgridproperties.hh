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
#include <ewoms/eclio/parser/eclipsestate/grid/gridproperty.hh>
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
    EclEpsGridProperties(const Ewoms::EclipseState& eclState,
                         bool useImbibition,
                         const std::vector<int>& compressedToCartesianElemIdx)
    {
        std::string kwPrefix = useImbibition?"I":"";

        const auto& ecl3dProps = eclState.get3DProperties();

        if (useImbibition)
            compressedSatnum = std::move(createCompressedPropertyData_(ecl3dProps.getIntGridProperty("IMBNUM").getData(), compressedToCartesianElemIdx));
        else
            compressedSatnum = std::move(createCompressedPropertyData_(ecl3dProps.getIntGridProperty("SATNUM").getData(), compressedToCartesianElemIdx));

        compressedSwl = std::move(retrieveCompressedPropertyData_(ecl3dProps, kwPrefix+"SWL", compressedToCartesianElemIdx));
        compressedSgl = std::move(retrieveCompressedPropertyData_(ecl3dProps, kwPrefix+"SGL", compressedToCartesianElemIdx));
        compressedSwcr = std::move(retrieveCompressedPropertyData_(ecl3dProps, kwPrefix+"SWCR", compressedToCartesianElemIdx));
        compressedSgcr = std::move(retrieveCompressedPropertyData_(ecl3dProps, kwPrefix+"SGCR", compressedToCartesianElemIdx));
        compressedSowcr = std::move(retrieveCompressedPropertyData_(ecl3dProps, kwPrefix+"SOWCR", compressedToCartesianElemIdx));
        compressedSogcr = std::move(retrieveCompressedPropertyData_(ecl3dProps, kwPrefix+"SOGCR", compressedToCartesianElemIdx));
        compressedSwu = std::move(retrieveCompressedPropertyData_(ecl3dProps, kwPrefix+"SWU", compressedToCartesianElemIdx));
        compressedSgu = std::move(retrieveCompressedPropertyData_(ecl3dProps, kwPrefix+"SGU", compressedToCartesianElemIdx));
        compressedPcw = std::move(retrieveCompressedPropertyData_(ecl3dProps, kwPrefix+"PCW", compressedToCartesianElemIdx));
        compressedPcg = std::move(retrieveCompressedPropertyData_(ecl3dProps, kwPrefix+"PCG", compressedToCartesianElemIdx));
        compressedKrw = std::move(retrieveCompressedPropertyData_(ecl3dProps, kwPrefix+"KRW", compressedToCartesianElemIdx));
        compressedKro = std::move(retrieveCompressedPropertyData_(ecl3dProps, kwPrefix+"KRO", compressedToCartesianElemIdx));
        compressedKrg = std::move(retrieveCompressedPropertyData_(ecl3dProps, kwPrefix+"KRG", compressedToCartesianElemIdx));

        // _may_ be needed to calculate the Leverett capillary pressure scaling factor
        compressedPoro = std::move(retrieveCompressedPropertyData_(ecl3dProps, "PORO", compressedToCartesianElemIdx));

        compressedPermx = std::move(createCompressedPropertyData_(ecl3dProps.getDoubleGridProperty("PERMX").getData(), compressedToCartesianElemIdx));

        if (ecl3dProps.hasDeckDoubleGridProperty("PERMY"))
            compressedPermy = std::move(createCompressedPropertyData_(ecl3dProps.getDoubleGridProperty("PERMY").getData(), compressedToCartesianElemIdx));
        else
            compressedPermy = compressedPermx;

        if (ecl3dProps.hasDeckDoubleGridProperty("PERMZ"))
            compressedPermz = std::move(createCompressedPropertyData_(ecl3dProps.getDoubleGridProperty("PERMZ").getData(), compressedToCartesianElemIdx));
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

private:
    // given the data of a ECL grid property in a logically Cartesian grid, return the
    // array for corresponding "compressed" field (i.e., for the grid where inactive
    // cells are removed).
    template <typename FieldType>
    std::vector<FieldType>
    createCompressedPropertyData_(const std::vector<FieldType>& cartesianPropertyData,
                                  const std::vector<int>& compressedToCartesianElemIdx)
    {
        std::vector<FieldType> compressedPropertyData(compressedToCartesianElemIdx.size());

        for (size_t activeIdx = 0; activeIdx < compressedToCartesianElemIdx.size(); activeIdx++) {
            auto cartesianIdx = compressedToCartesianElemIdx[activeIdx];
            compressedPropertyData[activeIdx] = cartesianPropertyData[cartesianIdx];
        }

        return compressedPropertyData;
    }

    // check if a given ECL grid property exists and if it does, create and return a
    // compressed copy of it
    std::vector<double>
    retrieveCompressedPropertyData_(const Eclipse3DProperties& props,
                                    const std::string& keyword,
                                    const std::vector<int>& compressedToCartesianElemIdx)
    {
        if (props.hasDeckDoubleGridProperty(keyword))
            return std::move(createCompressedPropertyData_(props.getDoubleGridProperty(keyword).getData(),
                                                           compressedToCartesianElemIdx));

        return std::move(std::vector<double>());
    }
};

} // namespace Ewoms

#endif

