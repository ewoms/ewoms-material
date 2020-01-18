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
 * \copydoc Ewoms::EclEpsConfig
 */
#ifndef EWOMS_ECL_EPS_CONFIG_HH
#define EWOMS_ECL_EPS_CONFIG_HH

#if HAVE_ECL_INPUT
#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/deck/deckitem.hh>
#include <ewoms/eclio/parser/deck/deckkeyword.hh>
#include <ewoms/eclio/parser/deck/deckrecord.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#endif

#include <ewoms/common/exceptions.hh>
#include <ewoms/common/unused.hh>

#include <string>
#include <cassert>
#include <algorithm>

namespace Ewoms {
/*!
 * \brief Specified which fluids are involved in a given twophase material law for
 *        endpoint scaling.
 */
enum EclTwoPhaseSystemType {
    EclGasOilSystem,
    EclOilWaterSystem,
};

/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Specifies the configuration used by the endpoint scaling code
 *
 * This means which quantities are scaled and how this is to be done.
 */
class EclEpsConfig
{
public:
    EclEpsConfig()
    {
        // by default, do not scale anything
        enableSatScaling_ = false;
        enablePcScaling_ = false;
        enableLeverettScaling_ = false;
        enableKrwScaling_ = false;
        enableKrnScaling_ = false;
        enableThreePointKrSatScaling_ = false;
    }

    /*!
     * \brief Specify whether saturation scaling is enabled.
     */
    void setEnableSatScaling(bool yesno)
    { enableSatScaling_ = yesno; }

    /*!
     * \brief Returns whether saturation scaling is enabled.
     */
    bool enableSatScaling() const
    { return enableSatScaling_; }

    /*!
     * \brief Specify whether three point saturation scaling is enabled for the relative
     *        permeabilities.
     */
    void setEnableThreePointKrSatScaling(bool yesno)
    { enableThreePointKrSatScaling_ = yesno; }

    /*!
     * \brief Returns whether three point saturation scaling is enabled for the relative
     *        permeabilities.
     */
    bool enableThreePointKrSatScaling() const
    { return enableThreePointKrSatScaling_; }

    /*!
     * \brief Specify whether relative permeability scaling is enabled for the wetting phase.
     */
    void setEnableKrwScaling(bool yesno)
    { enableKrwScaling_ = yesno; }

    /*!
     * \brief Returns whether relative permeability scaling is enabled for the wetting phase.
     */
    bool enableKrwScaling() const
    { return enableKrwScaling_; }

    /*!
     * \brief Specify whether relative permeability scaling is enabled for the non-wetting phase.
     */
    void setEnableKrnScaling(bool yesno)
    { enableKrnScaling_ = yesno; }

    /*!
     * \brief Returns whether relative permeability scaling is enabled for the non-wetting phase.
     */
    bool enableKrnScaling() const
    { return enableKrnScaling_; }

    /*!
     * \brief Specify whether capillary pressure scaling is enabled.
     */
    void setEnablePcScaling(bool yesno)
    { enablePcScaling_ = yesno; }

    /*!
     * \brief Returns whether capillary pressure scaling is enabled.
     */
    bool enablePcScaling() const
    { return enablePcScaling_; }

    /*!
     * \brief Specify whether the Leverett capillary pressure scaling is enabled.
     *
     * If setting this to true, Leverett capillary pressure scaling will be used instead
     * of the normal capillary pressure scaling and the value of enablePcScaling() will
     * not matter anymore.
     */
    void setEnableLeverettScaling(bool yesno)
    { enableLeverettScaling_ = yesno; }

    /*!
     * \brief Returns whether the Leverett capillary pressure scaling is enabled.
     *
     * If this returns true, Leverett capillary pressure scaling will be used instead of
     * the normal capillary pressure scaling and the value of enablePcScaling() does not
     * matter anymore.
     */
    bool enableLeverettScaling() const
    { return enableLeverettScaling_; }

    /*!
     * \brief Reads all relevant material parameters form a cell of a parsed ECL deck.
     *
     * This requires that the ewoms-eclio module is available.
     */
    void initFromDeck(const Ewoms::Deck& deck EWOMS_UNUSED,
                      const Ewoms::EclipseState& eclState,
                      Ewoms::EclTwoPhaseSystemType twoPhaseSystemType)
    {
        const auto& endscale = eclState.runspec().endpointScaling();
        // find out if endpoint scaling is used in the first place
        if (!endscale) {
            // it is not used, i.e., just set all enable$Foo attributes to 0 and be done
            // with it.
            enableSatScaling_ = false;
            enableThreePointKrSatScaling_ = false;
            enablePcScaling_ = false;
            enableLeverettScaling_ = false;
            enableKrwScaling_ = false;
            enableKrnScaling_ = false;
            return;
        }
        // endpoint scaling is used, i.e., at least saturation scaling needs to be enabled
        enableSatScaling_ = true;
        enableThreePointKrSatScaling_ = endscale.threepoint();

        if (twoPhaseSystemType == EclOilWaterSystem) {
            // check if Leverett capillary pressure scaling is requested
            if (eclState.getTableManager().useJFunc()) {
                const auto& jfunc = eclState.getTableManager().getJFunc();
                auto flag = jfunc.flag();
                if (flag == Ewoms::JFunc::Flag::BOTH || flag == Ewoms::JFunc::Flag::WATER)
                    enableLeverettScaling_ = true;
            }
        } else {
            // check if Leverett capillary pressure scaling is requested
            if (eclState.getTableManager().useJFunc()) {
                const auto& jfunc = eclState.getTableManager().getJFunc();
                auto flag = jfunc.flag();
                if (flag == Ewoms::JFunc::Flag::BOTH || flag == Ewoms::JFunc::Flag::GAS)
                    enableLeverettScaling_ = true;
            }
        }

        const auto& fieldProps = eclState.fieldProps();
        // check if we are supposed to scale the Y axis of the capillary pressure
        if (twoPhaseSystemType == EclOilWaterSystem) {
            enableKrnScaling_ = fieldProps.has<double>("KRO");
            enableKrwScaling_ = fieldProps.has<double>("KRW");
            enablePcScaling_  = fieldProps.has<double>("PCW") || fieldProps.has<double>("SWATINIT");
        }
        else {
            assert(twoPhaseSystemType == EclGasOilSystem);
            enableKrnScaling_ = fieldProps.has<double>("KRG");
            enableKrwScaling_ = fieldProps.has<double>("KRO");
            enablePcScaling_  = fieldProps.has<double>("PCG");
        }

        if (enablePcScaling_ && enableLeverettScaling_)
            throw std::runtime_error("Capillary pressure scaling and the Leverett scaling function are "
                                     "mutually exclusive: The deck contains the PCW/PCG property and the "
                                     "JFUNC keyword applies to the water phase.");

    }

private:
    // enable scaling of the input saturations (i.e., rescale the x-Axis)
    bool enableSatScaling_;

    // use three (instead of two) points to scale the saturations for the relative
    // permeabilities.
    //
    // This means that two piecewise linear functions are used for saturation scaling
    // instead of a single linear one
    bool enableThreePointKrSatScaling_;

    // enable the scaling of the capillary pressure and relative permeability outputs
    // (i.e., rescale the y-Axis)
    bool enablePcScaling_;
    bool enableLeverettScaling_;
    bool enableKrwScaling_;
    bool enableKrnScaling_;
};

} // namespace Ewoms

#endif
