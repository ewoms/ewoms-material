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
 * \copydoc Ewoms::EclHysteresisConfig
 */
#ifndef EWOMS_ECL_HYSTERESIS_CONFIG_HH
#define EWOMS_ECL_HYSTERESIS_CONFIG_HH

#if HAVE_ECL_INPUT
#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/deck/deckkeyword.hh>
#include <ewoms/eclio/parser/deck/deckrecord.hh>
#include <ewoms/eclio/parser/deck/deckitem.hh>
#endif

#include <ewoms/common/exceptions.hh>

#include <string>
#include <cassert>
#include <algorithm>

namespace Ewoms {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Specifies the configuration used by the ECL kr/pC hysteresis code
 */
class EclHysteresisConfig
{
public:
    EclHysteresisConfig()
    {
        enableHysteresis_ = false;
        pcHysteresisModel_ = 0;
        krHysteresisModel_ = 0;
    }

    /*!
     * \brief Specify whether hysteresis is enabled or not.
     */
    void setEnableHysteresis(bool yesno)
    { enableHysteresis_ = yesno; }

    /*!
     * \brief Returns whether hysteresis is enabled.
     */
    bool enableHysteresis() const
    { return enableHysteresis_; }

    /*!
     * \brief Set the type of the hysteresis model which is used for capillary pressure.
     *
     * -1: capillary pressure hysteresis is disabled
     * 0: use the Killough model for capillary pressure hysteresis
     */
    void setPcHysteresisModel(int value)
    { pcHysteresisModel_ = value; }

    /*!
     * \brief Return the type of the hysteresis model which is used for capillary pressure.
     *
     * -1: capillary pressure hysteresis is disabled
     * 0: use the Killough model for capillary pressure hysteresis
     */
    int pcHysteresisModel() const
    { return pcHysteresisModel_; }

    /*!
     * \brief Set the type of the hysteresis model which is used for relative permeability.
     *
     * -1: relperm hysteresis is disabled
     * 0: use the Carlson model for relative permeability hysteresis of the non-wetting
     *    phase and the drainage curve for the relperm of the wetting phase
     * 1: use the Carlson model for relative permeability hysteresis of the non-wetting
     *    phase and the imbibition curve for the relperm of the wetting phase
     */
    void setKrHysteresisModel(int value)
    { krHysteresisModel_ = value; }

    /*!
     * \brief Return the type of the hysteresis model which is used for relative permeability.
     *
     * -1: relperm hysteresis is disabled
     * 0: use the Carlson model for relative permeability hysteresis
     */
    int krHysteresisModel() const
    { return krHysteresisModel_; }

#if HAVE_ECL_INPUT
    /*!
     * \brief Reads all relevant material parameters form a cell of a parsed ECL deck.
     *
     * This requires that the ewoms-eclio module is available.
     */
    void initFromDeck(const Ewoms::Deck& deck)
    {
        enableHysteresis_ = false;

        if (!deck.hasKeyword("SATOPTS"))
            return;

        const auto& satoptsItem = deck.getKeyword("SATOPTS").getRecord(0).getItem(0);
        for (unsigned i = 0; i < satoptsItem.data_size(); ++i) {
            std::string satoptsValue = satoptsItem.get< std::string >(0);
            std::transform(satoptsValue.begin(),
                           satoptsValue.end(),
                           satoptsValue.begin(),
                           ::toupper);

            if (satoptsValue == "HYSTER")
                enableHysteresis_ = true;
        }

        // check for the (deprecated) HYST keyword
        if (deck.hasKeyword("HYST"))
            enableHysteresis_ = true;

        if (!enableHysteresis_)
            return;

        if (!deck.hasKeyword("EHYSTR"))
            throw std::runtime_error("Enabling hysteresis via the HYST parameter for SATOPTS requires the "
                                     "presence of the EHYSTR keyword");

        const auto& ehystrKeyword = deck.getKeyword("EHYSTR");
        if (deck.hasKeyword("NOHYKR"))
            krHysteresisModel_ = -1;
        else {
            krHysteresisModel_ = ehystrKeyword.getRecord(0).getItem("relative_perm_hyst").get<int>(0);

            if (krHysteresisModel_ != 0 && krHysteresisModel_ != 1)
                throw std::runtime_error(
                    "Only the Carlson relative permeability hystersis models (indicated by '0' or "
                    "'1' for the second item of the 'EHYSTR' keyword) are supported");
        }

        // this is slightly screwed: it is possible to specify contradicting hysteresis
        // models with HYPC/NOHYPC and the fifth item of EHYSTR. Let's ignore that for
        // now.
        std::string whereFlag =
            ehystrKeyword.getRecord(0).getItem("limiting_hyst_flag").getTrimmedString(0);
        if (deck.hasKeyword("NOHYPC") || whereFlag == "KR")
            pcHysteresisModel_ = -1;
        else {
            // if capillary pressure hysteresis is enabled, Eclipse always uses the
            // Killough model
            pcHysteresisModel_ = 0;

            throw std::runtime_error("Capillary pressure hysteresis is not supported yet");
        }
    }
#endif

    template<class Serializer>
    std::size_t packSize(Serializer& serializer) const
    {
        return serializer.packSize(enableHysteresis_) +
               serializer.packSize(pcHysteresisModel_) +
               serializer.packSize(krHysteresisModel_);
    }

    template<class Serializer>
    void pack(std::vector<char>& buffer, int& position,
              Serializer& serializer) const
    {
        serializer.pack(enableHysteresis_, buffer, position);
        serializer.pack(pcHysteresisModel_, buffer, position);
        serializer.pack(krHysteresisModel_, buffer, position);
    }

    template<class Serializer>
    void unpack(std::vector<char>& buffer, int& position,
                Serializer& serializer)
    {
        serializer.unpack(enableHysteresis_, buffer, position);
        serializer.unpack(pcHysteresisModel_, buffer, position);
        serializer.unpack(krHysteresisModel_, buffer, position);
    }

private:
    // enable hysteresis at all
    bool enableHysteresis_;

    // the capillary pressure and the relperm hysteresis models to be used
    int pcHysteresisModel_;
    int krHysteresisModel_;
};

} // namespace Ewoms

#endif
