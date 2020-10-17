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
 * \copydoc Ewoms::EclMultiplexerMaterialParams
 */
#ifndef EWOMS_ECL_MULTIPLEXER_MATERIAL_PARAMS_HH
#define EWOMS_ECL_MULTIPLEXER_MATERIAL_PARAMS_HH

#include "eclstone1material.hh"
#include "eclstone2material.hh"
#include "ecldefaultmaterial.hh"
#include "ecltwophasematerial.hh"

#include <type_traits>
#include <cassert>
#include <memory>

#include <ewoms/material/common/ensurefinalized.hh>

namespace Ewoms {

enum class EclMultiplexerApproach {
    EclDefaultApproach,
    EclStone1Approach,
    EclStone2Approach,
    EclTwoPhaseApproach,
    EclOnePhaseApproach
};

/*!
 * \brief Multiplexer implementation for the parameters required by the
 *        multiplexed three-phase material law.
 *
 * Essentially, this class just stores parameter object for the "nested" material law and
 * provides some methods to convert to it.
 */
template<class Traits, class GasOilMaterialLawT, class OilWaterMaterialLawT>
class EclMultiplexerMaterialParams : public Traits, public EnsureFinalized
{

    typedef Ewoms::EclStone1Material<Traits, GasOilMaterialLawT, OilWaterMaterialLawT> Stone1Material;
    typedef Ewoms::EclStone2Material<Traits, GasOilMaterialLawT, OilWaterMaterialLawT> Stone2Material;
    typedef Ewoms::EclDefaultMaterial<Traits, GasOilMaterialLawT, OilWaterMaterialLawT> DefaultMaterial;
    typedef Ewoms::EclTwoPhaseMaterial<Traits, GasOilMaterialLawT, OilWaterMaterialLawT> TwoPhaseMaterial;

    typedef typename Stone1Material::Params Stone1Params;
    typedef typename Stone2Material::Params Stone2Params;
    typedef typename DefaultMaterial::Params DefaultParams;
    typedef typename TwoPhaseMaterial::Params TwoPhaseParams;

    template <class ParamT>
    struct Deleter
    {
        inline void operator () ( void* ptr )
        {
            delete static_cast< ParamT* > (ptr);
        }
    };

    typedef std::shared_ptr< void > ParamPointerType;

public:
    using EnsureFinalized :: finalize;

    /*!
     * \brief The multiplexer constructor.
     */
    EclMultiplexerMaterialParams() : realParams_()
    {
    }

    EclMultiplexerMaterialParams(const EclMultiplexerMaterialParams& other)
        : realParams_()
    {
        setApproach( other.approach() );
    }

    EclMultiplexerMaterialParams& operator= ( const EclMultiplexerMaterialParams& other )
    {
        realParams_.reset();
        setApproach( other.approach() );
        return *this;
    }

    void setApproach(EclMultiplexerApproach newApproach)
    {
        assert(realParams_ == 0);
        approach_ = newApproach;

        switch (approach()) {
        case EclMultiplexerApproach::EclStone1Approach:
            realParams_ = ParamPointerType(new Stone1Params, Deleter< Stone1Params > () );
            break;

        case EclMultiplexerApproach::EclStone2Approach:
            realParams_ = ParamPointerType(new Stone2Params, Deleter< Stone2Params > () );
            break;

        case EclMultiplexerApproach::EclDefaultApproach:
            realParams_ = ParamPointerType(new DefaultParams, Deleter< DefaultParams > () );
            break;

        case EclMultiplexerApproach::EclTwoPhaseApproach:
            realParams_ = ParamPointerType(new TwoPhaseParams, Deleter< TwoPhaseParams > () );
            break;

        case EclMultiplexerApproach::EclOnePhaseApproach:
            // Do nothing, no parameters.
            break;
        }
    }

    EclMultiplexerApproach approach() const
    { return approach_; }

    // get the parameter object for the Stone1 case
    template <EclMultiplexerApproach approachV>
    typename std::enable_if<approachV == EclMultiplexerApproach::EclStone1Approach, Stone1Params>::type&
    getRealParams()
    {
        assert(approach() == approachV);
        return this->template castTo<Stone1Params>();
    }

    template <EclMultiplexerApproach approachV>
    typename std::enable_if<approachV == EclMultiplexerApproach::EclStone1Approach, const Stone1Params>::type&
    getRealParams() const
    {
        assert(approach() == approachV);
        return this->template castTo<Stone1Params>();
    }

    // get the parameter object for the Stone2 case
    template <EclMultiplexerApproach approachV>
    typename std::enable_if<approachV == EclMultiplexerApproach::EclStone2Approach, Stone2Params>::type&
    getRealParams()
    {
        assert(approach() == approachV);
        return this->template castTo<Stone2Params>();
    }

    template <EclMultiplexerApproach approachV>
    typename std::enable_if<approachV == EclMultiplexerApproach::EclStone2Approach, const Stone2Params>::type&
    getRealParams() const
    {
        assert(approach() == approachV);
        return this->template castTo<Stone2Params>();
    }

    // get the parameter object for the default case
    template <EclMultiplexerApproach approachV>
    typename std::enable_if<approachV == EclMultiplexerApproach::EclDefaultApproach, DefaultParams>::type&
    getRealParams()
    {
        assert(approach() == approachV);
        return this->template castTo<DefaultParams>();
    }

    template <EclMultiplexerApproach approachV>
    typename std::enable_if<approachV == EclMultiplexerApproach::EclDefaultApproach, const DefaultParams>::type&
    getRealParams() const
    {
        assert(approach() == approachV);
        return this->template castTo<DefaultParams>();
    }

    // get the parameter object for the twophase case
    template <EclMultiplexerApproach approachV>
    typename std::enable_if<approachV == EclMultiplexerApproach::EclTwoPhaseApproach, TwoPhaseParams>::type&
    getRealParams()
    {
        assert(approach() == approachV);
        return this->template castTo<TwoPhaseParams>();
    }

    template <EclMultiplexerApproach approachV>
    typename std::enable_if<approachV == EclMultiplexerApproach::EclTwoPhaseApproach, const TwoPhaseParams>::type&
    getRealParams() const
    {
        assert(approach() == approachV);
        return this->template castTo<TwoPhaseParams>();
    }

private:
    template <class ParamT>
    ParamT& castTo()
    {
        return *(static_cast<ParamT *> (realParams_.operator->()));
    }

    template <class ParamT>
    const ParamT& castTo() const
    {
        return *(static_cast<const ParamT *> (realParams_.operator->()));
    }

    EclMultiplexerApproach approach_;
    ParamPointerType realParams_;
};
} // namespace Ewoms

#endif
