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
 * \copydoc Ewoms::LinearMaterialParams
 */
#ifndef EWOMS_LINEAR_MATERIAL_PARAMS_HH
#define EWOMS_LINEAR_MATERIAL_PARAMS_HH

#include <cassert>

#include <ewoms/material/common/ensurefinalized.hh>

namespace Ewoms {

/*!
 * \brief Reference implementation of params for the linear M-phase
 *        material material.
 */
template<class TraitsT>
class LinearMaterialParams : public EnsureFinalized
{
    enum { numPhases = TraitsT::numPhases };

    typedef typename TraitsT::Scalar Scalar;

public:
    using EnsureFinalized :: finalize;

    typedef TraitsT Traits;

    /*!
     * \brief The default constructor.
     *
     * We set the capillary pressure to zero, if not specified otherwise.
     */
    LinearMaterialParams()
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            setPcMinSat(phaseIdx, 0.0);
            setPcMaxSat(phaseIdx, 0.0);
        }
    }

    /*!
     * \brief Return the relative phase pressure at the minimum saturation of a phase.
     *
     * This means \f$p_{c\alpha}\f$ at \f$S_\alpha=0\f$.
     */
    Scalar pcMinSat(unsigned phaseIdx) const
    { EnsureFinalized::check();return pcMinSat_[phaseIdx]; }

    /*!
     * \brief Set the relative phase pressure at the minimum saturation of a phase.
     *
     * This means \f$p_{c\alpha}\f$ at \f$S_\alpha=0\f$.
     */
    void setPcMinSat(unsigned phaseIdx, Scalar val)
    { pcMinSat_[phaseIdx] = val; }

    /*!
     * \brief Return the relative phase pressure at the maximum saturation of a phase.
     *
     * This means \f$p_{c\alpha}\f$ at \f$S_\alpha=1\f$.
     */
    Scalar pcMaxSat(unsigned phaseIdx) const
    { EnsureFinalized::check(); return pcMaxSat_[phaseIdx]; }

    /*!
     * \brief Set the relative phase pressure at the maximum saturation of a phase.
     *
     * This means \f$p_{c\alpha}\f$ at \f$S_\alpha=1\f$.
     */
    void setPcMaxSat(unsigned phaseIdx, Scalar val)
    { pcMaxSat_[phaseIdx] = val; }

private:
    Scalar pcMaxSat_[numPhases];
    Scalar pcMinSat_[numPhases];
};
} // namespace Ewoms

#endif
