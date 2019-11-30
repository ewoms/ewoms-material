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
 * \copydoc Ewoms::EffToAbsLawParams
 */
#ifndef EWOMS_EFF_TO_ABS_LAW_PARAMS_HH
#define EWOMS_EFF_TO_ABS_LAW_PARAMS_HH

#include <cassert>

#include <ewoms/material/common/ensurefinalized.hh>

namespace Ewoms {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief A default implementation of the parameters for the adapter
 *        class to convert material laws from effective to absolute
 *        saturations.
 */
template <class EffLawParamsT, int numPhases>
class EffToAbsLawParams : public EffLawParamsT
{
    typedef EffLawParamsT EffLawParams;
    typedef typename EffLawParams::Traits::Scalar Scalar;

public:
    typedef typename EffLawParams::Traits Traits;

    EffToAbsLawParams()
        : EffLawParams()
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            residualSaturation_[phaseIdx] = 0.0;
    }

    /*!
     * \brief Calculate all dependent quantities once the independent
     *        quantities of the parameter object have been set.
     */
    void finalize()
    {
        sumResidualSaturations_ = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            sumResidualSaturations_ += residualSaturation_[phaseIdx];

        EffLawParams::finalize();
    }

    /*!
     * \brief Return the residual saturation of a phase.
     */
    Scalar residualSaturation(unsigned phaseIdx) const
    { EnsureFinalized::check(); return residualSaturation_[phaseIdx]; }

    /*!
     * \brief Return the sum of the residual saturations.
     */
    Scalar sumResidualSaturations() const
    { EnsureFinalized::check(); return sumResidualSaturations_; }

    /*!
     * \brief Set the residual saturation of a phase.
     */
    void setResidualSaturation(unsigned phaseIdx, Scalar value)
    { residualSaturation_[phaseIdx] = value; }

private:

    Scalar residualSaturation_[numPhases];
    Scalar sumResidualSaturations_;
};

} // namespace Ewoms

#endif
