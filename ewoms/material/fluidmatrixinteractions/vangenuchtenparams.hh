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
 * \copydoc Ewoms::VanGenuchtenParams
 */
#ifndef VAN_GENUCHTEN_PARAMS_HH
#define VAN_GENUCHTEN_PARAMS_HH

#include <ewoms/material/common/ensurefinalized.hh>

namespace Ewoms {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Specification of the material parameters for the van
 *        Genuchten constitutive relations.
 *
 * In this implementation setting either the \f$n\f$ or \f$m\f$ shape
 * parameter automatically calculates the other. I.e. they cannot be
 * set independently.
 */
template<class TraitsT>
class VanGenuchtenParams : public EnsureFinalized
{
    typedef typename TraitsT::Scalar Scalar;

public:
    using EnsureFinalized :: finalize;

    typedef TraitsT Traits;

    VanGenuchtenParams()
    {
    }

    VanGenuchtenParams(Scalar alphaParam, Scalar nParam)
    {
        setVgAlpha(alphaParam);
        setVgN(nParam);
        finalize();
    }

    /*!
     * \brief Return the \f$\alpha\f$ shape parameter of van Genuchten's
     *        curve.
     */
    Scalar vgAlpha() const
    { EnsureFinalized::check(); return vgAlpha_; }

    /*!
     * \brief Set the \f$\alpha\f$ shape parameter of van Genuchten's
     *        curve.
     */
    void setVgAlpha(Scalar v)
    { vgAlpha_ = v; }

    /*!
     * \brief Return the \f$m\f$ shape parameter of van Genuchten's
     *        curve.
     */
    Scalar vgM() const
    { EnsureFinalized::check(); return vgM_; }

    /*!
     * \brief Set the \f$m\f$ shape parameter of van Genuchten's
     *        curve.
     *
     * The \f$n\f$ shape parameter is set to \f$n = \frac{1}{1 - m}\f$
     */
    void setVgM(Scalar m)
    { vgM_ = m; vgN_ = 1/(1 - vgM_); }

    /*!
     * \brief Return the \f$n\f$ shape parameter of van Genuchten's
     *        curve.
     */
    Scalar vgN() const
    { EnsureFinalized::check(); return vgN_; }

    /*!
     * \brief Set the \f$n\f$ shape parameter of van Genuchten's
     *        curve.
     *
     * The \f$n\f$ shape parameter is set to \f$m = 1 - \frac{1}{n}\f$
     */
    void setVgN(Scalar n)
    { vgN_ = n; vgM_ = 1 - 1/vgN_; }

private:
    Scalar vgAlpha_;
    Scalar vgM_;
    Scalar vgN_;
};
} // namespace Ewoms

#endif
