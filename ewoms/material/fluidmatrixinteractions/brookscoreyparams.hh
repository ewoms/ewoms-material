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
 * \copydoc Ewoms::BrooksCoreyParams
 */
#ifndef EWOMS_BROOKS_COREY_PARAMS_HH
#define EWOMS_BROOKS_COREY_PARAMS_HH

#include <ewoms/common/valgrind.hh>
#include <ewoms/material/common/ensurefinalized.hh>

#include <cassert>

namespace Ewoms {

/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Specification of the material parameters for the
 *        Brooks-Corey constitutive relations.
 *
 *\see BrooksCorey
 */
template <class TraitsT>
class BrooksCoreyParams : public EnsureFinalized
{
    typedef typename TraitsT::Scalar Scalar;
public:
    using EnsureFinalized :: finalize;

    typedef TraitsT Traits;

    BrooksCoreyParams()
    {
        Valgrind::SetUndefined(*this);
    }

    BrooksCoreyParams(Scalar ePressure, Scalar shapeParam)
        : entryPressure_(ePressure), lambda_(shapeParam)
    {
        finalize();
    }

    /*!
     * \brief Returns the entry pressure [Pa]
     */
    Scalar entryPressure() const
    { EnsureFinalized::check(); return entryPressure_; }

    /*!
     * \brief Set the entry pressure [Pa]
     */
    void setEntryPressure(Scalar v)
    { entryPressure_ = v; }

    /*!
     * \brief Returns the lambda shape parameter
     */
    Scalar lambda() const
    { EnsureFinalized::check(); return lambda_; }

    /*!
     * \brief Set the lambda shape parameter
     */
    void setLambda(Scalar v)
    { lambda_ = v; }

private:
    Scalar entryPressure_;
    Scalar lambda_;
};
} // namespace Ewoms

#endif
