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
 * \copydoc Ewoms::ConstantSolidHeatCapLawParams
 */
#ifndef EWOMS_CONSTANT_SOLID_HEAT_CAP_LAW_PARAMS_HH
#define EWOMS_CONSTANT_SOLID_HEAT_CAP_LAW_PARAMS_HH

#include <ewoms/material/common/ensurefinalized.hh>

namespace Ewoms {

/*!
 * \brief The default implementation of a parameter object for the
 *        solid energy storage law which assumes constant heat capacity.
 */
template <class ScalarT>
class ConstantSolidHeatCapLawParams : public EnsureFinalized
{
public:
    typedef ScalarT Scalar;

    ConstantSolidHeatCapLawParams(const ConstantSolidHeatCapLawParams&) = default;

    ConstantSolidHeatCapLawParams()
    { }

    /*!
     * \brief Set the specific heat capacity of the solid matrix [J/(m^3 K)].
     */
    void setSolidHeatCapacity(Scalar value)
    { solidHeatCapacity_ = value; }

    /*!
     * \brief Return the specific heat capacity of the solid matrix  [J/(m^3 K)].
     */
    Scalar solidHeatCapacity() const
    { EnsureFinalized::check(); return solidHeatCapacity_; }

private:
    Scalar solidHeatCapacity_;
};

} // namespace Ewoms

#endif
