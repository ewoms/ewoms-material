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
 * \copydoc Ewoms::SomertonThermalConductionLawParams
 */
#ifndef EWOMS_SOMERTON_THERMAL_CONDUCTION_LAW_PARAMS_HH
#define EWOMS_SOMERTON_THERMAL_CONDUCTION_LAW_PARAMS_HH

#include <cassert>

namespace Ewoms {

/*!
 * \brief The default implementation of a parameter object for the
 *        Somerton thermal conduction law.
 */
template <unsigned numPhases, class ScalarT>
class SomertonThermalConductionLawParams
{
    // do not copy!
    SomertonThermalConductionLawParams(const SomertonThermalConductionLawParams&)
    {}

public:
    typedef ScalarT Scalar;

    SomertonThermalConductionLawParams()
    { }

    /*!
     * \brief Return the "fully saturated" thermal conductivity of the
     *        porous medium [W/m^2 / (K/m)].
     *
     * In this context "fully saturated" means that the whole pore
     * space of the porous medium is filled by a given fluid phase.
     */
    Scalar fullySaturatedLambda(unsigned phaseIdx) const
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return fullySaturatedLambda_[phaseIdx];
    }

    /*!
     * \brief Set the "fully saturated" thermal conductivity of the
     *        porous medium [W/m^2 / (K/m)].
     *
     * In this context "fully saturated" means that the whole pore
     * space of the porous medium is filled by a given fluid phase.
     */
    void setFullySaturatedLambda(unsigned phaseIdx, Scalar value)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(value > 0);

        fullySaturatedLambda_[phaseIdx] = value;
    }

    /*!
     * \brief Return the thermal conductivity of the porous medium at
     *        vacuum [W/m^2 / (K/m)].
     */
    Scalar vacuumLambda() const
    {
        return vacuumLambda_;
    }

    /*!
     * \brief Set the "fully saturated" thermal conductivity of the
     *        porous medium [W/m^2 / (K/m)].
     *
     * In this context "fully saturated" means that the whole pore
     * space of the porous medium is filled by a given fluid phase.
     */
    void setVacuumLambda(Scalar value)
    {
        assert(value > 0);

        vacuumLambda_ = value;
    }

private:
    Scalar fullySaturatedLambda_[numPhases];
    Scalar vacuumLambda_;
};

} // namespace Ewoms

#endif
