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
 *
 * \brief Modules for the ModularFluidState which represent viscosity.
 */
#ifndef EWOMS_FLUID_STATE_VISCOSITY_MODULES_HH
#define EWOMS_FLUID_STATE_VISCOSITY_MODULES_HH

#include <ewoms/common/mathtoolbox.hh>
#include <ewoms/common/valgrind.hh>

#include <algorithm>

namespace Ewoms {
/*!
 * \brief Module for the modular fluid state which stores the
 *       viscosities explicitly.
 */
template <class Scalar,
          unsigned numPhases,
          class Implementation>
class FluidStateExplicitViscosityModule
{
public:
    FluidStateExplicitViscosityModule()
    { Valgrind::SetUndefined(viscosity_); }

    /*!
     * \brief The viscosity of a fluid phase [-]
     */
    const Scalar& viscosity(unsigned phaseIdx) const
    { return viscosity_[phaseIdx]; }

    /*!
     * \brief Set the dynamic viscosity of a phase [Pa s]
     */
    void setViscosity(unsigned phaseIdx, Scalar value)
    { viscosity_[phaseIdx] = value; }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            viscosity_[phaseIdx] = Ewoms::decay<Scalar>(fs.viscosity(phaseIdx));
        }
    }

    /*!
     * \brief Make sure that all attributes are defined.
     *
     * This method does not do anything if the program is not run
     * under valgrind. If it is, then valgrind will print an error
     * message if some attributes of the object have not been properly
     * defined.
     */
    void checkDefined() const
    {
        Valgrind::CheckDefined(viscosity_);
    }

protected:
    Scalar viscosity_[numPhases];
};

/*!
 * \brief Module for the modular fluid state which does not  the
 *        viscosities but throws std::logic_error instead.
 */
template <class Scalar,
          unsigned numPhases,
          class Implementation>
class FluidStateNullViscosityModule
{
public:
    FluidStateNullViscosityModule()
    { }

    /*!
     * \brief The viscosity of a fluid phase [-]
     */
    const Scalar& viscosity(unsigned /* phaseIdx */) const
    { throw std::logic_error("Viscosity is not provided by this fluid state"); }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& /* fs */)
    { }

    /*!
     * \brief Make sure that all attributes are defined.
     *
     * This method does not do anything if the program is not run
     * under valgrind. If it is, then valgrind will print an error
     * message if some attributes of the object have not been properly
     * defined.
     */
    void checkDefined() const
    { }
};

} // namespace Ewoms

#endif
