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
 * \brief Modules for the ModularFluidState which represent density.
 */
#ifndef EWOMS_FLUID_STATE_DENSITY_MODULES_HH
#define EWOMS_FLUID_STATE_DENSITY_MODULES_HH

#include <ewoms/common/mathtoolbox.hh>
#include <ewoms/common/valgrind.hh>

#include <algorithm>

namespace Ewoms {

/*!
 * \brief Module for the modular fluid state which stores the
 *       densities explicitly.
 */
template <class Scalar,
          unsigned numPhases,
          class Implementation>
class FluidStateExplicitDensityModule
{
public:
    FluidStateExplicitDensityModule()
    { Valgrind::SetUndefined(density_); }

    /*!
     * \brief The density of a fluid phase [kg/m^3]
     */
    const Scalar& density(unsigned phaseIdx) const
    { return density_[phaseIdx]; }

    /*!
     * \brief The molar density of a fluid phase [mol/m^3]
     */
    Scalar molarDensity(unsigned phaseIdx) const
    { return density_[phaseIdx]/asImp_().averageMolarMass(phaseIdx); }

    /*!
     * \brief The molar volume of a fluid phase [m^3/mol]
     */
    Scalar molarVolume(unsigned phaseIdx) const
    { return 1/molarDensity(phaseIdx); }

    /*!
     * \brief Set the density of a phase [kg/m^3]
     */
    void setDensity(unsigned phaseIdx, const Scalar&  value)
    { density_[phaseIdx] = value; }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs)
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            density_[phaseIdx] = Ewoms::decay<Scalar>(fs.density(phaseIdx));
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
        Valgrind::CheckDefined(density_);
    }

protected:
    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    Scalar density_[numPhases];
};

/*!
 * \brief Module for the modular fluid state which does not  the
 *        densities but throws std::logic_error instead.
 */
template <class Scalar,
          unsigned numPhases,
          class Implementation>
class FluidStateNullDensityModule
{
public:
    FluidStateNullDensityModule()
    { }

    /*!
     * \brief The density of a fluid phase [kg/m^3]
     */
    const Scalar& density(unsigned /* phaseIdx */) const
    { throw std::logic_error("Density is not provided by this fluid state"); }

    /*!
     * \brief The molar density of a fluid phase [mol/m^3]
     */
    const Scalar& molarDensity(unsigned /* phaseIdx */) const
    { throw std::logic_error("Molar density is not provided by this fluid state"); }

    /*!
     * \brief The molar volume of a fluid phase [m^3/mol]
     */
    const Scalar& molarVolume(unsigned /* phaseIdx */) const
    { throw std::logic_error("Molar volume is not provided by this fluid state"); }

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
