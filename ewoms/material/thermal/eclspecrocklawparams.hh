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
 * \copydoc Ewoms::EclSpecrockLawParams
 */
#ifndef EWOMS_ECL_SPECROCK_LAW_PARAMS_HH
#define EWOMS_ECL_SPECROCK_LAW_PARAMS_HH

#include <ewoms/material/common/ensurefinalized.hh>

#include <ewoms/common/tabulated1dfunction.hh>

namespace Ewoms {

/*!
 * \brief The default implementation of a parameter object for the
 *        ECL thermal law based on SPECROCK.
 */
template <class ScalarT>
class EclSpecrockLawParams : public EnsureFinalized
{
    typedef Tabulated1DFunction<ScalarT> InternalEnergyFunction;

public:
    typedef ScalarT Scalar;

    EclSpecrockLawParams(const EclSpecrockLawParams&) = default;

    EclSpecrockLawParams()
    { }

    /*!
     * \brief Specify the volumetric internal energy of rock via heat capacities.
     */
    template <class Container>
    void setHeatCapacities(const Container& temperature,
                           const Container& heatCapacity)
    {
        assert(temperature.size() == heatCapacity.size());

        // integrate the heat capacity to compute the internal energy
        Scalar curU = temperature[0]*heatCapacity[0];
        unsigned n = temperature.size();
        std::vector<Scalar> T(n);
        std::vector<Scalar> u(n);
        for (unsigned i = 0; i < temperature.size(); ++ i) {
            T[i] = temperature[i];
            u[i] = curU;

            if (i >= temperature.size() - 1)
                break;

            // integrate to the heat capacity from the current sampling point to the next
            // one. this leads to a quadratic polynomial.
            Scalar c_v0 = heatCapacity[i];
            Scalar c_v1 = heatCapacity[i + 1];
            Scalar T0 = temperature[i];
            Scalar T1 = temperature[i + 1];
            curU += 0.5*(c_v0 + c_v1)*(T1 - T0);
        }

        internalEnergyFunction_.setXYContainers(T, u);
    }

    /*!
     * \brief Return the function which maps temparature to the rock's volumetric
     *        internal energy
     *
     * Currently we assume this function to be piecewise linear. (Assuming piecewise
     * linear heat capacity, the real function is quadratic, but the difference should be
     * negligible.)
     */
    const InternalEnergyFunction& internalEnergyFunction() const
    { EnsureFinalized::check(); return internalEnergyFunction_; }

private:
    InternalEnergyFunction internalEnergyFunction_;
};

} // namespace Ewoms

#endif
