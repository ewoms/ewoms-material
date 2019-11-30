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
 *
 * \brief This test ensures that the API of the black-oil fluid state conforms to the
 *        fluid state specification
 */
#include "config.h"

#include <ewoms/common/densead/evaluation.hh>
#include <ewoms/common/densead/math.hh>
#include <ewoms/material/fluidstates/blackoilfluidstate.hh>
#include <ewoms/material/fluidsystems/blackoilfluidsystem.hh>
#include <ewoms/material/checkfluidsystem.hh>

#include <dune/common/parallel/mpihelper.hh>

int main()
{
    {
        typedef double Scalar;
        typedef double Evaluation;
        typedef typename Ewoms::BlackOilFluidSystem<Scalar> FluidSystem;
        typedef Ewoms::BlackOilFluidState<Scalar, FluidSystem> FluidState;

        FluidState fs;
        checkFluidState<Evaluation>(fs);
    }

    {
        typedef float Scalar;
        typedef float Evaluation;
        typedef typename Ewoms::BlackOilFluidSystem<Scalar> FluidSystem;
        typedef Ewoms::BlackOilFluidState<Scalar, FluidSystem> FluidState;

        FluidState fs;
        checkFluidState<Evaluation>(fs);
    }

    {
        typedef float Scalar;
        typedef Ewoms::DenseAd::Evaluation<Scalar, 2> Evaluation;
        typedef typename Ewoms::BlackOilFluidSystem<Scalar> FluidSystem;
        typedef Ewoms::BlackOilFluidState<Evaluation, FluidSystem> FluidState;

        FluidState fs;
        checkFluidState<Evaluation>(fs);
    }

    return 0;
}
