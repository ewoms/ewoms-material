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
 * \brief This test makes sure that mandated API is adhered to by all component classes
 */
#include "config.h"

#include <ewoms/common/densead/evaluation.hh>
#include <ewoms/common/densead/math.hh>

#include "checkcomponent.hh"

// include all components shipped with ewoms-material
#include <ewoms/material/components/unit.hh>
#include <ewoms/material/components/nullcomponent.hh>
#include <ewoms/material/components/component.hh>
#include <ewoms/material/components/dnapl.hh>
#include <ewoms/material/components/simpleh2o.hh>
#include <ewoms/material/components/lnapl.hh>
#include <ewoms/material/components/iapws/region2.hh>
#include <ewoms/material/components/iapws/region1.hh>
#include <ewoms/material/components/iapws/common.hh>
#include <ewoms/material/components/iapws/region4.hh>
#include <ewoms/material/components/h2o.hh>
#include <ewoms/material/components/co2.hh>
#include <ewoms/material/components/mesitylene.hh>
#include <ewoms/material/components/tabulatedcomponent.hh>
#include <ewoms/material/components/brine.hh>
#include <ewoms/material/components/n2.hh>
#include <ewoms/material/components/xylene.hh>
#include <ewoms/material/components/air.hh>
#include <ewoms/material/components/simpleco2.hh>

#include <ewoms/common/uniformtabulated2dfunction.hh>

namespace Ewoms {
namespace ComponentsTest {
#include <ewoms/material/components/co2tables.inc>
}}

#include <dune/common/parallel/mpihelper.hh>

template <class Scalar, class Evaluation>
void testAllComponents()
{
    typedef Ewoms::H2O<Scalar> H2O;

    checkComponent<Ewoms::Air<Scalar>, Evaluation>();
    checkComponent<Ewoms::Brine<Scalar, H2O>, Evaluation>();
    checkComponent<Ewoms::CO2<Scalar, Ewoms::ComponentsTest::CO2Tables>, Evaluation>();
    checkComponent<Ewoms::DNAPL<Scalar>, Evaluation>();
    checkComponent<Ewoms::H2O<Scalar>, Evaluation>();
    checkComponent<Ewoms::LNAPL<Scalar>, Evaluation>();
    checkComponent<Ewoms::Mesitylene<Scalar>, Evaluation>();
    checkComponent<Ewoms::N2<Scalar>, Evaluation>();
    checkComponent<Ewoms::NullComponent<Scalar>, Evaluation>();
    checkComponent<Ewoms::SimpleCO2<Scalar>, Evaluation>();
    checkComponent<Ewoms::SimpleH2O<Scalar>, Evaluation>();
    checkComponent<Ewoms::TabulatedComponent<Scalar, H2O>, Evaluation>();
    checkComponent<Ewoms::Unit<Scalar>, Evaluation>();
    checkComponent<Ewoms::Xylene<Scalar>, Evaluation>();
}

template <class Scalar>
inline void testAll()
{
    typedef Ewoms::DenseAd::Evaluation<Scalar, 3> Evaluation;

    // ensure that all components are API-compliant
    testAllComponents<Scalar, Scalar>();
    testAllComponents<Scalar, Evaluation>();
}

int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);

    testAll<double>();
    testAll<float>();

    return 0;
}
