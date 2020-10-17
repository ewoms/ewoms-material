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
 * \brief This is the unit test for the co2 brine PVT model
 *
 * This test requires the presence of opm-common.
 */
#include "config.h"

#if !HAVE_ECL_INPUT
#error "The test for the co2 brine PVT classes requires eclipse input support in opm-common"
#endif

//#include <ewoms/material/fluidsystems/blackoilpvt/co2gaspvt.hh>
//#include <ewoms/material/fluidsystems/blackoilpvt/brineco2pvt.hh>

#include <ewoms/material/fluidsystems/blackoilpvt/gaspvtmultiplexer.hh>
#include <ewoms/material/fluidsystems/blackoilpvt/oilpvtmultiplexer.hh>
#include <ewoms/material/fluidsystems/blackoilpvt/waterpvtmultiplexer.hh>

#include <ewoms/common/densead/evaluation.hh>
#include <ewoms/common/densead/math.hh>

#include <ewoms/eclio/parser/parser.hh>
#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>

#include <dune/common/parallel/mpihelper.hh>

namespace Ewoms {
namespace CO2DefaultTables {
#include <ewoms/material/components/co2tables.inc.cc>
}}

// values of strings based on the first SPE1 test case of opm-data.  note that in the
// real world it does not make much sense to specify a fluid phase using more than a
// single keyword, but for a unit test, this saves a lot of boiler-plate code.
static const char* deckString1 =
    "RUNSPEC\n"
    "\n"
    "DIMENS\n"
    "   10 10 3 /\n"
    "\n"
    "TABDIMS\n"
    " * 1 /\n"
    "\n"
    "OIL\n"
    "GAS\n"
    "CO2STOR\n"
    "\n"
    "DISGAS\n"
    "\n"
    "METRIC\n"
    "\n"
    "GRID\n"
    "\n"
    "DX\n"
    "   	300*1000 /\n"
    "DY\n"
    "	300*1000 /\n"
    "DZ\n"
    "	100*20 100*30 100*50 /\n"
    "\n"
    "TOPS\n"
    "	100*1234 /\n"
    "\n"
    "PORO\n"
    "  300*0.15 /\n"
    "PROPS\n"
    "\n";

template <class Evaluation, class BrinePvt, class Co2Pvt>
void ensurePvtApi(const BrinePvt& brinePvt, const Co2Pvt& co2Pvt)
{
    // we don't want to run this, we just want to make sure that it compiles
    while (0) {
        Evaluation temperature = 273.15 + 20.0;
        Evaluation pressure = 1e5;
        Evaluation Rs = 0.0;
        Evaluation Rv = 0.0;
        Evaluation So = 0.5;
        Evaluation maxSo = 1.0;
        Evaluation tmp;

        /////
        // brine PVT API
        /////
        tmp = brinePvt.viscosity(/*regionIdx=*/0,
                               temperature,
                               pressure,
                               Rs);
        tmp = brinePvt.inverseFormationVolumeFactor(/*regionIdx=*/0,
                                                  temperature,
                                                  pressure,
                                                  Rs);
        tmp = brinePvt.saturatedViscosity(/*regionIdx=*/0,
                                        temperature,
                                        pressure);
        tmp = brinePvt.saturatedInverseFormationVolumeFactor(/*regionIdx=*/0,
                                                           temperature,
                                                           pressure);
        tmp = brinePvt.saturationPressure(/*regionIdx=*/0,
                                        temperature,
                                        Rs);
        tmp = brinePvt.saturatedGasDissolutionFactor(/*regionIdx=*/0,
                                                   temperature,
                                                   pressure);
        tmp = brinePvt.saturatedGasDissolutionFactor(/*regionIdx=*/0,
                                                   temperature,
                                                   pressure,
                                                   So,
                                                   maxSo);

        /////
        // co2 PVT API
        /////
        tmp = co2Pvt.viscosity(/*regionIdx=*/0,
                               temperature,
                               pressure,
                               Rv);
        tmp = co2Pvt.inverseFormationVolumeFactor(/*regionIdx=*/0,
                                                  temperature,
                                                  pressure,
                                                  Rv);
        tmp = co2Pvt.saturatedViscosity(/*regionIdx=*/0,
                                        temperature,
                                        pressure);
        tmp = co2Pvt.saturatedInverseFormationVolumeFactor(/*regionIdx=*/0,
                                                           temperature,
                                                           pressure);
        tmp = co2Pvt.saturationPressure(/*regionIdx=*/0,
                                        temperature,
                                        Rv);
        tmp = co2Pvt.saturatedOilVaporizationFactor(/*regionIdx=*/0,
                                                    temperature,
                                                    pressure);
        tmp = co2Pvt.saturatedOilVaporizationFactor(/*regionIdx=*/0,
                                                    temperature,
                                                    pressure,
                                                    So,
                                                    maxSo);

        // prevent GCC from producing a "variable assigned but unused" warning
        tmp = 2.0*tmp;
    }
}

template <class Scalar>
inline void testAll()
{
    Ewoms::Parser parser;

    auto deck = parser.parseString(deckString1);
    Ewoms::EclipseState eclState(deck);
    Ewoms::Schedule schedule(deck, eclState);

    Ewoms::GasPvtMultiplexer<Scalar> co2Pvt;
    Ewoms::OilPvtMultiplexer<Scalar> brinePvt;

    co2Pvt.initFromEclState(eclState, schedule);
    brinePvt.initFromEclState(eclState, schedule);

    typedef Ewoms::DenseAd::Evaluation<Scalar, 1> FooEval;
    ensurePvtApi<Scalar>(brinePvt, co2Pvt);
    ensurePvtApi<FooEval>(brinePvt, co2Pvt);
}

int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);

    testAll<double>();
    testAll<float>();

    return 0;
}
