# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:
cmake_minimum_required(VERSION 3.5)

# set up project and specify the minimum cmake version
project("ewoms-material" C CXX)

# find the build system (i.e., ewoms-common) and set cmake's module path
find_package(ewoms-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${ewoms-common_MODULE_PATH})

# include the eWoms cmake macros
include(EwomsMacros NO_POLICY_SCOPE)

# do most of the book-keeping required
ewoms_project()

# either std::optional or std::experimental::optional needs to be supported
find_package(StdOptional REQUIRED)

# we want all features detected by the build system to be enabled,
# thank you!
dune_enable_all_packages()

# recursively mark all header files beneath the "ewoms" directory for
# installation.
ewoms_recusive_export_all_headers("ewoms")

# this module is headers-only, i.e., it does not feature a library
# which must be linked a against, so we skip this here

# add unit tests
ewoms_add_test(test_blackoilfluidstate)
ewoms_add_test(test_eclblackoilfluidsystem CONDITION ewoms-eclio_FOUND)
ewoms_add_test(test_eclblackoilpvt CONDITION ewoms-eclio_FOUND)
ewoms_add_test(test_co2brinepvt CONDITION ewoms-eclio_FOUND)
ewoms_add_test(test_eclmateriallawmanager CONDITION ewoms-eclio_FOUND)
ewoms_add_test(test_fluidmatrixinteractions)
ewoms_add_test(test_pengrobinson)
ewoms_add_test(test_ncpflash)
ewoms_add_test(test_tabulation)
ewoms_add_test(test_components)
ewoms_add_test(test_fluidsystems)
ewoms_add_test(test_immiscibleflash)

# finalize the project, e.g. generate the config.h etc.
finalize_ewoms_project()
