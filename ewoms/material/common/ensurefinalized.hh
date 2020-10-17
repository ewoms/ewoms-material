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
 * \copydoc Ewoms::EnsureFinalized
 */
#ifndef EWOMS_MATERIAL_ENSURE_FINALIZED_HH
#define EWOMS_MATERIAL_ENSURE_FINALIZED_HH

#include <cassert>
#include <stdexcept>

// TODO: move this variable to config.h
#define EWOMS_CHECK_PARAM_FINALIZED 1

#if ! defined(NDEBUG) && EWOMS_CHECK_PARAM_FINALIZED
#define USE_EWOMS_CHECK_PARAM_FINALIZED 1
#endif

namespace Ewoms {

/*!
 * \brief Default implementation for asserting finalization of parameter objects.
 *
 */
class EnsureFinalized
{
#if USE_EWOMS_CHECK_PARAM_FINALIZED
    bool finalized_;
#endif

protected:
    /*!
     * \brief The default constructor.
     */
    EnsureFinalized()
#if USE_EWOMS_CHECK_PARAM_FINALIZED
        : finalized_( false )
#endif
    {
    }

    void check() const
    {
#if USE_EWOMS_CHECK_PARAM_FINALIZED
        if (!finalized_)
            throw std::runtime_error("Parameter class has not been finalized before usage!");
#endif
    }

public:
    /*!
     * \brief Mark the object as finalized.
     */
    void finalize()
    {
#if USE_EWOMS_CHECK_PARAM_FINALIZED
        finalized_ = true;
#endif
    }
};

#undef USE_EWOMS_CHECK_PARAM_FINALIZED

} // namespace Ewoms
#endif
