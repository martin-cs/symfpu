/*
** Copyright (C) 2018 Martin Brain
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
** 
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
** sign.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 21/08/14
**
** The sign manipulating operations.
**
*/

#include "symfpu/core/unpackedFloat.h"

#ifndef SYMFPU_SIGN
#define SYMFPU_SIGN

namespace symfpu {

template <class t>
  unpackedFloat<t> negate (const typename t::fpt &format, const unpackedFloat<t> &uf) {

  PRECONDITION(uf.valid(format));

  unpackedFloat<t> result(uf, !uf.getSign());

  POSTCONDITION(result.valid(format));

  return result;
 }

template <class t>
  unpackedFloat<t> absolute (const typename t::fpt &format, const unpackedFloat<t> &uf) {

  PRECONDITION(uf.valid(format));

  unpackedFloat<t> result(uf, typename t::prop(false));

  POSTCONDITION(result.valid(format));

  return result;
 }


}

#endif
