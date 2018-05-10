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
** classify.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 21/08/14
**
** The classification functions for different classes of float.
**
*/

#include "symfpu/core/unpackedFloat.h"

#ifndef SYMFPU_CLASSIFY
#define SYMFPU_CLASSIFY

namespace symfpu {

template <class t>
  typename t::prop isNormal (const typename t::fpt &format, const unpackedFloat<t> &uf) {
  PRECONDITION(uf.valid(format));

  return !uf.getNaN() && !uf.getInf() && !uf.getZero() && uf.inNormalRange(format, typename t::prop(true));
 }


template <class t>
  typename t::prop isSubnormal (const typename t::fpt &format, const unpackedFloat<t> &uf) {
  PRECONDITION(uf.valid(format));

  return !uf.getNaN() && !uf.getInf() && !uf.getZero() && uf.inSubnormalRange(format, typename t::prop(true));
 }


template <class t>
  typename t::prop isZero (const typename t::fpt &format, const unpackedFloat<t> &uf) {
  PRECONDITION(uf.valid(format));

  return uf.getZero();
 }


template <class t>
  typename t::prop isInfinite (const typename t::fpt &format, const unpackedFloat<t> &uf) {
  PRECONDITION(uf.valid(format));

  return uf.getInf();
 }


template <class t>
  typename t::prop isNaN (const typename t::fpt &format, const unpackedFloat<t> &uf) {
  PRECONDITION(uf.valid(format));

  return uf.getNaN();
 }


// Note these are the SMT-LIB semantics, NaN is neither positive or negative

template <class t>
  typename t::prop isPositive (const typename t::fpt &format, const unpackedFloat<t> &uf) {
  PRECONDITION(uf.valid(format));

  return !uf.getNaN() && !uf.getSign();
 }

template <class t>
  typename t::prop isNegative (const typename t::fpt &format, const unpackedFloat<t> &uf) {
  PRECONDITION(uf.valid(format));

  return !uf.getNaN() && uf.getSign();
 }


// C semantics
 
 template <class t>
  typename t::prop isFinite (const typename t::fpt &format, const unpackedFloat<t> &uf) {
  PRECONDITION(uf.valid(format));

  return !uf.getNaN() && !uf.getInf();
 }


}

#endif

