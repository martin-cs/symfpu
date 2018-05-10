/*
** Copyright (C) 2018 Martin Brain
**
** See the file LICENSE for licensing information.
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

