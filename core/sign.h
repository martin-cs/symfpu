/*
** Copyright (C) 2018 Martin Brain
**
** See the file LICENSE for licensing information.
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
