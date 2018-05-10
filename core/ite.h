/*
** Copyright (C) 2018 Martin Brain
**
** See the file LICENSE for licensing information.
*/

/*
** ite.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 03/12/14
**
** The handling of if-then-else depends on the type of both the condition
** (i.e. is it concrete or is it symbolic) and the type of the data
** (i.e. is it a symbolic class that wraps a tree node, is it a structure 
** containing them, etc.).  Some of the handling will work for any condition,
** any of the data, etc.  So far not difficult in a language with multiple
** dispatch including type variables.
**
** However, we are using C++, so there are a few extra hoops to jump
** through.  We declare ITEs as a struct containing a static function
** as then we can partially specialise templates of them.
** Specialisations of these are then given when appropriate types are
** introduced.  Care must be taken with these to ensure that we don't
** get ambigious template instantiations.
**
*/

#ifndef SYMFPU_ITE
#define SYMFPU_ITE

namespace symfpu {

  template <class prop, class data>
    struct ite;

  // To avoid the need for putting the types *everywhere* we use a
  // helper function as C++ can perform type inference for functions
  // but not classes.
  template <class prop, class data>
    const data ITE (const prop &c, const data &l, const data &r) {
    return ite<prop, data>::iteOp(c, l, r);
  }

}

#endif
