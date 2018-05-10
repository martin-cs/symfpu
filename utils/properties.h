/*
** Copyright (C) 2018 Martin Brain
**
** See the file LICENSE for licensing information.
*/

/*
** properties.h
**
** Martin Brain
** martin.brain@cs.ox.ac.uk
** 07/08/14
**
** Macros for specifying invariants in a back-end specific way.
**
** Note that there are two kinds of assertions.
**  - Implementation assertions : should be statically or dynamically
**    resolvable to bools.  Intended to catch cases where the code is
**    buggy or being used incorrectly.
**  - Algorithm assertions : should be trait::prop and are used
**    to document and record properties and assumptions about the
**    floating-point computation.  Depending on the back-end these may
**    be concrete or symbolic and thus handled in different ways.
**
*/

#ifndef SYMFPU_PROPERTIES
#define SYMFPU_PROPERTIES

#define IMPLIES(X,Y) (!(X) || (Y))

#ifndef PRECONDITION
#define PRECONDITION(X) t::precondition(X)
#endif

#ifndef POSTCONDITION
#define POSTCONDITION(X) t::postcondition(X)
#endif

#ifndef INVARIANT
#define INVARIANT(X) t::invariant(X)
#endif

#endif



