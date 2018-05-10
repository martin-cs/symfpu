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



