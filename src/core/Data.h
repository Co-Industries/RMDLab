#ifndef RMD_Data_h
#define RMD_Data_h

#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/StaticArray.h>
#include <Corrade/Containers/String.h>

#include <Magnum/Math/Color.h>
#include <Magnum/Math/Constants.h>
#include <Magnum/Math/Range.h>
#include <Magnum/Math/Vector.h>
#include <Magnum/Math/Vector3.h>

#include <Magnum/Magnum.h>

namespace Magnum
{
    extern UnsignedLong NATOMS; // Number of Atoms
    extern Containers::Array<Double> qsfp, qsfv, qtfp, qtfv, qs;
    extern Containers::Array<Double> q; // Atom charge
    extern Int NMAXQEq;                 // Number of MAXimum iteration in QEq routine
    extern const Double rctap0;         // [A]
}

#endif