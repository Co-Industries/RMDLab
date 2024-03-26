#include "Data.h"

#include <Magnum/Magnum.h>

namespace Magnum
{
    UnsignedLong NATOMS;
    Int NMAXQEq = 500; // ! from rxmd.in
    Containers::Array<Double> qsfp, qsfv, qtfp, qtfv, qs;
    Containers::Array<Double> q;
    const Double rctap0 = 10.0;
}