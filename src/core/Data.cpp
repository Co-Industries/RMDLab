#include "Data.h"

#include <Magnum/Magnum.h>

namespace Magnum
{
    // * Parameters
    std::size_t NATOMS = 0;
    Float atomRadius = 0.1;
    Float randomVelocity = 0.02;
    bool drawOctreeBounds = true;

    // * Simulation constants
    const Float atomRange = 0.1;

    const Int NMAXQEq = 500; // ! from rxmd.in
    const std::size_t NTABLE = 5000;
    const Float rctap0 = 10.0;
    Float rctap, rctap2;
    Float UDR, UDRi;
    Containers::StaticArray<8, Double> CTap;
    Containers::Array<AtomData> atomData;
}