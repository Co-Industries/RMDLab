#ifndef RMD_Data_h
#define RMD_Data_h

#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/StaticArray.h>
#include <Corrade/Containers/String.h>
#include <Corrade/Utility/Arguments.h> //std::string

#include <Magnum/Math/Color.h>
#include <Magnum/Math/Constants.h>
#include <Magnum/Math/Range.h>
#include <Magnum/Math/Vector.h>
#include <Magnum/Math/Vector3.h>

#include <Magnum/Magnum.h>

namespace Magnum
{
    // * Parameters
    extern std::size_t NATOMS;   // Number of Atoms
    extern Float atomRadius;     // rendering
    extern Float randomVelocity; // ! not needed
    extern bool drawOctreeBounds;

    // * Simulation constants
    extern const Float atomRange; // [rctap0]
    extern const Int NMAXQEq;     // Number of MAXimum iteration in QEq routine
    extern const Float rctap0;    // [10A]
    extern Float rctap, rctap2;
    extern Float UDR, UDRi;
    extern const std::size_t NTABLE;
    extern Containers::StaticArray<8, Double> CTap;

    struct AtomData
    {
        std::string name;
        std::size_t type;
        Vector3d position, velocity;
        Double q, qs, qsfp, qsfv;
        Containers::Array<Double> hessian;
        // Containers::StaticArray<5000, Double>
    };
    extern Containers::Array<AtomData> atomData;
}

#endif