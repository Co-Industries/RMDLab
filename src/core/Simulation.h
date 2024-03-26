#ifndef RMD_Simulation_h
#define RMD_Simulation_h

#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/StaticArray.h>
#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/String.h>

#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Math/Constants.h>
#include <Magnum/Math/Range.h>
#include <Magnum/Math/Vector.h>
#include <Magnum/Math/Vector3.h>

namespace Magnum
{
    using namespace Math::Literals;
    // class Atom
    //{
    // protected:
    //     Containers::String name;
    //     Double mass;
    //     Double q;
    //     Color3 color;
    // };

    // TODO CODE CLEANUP
    class Simulation
    {
    public:
        // TODO Make these parameters into a struct
        explicit Simulation();
        void UPDATE();
        void RUN();

    private:
        void GETPARAMS();
        void INITSYSTEM();
        void ALLOCATE();
    };
}

#endif