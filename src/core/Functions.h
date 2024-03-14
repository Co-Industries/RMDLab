#ifndef RMD_Functions_h
#define RMD_Functions_h

#include <Corrade/Containers/Array.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Magnum.h>

namespace Magnum
{
    using namespace Math::Literals;

    class Functions
    {
    public:
        explicit Functions(Containers::Array<UnsignedInt> &atype, Containers::Array<Vector3i> &pos);

    private:
        /* Calculate atom bond order*/
        void BOCALC();
        /* [BO prime] Calculates the BOp(0:3,i,j) and the deltap(i) */
        void BOPRIM();
        /* [full BO] Subroutine calculates the Bond Order and its derivatives */
        void BOFULL();

    protected:
        Containers::Array<UnsignedInt> &_atype; /* Atom type [H, He, Na, C, O, ...]*/
        Containers::Array<Vector3i> &_pos;      /* Atom position [x, y, z] {int}*/
    };
}

#endif