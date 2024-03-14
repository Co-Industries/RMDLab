/* A B C D E F G H I J K L M N O P Q R S T U V W X Y Z */

#include "Functions.h"

namespace Magnum
{
    Functions::Functions(Containers::Array<UnsignedInt> &atype,
                         Containers::Array<Vector3i> &pos) : _atype(atype),
                                                             _pos(pos) {}
    void Functions::BOCALC()
    {
        BOPRIM();
        BOFULL();
    }

    void Functions::BOPRIM()
    {
        Int n, i, j, j1, i1, ity, jty, inxn, c1, c2, c3;
        UnsignedInt dr2;
        Containers::StaticArray<3, UnsignedInt> dr;
        Containers::StaticArray<3, UnsignedInt> arg_BOpij;

        /* initialize deltap(1:,1) to -Val(atype(i)) */
        
    }
}