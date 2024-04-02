#ifndef RMD_Functions_h
#define RMD_Functions_h

#include <Corrade/Containers/Array.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>

namespace Magnum
{
    using namespace Math::Literals;

    // ? QEq qeq.F90
    void QEq();
    Double dot_product(const Int &mode);

    // ? FORCE pot.F90
    void FORCE();

    // ? bo.F90
    void BOCALC();
    void BOPRIM();
    void BOFULL();
}

#endif