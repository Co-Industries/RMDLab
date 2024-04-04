#ifndef RMD_Functions_h
#define RMD_Functions_h

#include <Corrade/Containers/Array.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>

namespace Magnum
{
    using namespace Math::Literals;

    // ? QEq qeq.F90
    void QEq(); /* Two vector electronegativity equilization routine */
    Double dot_product(const Int &mode);

    // ? FORCE pot.F90
    void FORCE();
    void ENbond(); /* Calculates the energy and the forces due to the Van der Waals and Coulomb terms */

    // ? bo.F90
    void BOCALC();
    void BOPRIM(); /* Calculates the BOp(0:3,i,j) and the deltap(i) */
    void BOFULL(); /* Calculates the Bond Order and its derivatives */
}

#endif