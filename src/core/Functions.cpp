/* A B C D E F G H I J K L M N O P Q R S T U V W X Y Z */

#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>

#include "Data.h"
#include "Functions.h"

namespace Magnum
{
    void QEq()
    {
        // Method 1 (Original QEq)
        for (UnsignedLong i = 0; i < NATOMS; ++i)
        {
            qsfp[i] = q[i];
            qs[i] = q[i];
        }
        // TODO Implement case(2) Extender Lagrangian method (prob better performance)
    }
}
