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
        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            atomData[i].qsfp = atomData[i].q;
            atomData[i].qs = atomData[i].q;
        }
        // TODO Implement case(2) Extender Lagrangian method (prob better performance)
    }
}
