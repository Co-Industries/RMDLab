/* A B C D E F G H I J K L M N O P Q R S T U V W X Y Z */

#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>

#include "Data.h"
#include "Functions.h"

namespace Magnum
{
    Double dot_product(const Int &mode)
    {
        Double result = 0.0;
        switch (mode)
        {
        case 1:
            for (std::size_t i = 0; i < NATOMS; ++i)
            {
                result += atomData[i].gs * atomData[i].gs;
            }
            return result;
        case 2:
            for (std::size_t i = 0; i < NATOMS; ++i)
            {
                result += atomData[i].gt * atomData[i].gt;
            }
            return result;
        default:
            return result;
        }
        return result;
    }

    void QEq()
    {
        Double gssum = 0.0;
        Double gtsum = 0.0;

        // Method 1 (Original QEq) t
        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            atomData[i].qsfp = atomData[i].q;
            atomData[i].qsfv = 0.0;
            atomData[i].qs = atomData[i].q;
            atomData[i].qt = 0.0;
            // TODO Implement case(2) Extender Lagrangian method (prob better performance) only works for a few atoms i guess

            // ? get_gradient [Gnew]

            for (std::size_t j0 = 0; j0 < atomData[i].neighbors.size(); ++j0)
            {
                std::size_t j = atomData[i].neighbors[j0];
                gssum += atomData[i].hessian[j0] * atomData[j].qs;
                gtsum += atomData[i].hessian[j0] * atomData[j].qt;
            }
            std::size_t _type = atomData[i].type;
            Double _eta = atom[_type].eta;
            atomData[i].gs = atom[_type].chi - _eta * atomData[i].qs - gssum;
            atomData[i].gt = 1.0 - _eta * atomData[i].qt - gtsum;
        }

        Gnew[0] = dot_product(1);
        Gnew[1] = dot_product(2);
    }

    void BOPRIM()
    {
        //for (std::size_t i = 0; i < nso; ++i)
        //{
        //     atomData[i].
        //}
    }

    void BOCALC()
    {
        BOPRIM();
    }

    void FORCE()
    {
        BOCALC();
    }
}
