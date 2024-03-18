/* A B C D E F G H I J K L M N O P Q R S T U V W X Y Z */

#include "Simulation.h"

namespace Magnum
{
    Simulation::Simulation(Int &test) : _test(test)
    {
        Int i, ity, it1, it2, irt;
        Double ctmp;
        Containers::StaticArray<3, Double> dr;

        Debug{} << "Simulation has started";
    }

    void Simulation::ImGuiTest()
    {
        Debug{} << _test;
    }

    void Simulation::INITSYSTEM()
    {
        Int i, j, k, ity, ist = 0;
        Containers::StaticArray<3, Int> l;
        Double mm, gmm, dns;
        Containers::StaticArray<3, Containers::Array<Double>> mat; // mat(3)
        Long i8;
        Double maxrcell;
        Containers::StaticArray<3, Double> rcsize;

        // for QEq (simulation should be Qeq instead of PQeq)
        rctap = rctap0;
    }
}