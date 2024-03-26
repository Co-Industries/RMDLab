/* A B C D E F G H I J K L M N O P Q R S T U V W X Y Z */

#include <Magnum/Math/Functions.h>

#include "Data.h"
#include "Functions.h"
#include "Simulation.h"

namespace Magnum
{
    // TODO CODE CLEANUP (struct)
    Simulation::Simulation()
    {
        NATOMS = 100;
        Debug{} << "Simulation has loaded";
    }

    void Simulation::ALLOCATE()
    {
        arrayResize(qsfp, NATOMS, 0.0);
        arrayResize(qsfv, NATOMS, 0.0);
        arrayResize(q, NATOMS, 0.0);
        arrayResize(qs, NATOMS);
    }

    void Simulation::RUN()
    {
        // Int i, ity, it1, it2, irt;
        // Double ctmp;
        Containers::StaticArray<3, Double> dr;

        // GETPARAMS();
        // INITSYSTEM();

        // Debug{} << "Parameters -------";
        // Debug{} << "<pvdW1>:" << pvdW1;
        // Debug{} << "<cutoff_vpar30>:" << cutoff_vpar30;
        // Debug{} << "<nso>:" << nso;
        // Debug{} << "<nboty>:" << nboty;
        // Debug{} << "<plp1>:" << plp1;
        // Debug{} << "<povun3>:" << povun3;
        // Debug{} << "<povun4>:" << povun4;
        // Debug{} << "<povun6>:" << povun6;
        // Debug{} << "<povun7>:" << povun7;
        // Debug{} << "<povun8>:" << povun8;
        ALLOCATE();
        Debug{} << "Simulation running";
        QEq();
    }
}