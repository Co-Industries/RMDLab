/* A B C D E F G H I J K L M N O P Q R S T U V W X Y Z */

#include "Simulation.h"
#include <Magnum/Math/Functions.h>

namespace Magnum
{
    Simulation::Simulation(
        Double &_pvdW1,
        Double &_cutoff_vpar30,
        UnsignedInt &_nso,
        Double &_plp1param,
        Double &_povun3param,
        Double &_povun4param,
        Double &_povun6param,
        Double &_povun7param,
        Double &_povun8param) : pvdW1(_pvdW1),
                                cutoff_vpar30(_cutoff_vpar30),
                                nso(_nso),
                                plp1param(_plp1param),
                                povun3param(_povun3param),
                                povun4param(_povun4param),
                                povun6param(_povun6param),
                                povun7param(_povun7param),
                                povun8param(_povun8param)
    {
        Debug{} << "Simulation has loaded";
    }

    void Simulation::run()
    {
        Debug{} << "Simulation running";
        Int i, ity, it1, it2, irt;
        Double ctmp;
        Containers::StaticArray<3, Double> dr;

        GETPARAMS();
        INITSYSTEM();

        Debug{} << "<pvdW1>:" << pvdW1;
        Debug{} << "<cutoff_vpar30>:" << cutoff_vpar30;
        Debug{} << "<nso>:" << nso;
        Debug{} << "<plp1>:" << plp1;
        Debug{} << "<povun3>:" << povun3;
        Debug{} << "<povun4>:" << povun4;
        Debug{} << "<povun6>:" << povun6;
        Debug{} << "<povun7>:" << povun7;
        Debug{} << "<povun8>:" << povun8;
    }

    void Simulation::GETPARAMS()
    {
        // vdWaals
        pvdW1h = 0.5 * pvdW1;
        pvdW1inv = 1.0 / pvdW1;

        // Variable allocation with values
        arrayResize(plp1, nso, plp1param);
        arrayResize(povun3, nso, povun3param);
        arrayResize(povun4, nso, povun4param);
        arrayResize(povun6, nso, povun6param);
        arrayResize(povun7, nso, povun7param);
        arrayResize(povun8, nso, povun8param);
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
        rctap2 = pow(rctap, 2);

        // CTap
        CTap[0] = 1.0;
        CTap[1] = 0.0;
        CTap[2] = 0.0;
        CTap[3] = 0.0;
        CTap[4] = -35.0 / pow(rctap, 4);
        CTap[5] = 84.0 / pow(rctap, 5);
        CTap[6] = -70.0 / pow(rctap, 6);
        CTap[7] = 20.0 / pow(rctap, 7);
    }
}