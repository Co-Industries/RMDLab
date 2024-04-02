#ifndef RMD_Data_h
#define RMD_Data_h

#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/StaticArray.h>
#include <Corrade/Containers/String.h>
#include <Corrade/Utility/Arguments.h> //std::string

#include <Magnum/Math/Color.h>
#include <Magnum/Math/Constants.h>
#include <Magnum/Math/Range.h>
#include <Magnum/Math/Vector.h>
#include <Magnum/Math/Vector3.h>

#include <Magnum/Magnum.h>

namespace Magnum
{
    struct AtomData
    {
        std::size_t type;
        Vector3d position, velocity;
        Double q, qs, qt, gs, gt, qsfp, qsfv, deltap0; // deltap1 instead of deltap[0]
        Double delta;
        Containers::Array<Double> hessian;
        Containers::Array<std::size_t> neighbors;
        Containers::Array<std::size_t> bonds;
        Containers::Array<Vector3d> bo, dln_BOp;
        Containers::Array<Double> dBOp, bo_sum, A0, A1, A2, A3;
        Containers::StaticArray<3, Double> deltap;
    };
    extern Containers::Array<AtomData> atomData;

    struct Atom
    {
        std::string name;
        Double vop, gam, eta, chi, rat, rapt, vnq, Val, Valval;
        Containers::Array<std::size_t> inxn2;
        Containers::Array<Double> gamW, gamij, r0s, r0p, r0pp;
        Double bo131, bo132, bo133;
    };
    extern Containers::StaticArray<3, Atom> atom; // 3 = nso

    struct Bond
    {
        Double rc, rc2;                            /* <RCUT>: cutoff length for sigma-bonding */
        Double pbo1, pbo2, pbo3, pbo4, pbo5, pbo6; /* Bond Order terms */
        Double pbo2h, pbo4h, pbo6h;                /* Half of pbo */
        Double pboc1, pboc2, pboc3, pboc4, pboc5;  /* Bond Order correction terms (f1-5) */
        Double cBOp1, cBOp3, cBOp5;                /* Bond order calculations */
        Double ovc, v13cor;                        /* a flag to apply fn4 and fn5 */
        Containers::StaticArray<3, Double> swh;    /* <switch> flag to omit pi and double pi bond in bond-order prime calculation */
    };
    extern Containers::StaticArray<3, Bond> bond; // 3 = nboty

    // * Parameters
    extern std::size_t NATOMS;   // Number of Atoms
    extern Float atomRadius;     // rendering
    extern Float randomVelocity; // ! not needed
    extern bool drawOctreeBounds;

    // * Simulation constants
    extern const Float atomRange; // [rctap0]
    extern const Int NMAXQEq;     // Number of MAXimum iteration in QEq routine
    extern const Float rctap0;    // [10A]
    extern const std::size_t NTABLE;
    extern const Double MINBOSIG; /* <minBOsig>: criterion to decide <rc> */
    extern const Double cutof2_bo;
    extern const Double vpar30;       /* Cutoff for bond order (*100) */
    extern const Double vpar1, vpar2; /* Overcoordination parameter [ffield]*/

    // ? Coulomb Energy (eq. 22)
    extern const Double Cclmb0_qeq;

    extern Float rctap, rctap2;
    extern Float UDR, UDRi;
    extern Containers::StaticArray<5000, Containers::Array<Double>> TBL_Eclmb_QEq; // 5000 = NTABLE
    extern Containers::StaticArray<8, Double> CTap;
    extern Double cutoff_vpar30; /* cutoff_vpar30 = cutof2_bo*vpar30, used in BOPRIM() */

    extern std::size_t nso;
    extern std::size_t nboty;

    //? QEq
    extern Containers::StaticArray<2, Double> Gnew;

    //? BOPRIM
    extern Containers::StaticArray<3, Double> arg_BOpij;
}

#endif