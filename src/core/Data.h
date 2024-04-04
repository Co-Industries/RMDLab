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
        Vector3d position, velocity, force;
        Double q, qs, qt, qsfp, qsfv, deltap0; // deltap1 instead of deltap[0]
        Double gs, gt;                         /* gradient */
        Double hs, ht;                         /* conjugate direction */
        Double delta;
        Double ccbnd; /* coefficient of bonding energy derivative  */
        Containers::Array<Double> hessian, dpq2;
        Containers::Array<std::size_t> neighbors, bonds;
        Containers::Array<Vector3d> bo, BO, dln_BOp;
        Containers::Array<Double> dBOp, bo_sum, BO_sum, A0, A1, A2, A3;
        Containers::StaticArray<3, Double> deltap;
    };
    extern Containers::Array<AtomData> atomData;

    struct Atom
    {
        std::string name;
        Color3 color;
        Double vop, gam, eta, chi, rat, rapt, vnq, Val, Valval, eps, alf, rvdw1;
        Containers::Array<std::size_t> inxn2;
        Containers::Array<Double> gamW, gamij, r0s, r0p, r0pp, Dij, alpij, rvdW;
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
        Double pbe1, pbe2, Desig, Depi, Depipi;                  /* Bond Energy parameters (eq. 6) */
        Containers::StaticArray<3, Double> swh;    /* <switch> flag to omit pi and double pi bond in bond-order prime calculation */
    };
    extern Containers::StaticArray<3, Bond> bond; // 3 = nboty

    // * Parameters
    extern std::size_t NATOMS;   // Number of Atoms
    extern Float atomRadius;     // rendering
    extern Float randomVelocity; // ! not needed
    extern bool drawOctreeBounds;

    // * Simulation constants
    extern const Float atomRange;     // [rctap0]
    extern const std::size_t NMAXQEq; // Number of MAXimum iteration in QEq routine
    extern const Float rctap0;        // [10A]
    extern const std::size_t NTABLE;
    extern const Double MINBOSIG; /* <minBOsig>: criterion to decide <rc> */
    extern const Double cutof2_bo;
    extern const Double vpar30;       /* Cutoff for bond order (*100) */
    extern const Double vpar1, vpar2; /* Overcoordination parameter [ffield]*/
    extern const Double QEq_tol;      /* <QEq_thrsld> energy criterion in QEq routine */
    extern const Double pvdW1;
    extern Double pvdW1h, pvdW1inv;

    // ? Coulomb Energy (eq. 22)
    extern const Double Cclmb0_qeq; /* [ev] */
    extern const Double Cclmb0;     // [kcal/mol/A] line 2481 in poten.f
    extern const Double CEchrge;    /* [ev] */

    extern Float rctap, rctap2;
    extern Float UDR, UDRi;
    extern Containers::StaticArray<5001, Containers::Array<Double>> TBL_Eclmb_QEq;           // 5000 = NTABLE (+1 for itb + 1 = 0)
    extern Containers::StaticArray<5001, Containers::Array<Double>> TBL_Evdw_p, TBL_Eclmb_p; // 0: potencial
    extern Containers::StaticArray<5001, Containers::Array<Double>> TBL_Evdw_d, TBL_Eclmb_d; // 1: derivative of potencial

    extern Containers::StaticArray<8, Double> CTap;
    extern Double cutoff_vpar30; /* cutoff_vpar30 = cutof2_bo*vpar30, used in BOPRIM() */

    extern std::size_t nso;
    extern std::size_t nboty;

    //? QEq
    extern Containers::StaticArray<2, Double> Gnew;

    //? BOPRIM
    extern Containers::StaticArray<3, Double> arg_BOpij;

    //? ENbond
    // 0-Esystem, 1-Ebond, 2-Elp, 3-Eover, 4-Eunder, 5-Eval, 6-Epen
    // 7-Ecoa,  8-Etors, 9-Econj, 10-Ehbond, 11-Evdwaals, 12-Ecoulomb 13-Echarge
    extern Containers::StaticArray<14, Double> PE; /* Potential Energies */
}

#endif