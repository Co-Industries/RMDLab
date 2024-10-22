#include "Data.h"

#include <Magnum/Magnum.h>

namespace Magnum
{
    Double dt = 0.25;
    Containers::Array<AtomData> atomData;
    Containers::StaticArray<3, Atom> atom;
    Containers::StaticArray<3, Bond> bond;
    Containers::StaticArray<1, H_Bond> h_bond;
    Containers::StaticArray<6, Angle> angle;
    Containers::StaticArray<6, Torsion> torsion;

    // * Parameters
    std::size_t NATOMS = 0;
    Float atomRadius = 0.1;
    Float randomVelocity = 0.02;
    bool drawOctreeBounds = true;

    // * Simulation constants
    const Float atomRange = 0.1;

    const std::size_t NMAXQEq = 500; // ! from rxmd.in
    const std::size_t NTABLE = 5000;
    const Float rctap0 = 10.0;
    const Double MINBOSIG = 1.0e-3;
    const Double MINBO0 = 1.0e-4;
    const Double cutof2_esub = 1.0e-4;
    const Double cutof2_bo = 1.0e-3;
    const Double vpar30 = 0.1;
    const Double vpar1 = 50.0, vpar2 = 9.5469;
    const Double QEq_tol = 1.0e-7;
    const Double pvdW1 = 1.5591;
    Double pvdW1h, pvdW1inv;
    
    const Double MAXANGLE = 0.999999999999;
    const Double MINANGLE = -0.999999999999;
    const Double NSMALL = 1.0e-10;

    // ? Coulomb Energy (eq. 22)
    const Double Cclmb0_qeq = 14.4;
    const Double Cclmb0 = 332.0638;
    const Double CEchrge = 23.02;
    Double rchb = 10.0;
    Double rchb2 = rchb * rchb;

    Float rctap, rctap2;
    Float UDR, UDRi;
    Containers::StaticArray<5001, Containers::Array<Double>> TBL_Eclmb_QEq;
    Containers::StaticArray<5001, Containers::Array<Double>> TBL_Evdw_p, TBL_Eclmb_p;
    Containers::StaticArray<5001, Containers::Array<Double>> TBL_Evdw_d, TBL_Eclmb_d;
    Containers::StaticArray<8, Double> CTap;
    Double cutoff_vpar30;

    std::size_t nso = 3;
    std::size_t nboty = 3;
    std::size_t nhbty = 1;
    std::size_t nvaty = 6;
    std::size_t ntoty = 6;

    // ? QEq
    Containers::StaticArray<2, Double> Gnew;

    // ? BOPRIM
    Containers::StaticArray<3, Double> arg_BOpij;

    // ? ENbond
    Containers::StaticArray<14, Double> PE;

    //? Elnpr
    Containers::StaticArray<7, Double> CEover;
    Containers::StaticArray<6, Double> CEunder;

    //? Force
    Containers::StaticArray<6, Double> astr;

    const Double UTIME = 1.0e3 / 20.455;
    Double Lex_w2 = 1.0;
    Double Lex_k = 2.0;
    std::size_t nstep = 0, qstep = 1;

    Double BORDER = 100.0;
    Double BORDER2 = 2 * BORDER;

    Double BO_MAX = 1.0;
}