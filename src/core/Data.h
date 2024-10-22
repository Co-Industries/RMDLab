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
    extern Double dt; /* velocity scaling factor, <dt> one time step */
    struct AtomData
    {
        std::size_t type;
        Vector3d position, velocity, force;
        Double q, qs, qt, qsfp, qsfv, deltap0; // deltap1 instead of deltap[0]
        Double gs, gt;                         /* gradient */
        Double hs, ht;                         /* conjugate direction */
        Double delta;
        Double ccbnd, cdbnd; /* coefficient of bonding energy derivative  */

        Double nlp, dDlp; /* Number of Lone Pairs, its derivatives */
        Double deltalp;

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
        Float size;
        Color3 color;
        Double vop, gam, eta, chi, rat, rapt, vnq, Val, Valangle, Valboc;
        Double Valval, eps, alf, rvdw1;
        // ? Elnpr()
        Double Vale, plp1, plp2, nlpopt, mass;
        Double povun2, povun3, povun4, povun5, povun6, povun7, povun8; /* Over / under coordination Energy */
        // ? angle
        Double pval3, pval5;
        
        Containers::Array<std::size_t> inxn2;
        Containers::Array<Containers::Array<std::size_t>> inxn3;
        Containers::Array<Containers::Array<std::size_t>> inxn3hb;
        Containers::Array<Containers::Array<Containers::Array<std::size_t>>> inxn4;
        
        Containers::Array<Double> gamW, gamij, r0s, r0p, r0pp, Dij, alpij, rvdW;
        Double bo131, bo132, bo133;
        std::size_t natoms_per_type;
        
        // ? vkick
        Double dthm, hmas;
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
        Double pbe1, pbe2, Desig, Depi, Depipi;    /* Bond Energy parameters (eq. 6) */
        Double povun1;     /* Overcoordination Energy (eq. 11) */
        Containers::StaticArray<3, Double> swh;    /* <switch> flag to omit pi and double pi bond in bond-order prime calculation */
    };
    extern Containers::StaticArray<3, Bond> bond; // 3 = nboty

    struct H_Bond
    {
        Double r0hb, phb1, phb2, phb3; /* Hydrogren Bond Energy (eq. 18) */
    };
    extern Containers::StaticArray<1, H_Bond> h_bond; // 1 = nhbty

    struct Angle
    {
        Double pval1, pval2, pval3, pval4, pval6, pval7, pval8, pval9, pval10; // Valency Angle Energy (eq. 13a-g)
        Double theta00;
        Double ppen1, ppen2, ppen3, ppen4; // Penalty Energy (eq. 14ab)
        Double pcoa1, pcoa2, pcoa3, pcoa4; // Conjugation (3 body) Energy (eq.15)
    };
    extern Containers::StaticArray<6, Angle> angle; // 6 = nvaty;

    struct Torsion
    {
        Double ptor1, ptor2, ptor3, ptor4, V1, V2, V3;
        Double pcot1, pcot2;
    };
    extern Containers::StaticArray<6, Torsion> torsion;
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
    extern const Double MINBO0; /* <minBO0>: cutoff bond order  */
    extern const Double cutof2_esub;
    extern const Double cutof2_bo;
    extern const Double vpar30;       /* Cutoff for bond order (*100) */
    extern const Double vpar1, vpar2; /* Overcoordination parameter [ffield]*/
    extern const Double QEq_tol;      /* <QEq_thrsld> energy criterion in QEq routine */
    extern const Double pvdW1;
    extern Double pvdW1h, pvdW1inv;

    extern const Double MAXANGLE, MINANGLE, NSMALL;

    // ? Coulomb Energy (eq. 22)
    extern const Double Cclmb0_qeq; /* [ev] */
    extern const Double Cclmb0;     // [kcal/mol/A] line 2481 in poten.f
    extern const Double CEchrge;    /* [ev] */
    extern Double rchb, rchb2;

    extern Float rctap, rctap2;
    extern Float UDR, UDRi;
    extern Containers::StaticArray<5001, Containers::Array<Double>> TBL_Eclmb_QEq;           // 5000 = NTABLE (+1 for itb + 1 = 0)
    extern Containers::StaticArray<5001, Containers::Array<Double>> TBL_Evdw_p, TBL_Eclmb_p; // 0: potencial
    extern Containers::StaticArray<5001, Containers::Array<Double>> TBL_Evdw_d, TBL_Eclmb_d; // 1: derivative of potencial

    extern Containers::StaticArray<8, Double> CTap;
    extern Double cutoff_vpar30; /* cutoff_vpar30 = cutof2_bo*vpar30, used in BOPRIM() */

    extern std::size_t nso;
    extern std::size_t nboty;
    extern std::size_t nhbty;
    extern std::size_t nvaty;
    extern std::size_t ntoty;

    //? QEq
    extern Containers::StaticArray<2, Double> Gnew;

    //? BOPRIM
    extern Containers::StaticArray<3, Double> arg_BOpij;

    //? ENbond
    // 0-Esystem, 1-Ebond, 2-Elp, 3-Eover, 4-Eunder, 5-Eval, 6-Epen
    // 7-Ecoa,  8-Etors, 9-Econj, 10-Ehbond, 11-Evdwaals, 12-Ecoulomb 13-Echarge
    extern Containers::StaticArray<14, Double> PE; /* Potential Energies */

    //? Elnpr
    extern Containers::StaticArray<7, Double> CEover;
    extern Containers::StaticArray<6, Double> CEunder;

    //? FORCE
    extern Containers::StaticArray<6, Double> astr;

    extern const Double UTIME; /* 1 = 1/20.445[ps] = 48.88780[fs] */
    extern Double Lex_w2, Lex_k;
    extern std::size_t nstep, qstep;

    extern Double BORDER;
    extern Double BORDER2;

    extern Double BO_MAX;
}

#endif