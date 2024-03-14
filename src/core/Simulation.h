#ifndef RMD_Simulation_h
#define RMD_Simulation_h

/* A B C D E F G H I J K L M N O P Q R S T U V W X Y Z */
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/Array.h>

#include <Magnum/Math/Range.h>
#include <Magnum/Math/Constants.h>
#include <Magnum/Math/Vector.h>
#include <Magnum/Magnum.h>

namespace Magnum
{
    using namespace Math::Literals;

    class Simulation
    {
    public:
        explicit Simulation();

    protected:
        Containers::Array<UnsignedInt> _atype; /* Atom type [H, He, Na, C, O, ...] */
        Containers::Array<Int> _q;             /* Atom charge */
        Containers::Array<Vector3i> _pos;      /* Atom position [x, y, z] {int} */
        Containers::Array<Vector3i> _v;        /* Atom velocity */
        Containers::Array<Vector3i> f;         /* Atom force */

        // Parameters
        UnsignedInt nso;   // Number of different types of atoms
        UnsignedInt nboty; // Number of different bonds given

        // Atom Dependant
        Containers::Array<std::string> atmname;                                                          /* Chemical Abbrev for each atomtype */
        Containers::Array<Double> Val, Valboc, mass;                                                     /* Valency of atomtype (norm, boc), mass of atomtype */
        Containers::Array<Double> pbo1, pbo2, pbo3, pbo4, pbo5, pbo6;                                    /* Bond Order terms */
        Containers::Array<Double> pboc1, pboc2, pboc3, pboc4, pboc5;                                     /* Bond Order correction terms (f1-5) */
        Containers::Array<Double> v13cor;                                                                /* <kn> */
        Containers::Array<Double> rat, rapt, vnq;                                                        /* r0s/r0p/r0pp for like bonds */
        Containers::Array<Double> ovc;                                                                   /* a flag to apply fn4 and fn5 */
        Containers::Array<Double> Desig, Depi, Depipi;                                                   /* Bond Energy parameters (eq. 6) */
        Containers::Array<Double> pbe1, pbe2;                                                            /* Bond Energy parameters (eq. 6) */
        Containers::Array<Double> Vale;                                                                  /* Lone Pair Energy parameters (eq. 7) */
        Containers::Array<Double> plp1, nlpopt, plp2;                                                    /* Lone Pair Energy parameters (eq.8-10) */
        Containers::Array<Double> povun1, povun2, povun3, povun4, povun5, povun6, povun7, povun8;        /* Over/Undercoordination Energy (eq. 11-12) */
        Containers::Array<Double> pval1, pval2, pval3, pval4, pval5, pval6, pval7, pval8, pval9, pval10; /* Valency Angle Energy (eq. 13a-g) */
        Containers::Array<Double> Valangle, theta00;                                                     //
        Containers::Array<Double> ppen1, ppen2, ppen3, ppen4;                                            /* Penalty Energy (eq. 14ab) */
        Containers::Array<Double> pcoa1, pcoa2, pcoa3, pcoa4;                                            /* Conjugation (3 body) Energy (eq.15) */
        Containers::Array<Double> Valval;                                                                //
        Containers::Array<Double> ptor1, ptor2, ptor3, ptor4;                                            /* Torsional Energy Terms (eq.16abc) */
        Containers::Array<Double> V1, V2, V3;                                                            //
        Containers::Array<Double> pcot1, pcot2;                                                          /* Conjugation (4body) Energy (eq. 17ab) */
        Containers::Array<Double> phb1, phb2, phb3, r0hb;                                                /* Hydrogren Bond Energy (eq. 18) */
        Containers::Array<Containers::Array<Double>> Dij, alpij, rvdW, gamW;                             /* Van der Waals Energy (eq. 21ab) */

        Double pvdW1, pvdW1h, pvdW1inv; //
        Double rchb = 10.0;             /* hydrogen bonding interaction cutoff [A] */
        Double rchb2 = rchb * rchb;     //
        Double Cclmb0 = 332.0638;       /* Coulomb Energy (eq. 22) */
        Double Cclmb0_qeq = 14.4;       //
        Double CEchrge = 23.02;         //
        Double Cclmb = Cclmb0;          //

        Containers::Array<Double> gam;                                                         //
        Containers::Array<Containers::Array<Double>> gamij;                                    //
        Containers::Array<Double> chi, eta;                                                    /* Charge Equilibration part, <chi> electronegativity, <eta> stiffness */
        Containers::Array<Double> bom;                                                         /* Not Understood Parameters */
        Containers::Array<Containers::Array<Double>> r0s, r0p, r0pp;                           /* Bond Order terms (2-atom combo dependant) */
        Containers::Array<Containers::Array<Int>> inxn2;                                       /* 2-body interaction type 1=C-C, 2=H-C, 3=H-H */
        Containers::Array<Containers::Array<Containers::Array<Int>>> inxn3, inxn3hb;           /* 3-body interaction */
        Containers::Array<Containers::Array<Containers::Array<Containers::Array<Int>>>> inxn4; /* 4-body interaction */
        Containers::Array<Double> cBOp1, cBOp3, cBOp5;                                         /* Saved calculations */
        Containers::Array<Double> pbo2h, pbo4h, pbo6h;                                         //
        Containers::Array<Containers::Array<Double>> _switch;                                  /* For debugging purpose variables */
        Containers::Array<Containers::Array<Double>> C_lg, Re_lg;                              /* LG params */
        Containers::Array<Double> rcore2, ecore2, acore2;                                      //
        Containers::Array<Containers::Array<Double>> rcore, ecore, acore;                      //

        // INFO constexpr calculates at compile time (could break)
        // Atoms

        /* For array size statistics
        1-NATOMS
        2-nbrlist
        3-nbrlist for qeq
        4-NBUFFER for move
        5-NBUFFER for copy
        6-NBUFFER for qeq */
        Int nmaxas = 5;
        Containers::StaticArray<5, Containers::Array<Int>> maxas; /* 5 - nmaxas */

        Double lata, latb, latc, lalpha, lbeta, lgamma; /* Lattice parameters */
        Int myid, nprocs, ierr;                         //
        Containers::StaticArray<3, Int> myparity, vID;  //
        Containers::Array<Double> sbuffer, rbuffer;     /* Send and receive buffers */

        Int ns; /* <ns> # of atoms to be sent */
        Int nr; /* <nr> # of atoms to be received */
        Int na; /* <na> # of all of transfered atoms */
        Int ne; /* <ne> # of elements for one atom */
        // Example In case atom type, position and velocity to be sent,  ne = 1+3+3 = 7

        /* constexpr int MODE_COPY = 1, MODE_MOVE = 2, MODE_CPBK = 3, MODE_QCOPY1 = 4, MODE_QCOPY2 = 5; */

        /* <NE_COPY>,<NE_MOVE>,<NE_CPBK> :: Number of Elements to COPY, MOVE atoms and copy back force */
        enum MODE
        {
            MODE_COPY = 1,
            MODE_MOVE,
            MODE_CPBK,
            MODE_QCOPY1,
            MODE_QCOPY2
        };
        Int NE_CPBK = 4;                             //
        Int MAXLAYERS = 5;                           /* MAXimum # of linkedlist cell LAYERS. */
        Int MAXLAYERS_NB = 10;                       //
        Containers::StaticArray<6, Int> target_node; /* stores partner node ID in the 6-communications, if targe_node(i)==-1, the node doesn't have a partner in i-direction. */
        Containers::StaticArray<3, Int> vprocs;      /* For benchmarking, <vprocs> and <mc> will be read from vprocs.in */
        Containers::StaticArray<3, Int> mc;          /* # of unit cells in each directions */
        Containers::Array<Double> rc, rc2;           /* cutoff length for sigma-bonding. */
        Containers::Array<Double> rcpi, rcpp;        /* cutoff length for other bonding. */
        Double MINBOSIG = 1e-3;                      /* criterion to decide <rc> */
        Double MINBO0 = 1e-4;                        /* cutoff bond order */
        Double cutof2_esub = 1e-4, cutof2_bo = 1e-3; //
        Int is_idEh = 1;                             //
        Double cutoff_vpar30;                        /* cutoff_vpar30 = cutof2_bo*vpar30, used in BOPRIM() */

        Int NBUFFER = 30000;               /* Atom buffer / count */
        Int MAXNEIGHBS = 30;               /* Max # of Ngbs one atom may have */
        Int MAXNEIGHBS10 = 1500;           /* Max # of Ngbs within 10[A] */
        Int NMINCELL = 4;                  /* Nr of minimum linkedlist cell <-> minimum grain size */
        Double MAXANGLE = 0.999999999999;  //
        Double MINANGLE = -0.999999999999; //
        Double NSMALL = 1e-10;             //
        Double maxrc;                      /* Max cutoff length. used to decide lcsize */
        Double pi = pi();
        Double sqrtpi_inv = 1.0 / pi;

        /* stress tensor */
        Containers::StaticArray<6, Double> astr;

        /* coefficient of bonding energy derivative */
        Containers::Array<Double> ccbnd, cdbnd;

        Containers::StaticArray<3, Containers::StaticArray<3, Containers::StaticArray<2, Double>>> HH; /* [i][j][k]*/
        Containers::StaticArray<3, Containers::StaticArray<3, Double>> HHi;                            /* [i][j]*/

        Double MDBOX;                            /* MD box */
        Containers::StaticArray<4, Double> LBOX; /* local MD box */
        Containers::StaticArray<3, Double> OBOX; /* origin of box */
        Int NATOMS;                              /* local # of atoms */
        Long GNATOMS;                            /* global # of atoms */
        Int ALLATOMS;
        /* <llist> Linked List */
        /* <header> header atom of linkedlist cell. */
        /* <nacell> Nr of atoms in a likedlist cell. */
        Containers::Array<Int> llist;
        Containers::Array<Containers::Array<Containers::Array<Int>>> header, nacell;

        /* <nbllist> Linked List for non-bonding interaction */
        /* <nbheader> header atom of linkedlist cell for non-bonding interaction */
        /* <nbnacell> Nr of atoms in a likedlist cell for non-bonding interaction */
        Containers::Array<Int> nbllist;
        Containers::Array<Containers::Array<Containers::Array<Int>>> nbheader, nbnacell;

        /* <nbrlist> neighbor list, <nbrindx> neighbor index */
        Containers::Array<Containers::Array<Int>> nbrlist, nbrindx;

        /* <nbplist> neighbor list of nonbonding interaction, non-bonding pair list */
        Containers::Array<Containers::Array<Int>> nbplist;

        /* <BO> Bond Order of atoms i-j (nearest neighb only) - (Eq 3a-3d) */
        Containers::Array<Containers::Array<Containers::Array<Double>>> BO;

        Containers::Array<Double> delta; /* Delta */

        // Output variables from the BOp_CALC() subroutine
        Containers::Array<Containers::Array<Double>> deltap;
        Containers::Array<Containers::Array<Containers::Array<Double>>> dln_BOp;
        Containers::Array<Containers::Array<Double>> dBOp;

        /* For NEW DBO calc */
        Containers::Array<Containers::Array<Double>> exp_delt1, exp_delt2; // exp( -pboc#(inxn) * deltap(i,1) ) - {inxn, i}

        /* A[0123] coefficients for force calculation */
        Containers::Array<Containers::Array<Double>> A0, A1, A2, A3;

        /* Passed between Elnpr and E3body */
        Containers::Array<Double> nlp, dDlp; // Number of Lone Pairs, its derivatives
        Containers::Array<Double> deltalp;

        // Total Energy, Kinetic Energy, Potential Energies
        // 0-Esystem, 1-Ebond, 2-Elp, 3-Eover, 4-Eunder, 5-Eval, 6-Epen
        // 7-Ecoa,  8-Etors, 9-Econj, 10-Ehbond, 11-Evdwaals, 12-Ecoulomb 13-Echarge
        Double TE, KE;                           /* TE: Total Energy,  KE: Kinetic Energy */
        Containers::StaticArray<14, Double> PE;  /* PE :: Potential Energies */
        Double GTE, GKE;                         /* TE: Total Energy,  KE: Kinetic Energy */
        Containers::StaticArray<14, Double> GPE; /* PE :: Potential Energies */

        // Output file format flags
        bool isBinary = false, isBondFile = false, isPDB = false, isXYZ = false;
    };
}

#endif