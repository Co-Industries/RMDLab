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
        Debug{} << "Simulation running";
        QEq();
    }

    void Simulation::GETPARAMS()
    {
        // NULL Transfer Fields not needed in program (used to calc other values)
        Containers::Array<Double> rvdw1, eps, alf, vop, bo131, bo132, bo133;

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

        // Atom dependant variables (line 1)
        Containers::StaticArray<7, Containers::String> _atmname{"C", "H", "O", "N", "S", "Si", "X"};
        arrayAppend(atmname, _atmname);
        Containers::StaticArray<7, Double> _rat{1.3742, 0.6867, 1.3142, 1.245, 1.9647, 2.0276, -0.1};
        arrayAppend(rat, _rat);
        Containers::StaticArray<7, Double> _Val{4.0, 1.0, 2.0, 3.0, 2.0, 4.0, 2.0};
        arrayAppend(Val, _Val);
        Containers::StaticArray<7, Double> _mass{12.0, 1.008, 15.999, 14.0, 32.06, 28.06, 1.0};
        arrayAppend(mass, _mass);
        Containers::StaticArray<7, Double> _rvdw1{1.9684, 1.3525, 1.9741, 1.9951, 2.0783, 2.2042, 2.0};
        arrayAppend(rvdw1, _rvdw1);
        Containers::StaticArray<7, Double> _eps{0.1723, 0.0616, 0.088, 0.1088, 0.2176, 0.1322, 0.0};
        arrayAppend(eps, _eps);
        Containers::StaticArray<7, Double> _gam{0.8712, 0.891, 0.8712, 1.0512, 1.0336, 0.8218, 1.0};
        arrayAppend(gam, _gam);
        Containers::StaticArray<7, Double> _rapt{1.2385, -0.1, 1.1139, 1.1911, 1.5386, 1.5758, -0.1};
        arrayAppend(rapt, _rapt);
        Containers::StaticArray<7, Double> _Vale{4.0, 1.0, 6.0, 5.0, 6.0, 4.0, 6.0};
        arrayAppend(Vale, _Vale);
        // line 2
        Containers::StaticArray<7, Double> _alf{9.406, 9.3858, 10.2186, 9.9303, 9.9676, 11.9413, 10.0};
        arrayAppend(alf, _alf);
        Containers::StaticArray<7, Double> _vop{2.1346, 5.0013, 7.7719, 7.8431, 5.0812, 2.0618, 2.5};
        arrayAppend(vop, _vop);
        Containers::StaticArray<7, Double> _Valboc{4.0, 1.0, 4.0, 4.0, 4.0, 4.0, 4.0};
        arrayAppend(Valboc, _Valboc);
        Containers::StaticArray<7, Double> _povun5{31.0823, 0.0, 29.5271, 32.4758, 35.1648, 11.8211, 0.0};
        arrayAppend(povun5, _povun5);
        Containers::StaticArray<7, Double> _chi{5.7254, 3.8446, 8.5, 6.7768, 6.5, 1.8038, 8.5};
        arrayAppend(chi, _chi);
        Containers::StaticArray<7, Double> _eta{6.9235, 10.0839, 7.1412, 6.8035, 8.2545, 7.3852, 1.5};
        arrayAppend(eta, _eta);
        // line 3
        Containers::StaticArray<7, Double> _vnq{1.2104, -0.1, 0.9909, 1.0636, 1.4703, -1.0, -0.1};
        arrayAppend(vnq, _vnq);
        Containers::StaticArray<7, Double> _plp2{0.0, 0.0, 14.9473, 0.1045, 9.4922, 0.0, 0.0};
        arrayAppend(plp2, _plp2);
        Containers::StaticArray<7, Double> _bo131{5.7419, 3.8461, 9.1371, 2.1604, 8.5146, 6.4918, 8.741};
        arrayAppend(bo131, _bo131);
        Containers::StaticArray<7, Double> _bo132{33.3951, 3.2540, 1.6258, 2.9464, 28.0801, 8.5961, 13.3640};
        arrayAppend(bo132, _bo132);
        Containers::StaticArray<7, Double> _bo133{11.9957, 1.0, 0.1863, 2.5181, 8.501, 0.2368, 0.669};
        arrayAppend(bo133, _bo133);
        // line 4
        Containers::StaticArray<7, Double> _povun2{-2.8983, -15.7683, -3.5965, -4.0959, -10.0773, -3.8112, -11.0};
        arrayAppend(povun2, _povun2);
        Containers::StaticArray<7, Double> _pval3{2.5, 2.1504, 2.5, 2.0047, 2.7466, 3.1873, 2.7466};
        arrayAppend(pval3, _pval3)

            //
            ;
        Containers::StaticArray<7, Double> _Valval{4.0, 1.0, 4.0, 4.0, 6.2298, 4.0, 4.0};
        arrayAppend(Valval, _Valval);
        Containers::StaticArray<7, Double> _pval5{2.9663, 2.8793, 2.9225, 2.8793, 2.8793, 2.5791, 2.8793};
        arrayAppend(pval5, _pval5);

        // Temp
        Debug{} << "[ffield]";
        Debug{} << "line 1------------";
        Debug{} << "<atmname>:" << atmname;
        Debug{} << "<rat>:" << rat;
        Debug{} << "<Val>:" << Val;
        Debug{} << "<mass>:" << mass;
        Debug{} << "(temp) <rvdw1>:" << rvdw1;
        Debug{} << "(temp) <eps>:" << eps;
        Debug{} << "<gam>:" << gam;
        Debug{} << "<rapt>:" << rapt;
        Debug{} << "<Vale>:" << Vale;
        Debug{} << "line 2------------";
        Debug{} << "(temp) <alf>:" << alf;
        Debug{} << "(temp) <vop>:" << vop;
        Debug{} << "<Valboc>:" << Valboc;
        Debug{} << "<povun5>:" << povun5;
        Debug{} << "<chi>:" << chi;
        Debug{} << "<eta>:" << eta;
        Debug{} << "line 3 ------------";
        Debug{} << "<vnq>:" << vnq;
        Debug{} << "<plp2>:" << plp2;
        Debug{} << "(temp) <bo131>:" << bo131;
        Debug{} << "(temp) <bo132>:" << bo132;
        Debug{} << "(temp) <bo133>:" << bo133;
        Debug{} << "line 4 ------------";
        Debug{} << "<povun2>:" << povun2;
        Debug{} << "<pval3>:" << pval3;
        Debug{} << "<Valval>:" << Valval;
        Debug{} << "<pval5>:" << pval5;

        arrayResize(nlpopt, nso);
        arrayResize(Valangle, nso);
        arrayResize(r0s, nso);
        arrayResize(r0p, nso);
        arrayResize(r0pp, nso);

        arrayResize(rvdW, nso);
        arrayResize(Dij, nso);
        arrayResize(alpij, nso);
        arrayResize(gamW, nso);
        arrayResize(gamij, nso);

        // Mass update
        for (UnsignedInt i = 0; i < nso; ++i)
        {
            if (mass[i] < 21.0 && Valboc[i] != Valval[i])
                Valboc[i] = Valval[i];
            nlpopt[i] = 0.5 * (Vale[i] - Val[i]);
            // duplicate values
            Valangle[i] = Valboc[i];

            // calc default r0s, r0p, r0pp:
            arrayResize(r0s[i], nso);
            arrayResize(r0p[i], nso);
            arrayResize(r0pp[i], nso);

            arrayResize(rvdW[i], nso);
            arrayResize(Dij[i], nso);
            arrayResize(alpij[i], nso);
            arrayResize(gamW[i], nso);
            arrayResize(gamij[i], nso);
            for (UnsignedInt j = 0; j < nso; ++j)
            {
                // Terms for the Bond Order calculation:
                r0s[i][j] = 0.5 * (rat[i] + rat[j]);
                r0p[i][j] = 0.5 * (rapt[i] + rapt[j]);
                r0pp[i][j] = 0.5 * (vnq[i] + vnq[j]);
                // Terms used in van der Waals calculation:
                rvdW[i][j] = sqrt(4.0 * rvdw1[i] * rvdw1[j]);
                Dij[i][j] = sqrt(eps[i] * eps[j]);
                alpij[i][j] = sqrt(alf[i] * alf[j]);
                gamW[i][j] = sqrt(vop[i] * vop[j]);
                gamij[i][j] = pow(gam[i] * gam[j], -1.5); // <- gamco in reac.f
            }
        }

        Debug{} << "Loop ---------------";
        Debug{} << "<Valboc>:" << Valboc;
        Debug{} << "<nlpopt>:" << nlpopt;
        Debug{} << "<Valangle>:" << Valangle;
        Debug{} << "<r0s>:" << r0s;
        Debug{} << "<r0p>:" << r0p;
        Debug{} << "<r0pp>:" << r0pp;
        Debug{} << "<rvdW>:" << rvdW;
        Debug{} << "<Dij>:" << Dij;
        Debug{} << "<alpij>:" << alpij;
        Debug{} << "<gamW>:" << gamW;
        Debug{} << "<gamij>:" << gamij;

        // Bonds
        Containers::Array<Containers::Array<Int>> inxn2{nso};
        for (UnsignedInt i = 0; i < nso; ++i)
        {
            arrayResize(inxn2[i], nso);
        }

        // line 1
        Containers::StaticArray<18, Int> typea{1, 1, 2, 1, 3, 1, 3, 4, 2, 2, 1, 2, 3, 4, 5, 6, 2, 3};
        Containers::StaticArray<18, Int> typeb{1, 2, 2, 3, 3, 4, 4, 4, 3, 4, 5, 5, 5, 5, 5, 6, 6, 6};

        Containers::StaticArray<18, Double> _Desig{141.9346, 163.6889, 169.8421, 164.0476, 110.4748, 130.7147, 85.495, 157.7518, 224.3076, 212.1772, 128.7959, 128.6090, 0.0, 0.0, 96.1871, 109.1904, 137.1002, 191.1743};
        Containers::StaticArray<18, Double> _Depi{113.4487, 0.0, 0.0, 117.4881, 155.6441, 175.2276, 114.0081, 67.1322, 0.0, 0.0, 56.4134, 0.0, 0.0, 0.0, 93.7006, 70.8314, 0.0, 52.0733};
        Containers::StaticArray<18, Double> _Depipi{67.6027, 0.0, 0.0, 72.1261, 40.0, 97.2523, 70.1453, 160.9732, 0.0, 0.0, 39.0716, 0.0, 0.0, 0.0, 68.686, 30.0, 0.0, 43.3991};
        Containers::StaticArray<18, Double> _pbe1{0.1554, -0.4525, -0.3591, -0.6031, 0.115, -0.0368, 0.5778, -0.5869, -0.628, -0.3585, 0.0688, -0.5555, 0.5563, 0.4438, 0.0955, 0.2765, -0.1902, -0.2584};
        Containers::StaticArray<18, Double> _pbo5{-0.3045, 0.0, 0.0, -0.1795, -0.1054, -0.4942, -0.107, -0.1824, 0.0, 0.0, -0.4463, 0.0, -0.4038, -0.2034, -0.4781, -0.3, 0.0, -0.3};
        Containers::StaticArray<18, Double> _v13cor{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        Containers::StaticArray<18, Double> _pbo6{30.4515, 6.0, 6.0, 14.9755, 28.5221, 26.7545, 16.6611, 12.0, 6.0, 6.0, 31.1766, 6.0, 49.5611, 40.3399, 17.8574, 16.0, 6.0, 36.0};
        Containers::StaticArray<18, Double> _povun1{0.4283, 0.5921, 0.7503, 0.5413, 0.2, 0.5133, 0.2339, 0.7136, 1.0, 0.3316, 0.4530, 0.4721, 0.6, 0.6, 0.6, 0.1583, 0.4256, 0.8764};
        arrayAppend(Desig, _Desig);
        arrayAppend(Depi, _Depi);
        arrayAppend(Depipi, _Depipi);
        arrayAppend(pbe1, _pbe1);
        arrayAppend(pbo5, _pbo5);
        arrayAppend(v13cor, _v13cor);
        arrayAppend(pbo6, _pbo6);
        arrayAppend(povun1, _povun1);

        Int ih = 0;
        for (UnsignedInt i = 0; i < nboty; ++i)
        {
            ++ih;

            inxn2[typea[i]][typeb[i]] = ih;
            inxn2[typeb[i]][typea[i]] = ih;
        }
        Debug{} << "<inxn2>:" << inxn2;
    }

    void Simulation::INITSYSTEM()
    {
        // Int i, j, k, ity, ist = 0;
        // Containers::StaticArray<3, Int> l;
        // Double mm, gmm, dns;
        // Containers::StaticArray<3, Containers::Array<Double>> mat; // mat(3)
        // Long i8;
        // Double maxrcell;
        // Containers::StaticArray<3, Double> rcsize;

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

        // time unit conversion from [fs] -> time unit
        dt = dt / UTIME;
        // square the spring const in the extended Lagrangian method
        Lex_w2 = 2.0 * Lex_k / dt / dt;
        // get reduced temperature from [K]
        treq = treq / UTEMP0;
        // setup the vector ID and parity for processes, in x, y and z order
        // vID[0] = ;

        arrayResize(dthm, nso);
        arrayResize(hmas, nso);

        for (UnsignedInt i = 0; i < nso; ++i)
        {
            dthm[i] = dt * 0.5 / mass[i];
            hmas[i] = 0.5 * mass[i];
        }
    }

    //! MAIN LOOP
    void Simulation::UPDATE()
    {
        // output
        // if (nstep % pstep == 0)
        //     PRINTE(...);
        // if(nstep % fstep == 0)
        //     OUTPUT(...);

        // if (nstep % sstep == 0)
        //{
        //     for (int i = 0)
        // }
    }

}