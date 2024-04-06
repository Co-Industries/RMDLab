/* A B C D E F G H I J K L M N O P Q R S T U V W X Y Z */

#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Math/Constants.h>

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
        case 3:
            for (std::size_t i = 0; i < NATOMS; ++i)
            {
                result += atomData[i].gs * atomData[i].hs;
            }
            return result;
        case 4:
            for (std::size_t i = 0; i < NATOMS; ++i)
            {
                result += atomData[i].gt * atomData[i].ht;
            }
            return result;
        default:
            return result;
        }
        return result;
    }

    void QEq()
    {
        Double gssum, gtsum, drtb, _eta;
        std::size_t j, itb, itb1, inxn, itype;

        Double GEst1, GEst2, Est, Est1, _hshs, _hsht, hshs_sum, hsht_sum;
        Vector2d g_h, h_hsh, lmin;

        Double ssum, tsum, mu;

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            // Method 1 (Original QEq)
            atomData[i].qsfp = atomData[i].q;
            atomData[i].qsfv = 0.0;
            atomData[i].qs = atomData[i].q;
            atomData[i].qt = 0.0;
            // TODO Implement case(2) Extender Lagrangian method (prob better performance) only works for a few atoms i guess

            for (std::size_t n = 0; n < atomData[i].neighbors.size(); ++n)
            {
                j = atomData[i].neighbors[n];

                // ? qeq_initialize()
                // if (j < i)
                //    break;

                itb = Int(atomData[i].dpq2[n] * Double(UDRi));
                itb1 = itb + 1;
                drtb = atomData[i].dpq2[n] - itb * Double(UDR);
                drtb = drtb * Double(UDRi);

                inxn = atom[atomData[i].type].inxn2[atomData[j].type] - 1;
                arrayAppend(atomData[i].hessian, InPlaceInit, (1.0 - drtb) * TBL_Eclmb_QEq[itb][inxn] + drtb * TBL_Eclmb_QEq[itb1][inxn]);
                // arrayAppend(atomData[j].hessian, InPlaceInit, (1.0 - drtb) * TBL_Eclmb_QEq[itb][inxn] + drtb * TBL_Eclmb_QEq[itb1][inxn]);

                // ? get_gradient [Gnew]

                gssum += atomData[i].hessian[n] * atomData[j].qs;
                gtsum += atomData[i].hessian[n] * atomData[j].qt;

                // TODO might need to do atomData[j] - idk
            }
            itype = atomData[i].type;
            _eta = atom[itype].eta;
            atomData[i].gs = atom[itype].chi - _eta * atomData[i].qs - gssum;
            atomData[i].gt = 1.0 - _eta * atomData[i].qt - gtsum;

            atomData[i].hs = atomData[i].gs;
            atomData[i].ht = atomData[i].gt;
        }

        Gnew[0] = dot_product(1);
        Gnew[1] = dot_product(2);

        GEst2 = 1.0e99;

        for (std::size_t nstep_qeq = 0; nstep_qeq < NMAXQEq; ++nstep_qeq)
        {
            // ? get_hsh()
            Est = 0.0, hshs_sum = 0.0, hsht_sum = 0.0;

            for (std::size_t i = 0; i < NATOMS; ++i)
            {
                itype = atomData[i].type;
                _eta = atom[itype].eta;
                _hshs = _eta * atomData[i].hs;
                _hsht = _eta * atomData[i].ht;
                Est += atom[itype].chi * atomData[i].q + 0.5 * _eta * atomData[i].q * atomData[i].q;

                for (std::size_t n = 0; n < atomData[i].neighbors.size(); ++n)
                {
                    j = atomData[i].neighbors[n];

                    // if (j < i)
                    //     break;

                    _hshs += atomData[i].hessian[n] * atomData[j].hs;
                    _hsht += atomData[i].hessian[n] * atomData[j].ht;
                    // get half of potential energy, then sum it up if atoms are resident
                    // Double Est1 = 0.5 * atomData[i].hessian[n] * atomData[i].q * atomData[j].q;
                    Est1 = atomData[i].hessian[n] * atomData[i].q * atomData[j].q;
                    Est += Est1;
                    // Est += Est1;
                }

                hshs_sum += _hshs * atomData[i].hs;
                hsht_sum += _hsht * atomData[i].ht;
            }

            // ? QEq
            GEst1 = Est;
            if (0.5 * (abs(GEst2) + abs(GEst1)) < QEq_tol)
                break;
            if (abs(GEst2) > 0.0 && (abs(GEst1 / GEst2 - 1.0) < QEq_tol))
                break;
            GEst2 = GEst1;
            // line minimization factor of <s> vector
            g_h[0] = dot_product(3);
            h_hsh[0] = hshs_sum;

            // line minimization factor of <t> vector
            g_h[1] = dot_product(4);
            h_hsh[1] = hsht_sum;

            // TODO MPI_ALLREDUCE does some buffer stuff, might be important

            lmin = g_h / h_hsh;
            ssum = 0.0, tsum = 0.0;

            for (std::size_t i = 0; i < NATOMS; ++i)
            {
                atomData[i].qs += lmin[0] * atomData[i].hs;
                atomData[i].qt += lmin[1] * atomData[i].ht;
                ssum += atomData[i].qs;
                tsum += atomData[i].qt;
            }

            mu = ssum / tsum;

            gssum = 0.0;
            gtsum = 0.0;

            for (std::size_t i = 0; i < NATOMS; ++i)
            {
                atomData[i].q = atomData[i].qs - mu * atomData[i].qt;
                for (std::size_t n = 0; n < atomData[i].neighbors.size(); ++n)
                {
                    j = atomData[i].neighbors[n];

                    // ? get_gradient [Gnew]
                    // if (j < i)
                    //    break;

                    gssum += atomData[i].hessian[n] * atomData[j].qs;
                    gtsum += atomData[i].hessian[n] * atomData[j].qt;
                }
                itype = atomData[i].type;
                _eta = atom[itype].eta;
                atomData[i].gs = atom[itype].chi - _eta * atomData[i].qs - gssum;
                atomData[i].gt = 1.0 - _eta * atomData[i].qt - gtsum;
            }

            const Containers::StaticArray<2, Double> Gold(Gnew);

            Gnew[0] = dot_product(1);
            Gnew[1] = dot_product(2);

            for (std::size_t i = 0; i < NATOMS; ++i)
            {
                atomData[i].hs = atomData[i].gs + (Gnew[0] / Gold[0]) * atomData[i].hs;
                atomData[i].ht = atomData[i].gt + (Gnew[1] / Gold[1]) * atomData[i].ht;
            }
        }
    }

    void BOPRIM()
    {
        std::size_t j, inxn;
        Double dr2, _bo_sum, dBOp;
        Vector3d dr, _bo, _dln_BOp;

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            atomData[i].deltap[0] = -atom[atomData[i].type].Val;
        }

        for (std::size_t i = 0; i < NATOMS; ++i)
        {

            for (std::size_t b = 0; b < atomData[i].bonds.size(); ++b)
            {
                j = atomData[i].bonds[b];

                if (j < i)
                    continue;

                inxn = atom[atomData[i].type].inxn2[atomData[j].type] - 1;
                dr = atomData[i].position - atomData[j].position;
                // TODO idk if 0.5 * dr is correct.
                dr2 = (0.5 * dr).dot();

                arg_BOpij[0] = bond[inxn].cBOp1 * pow(Double(dr2), bond[inxn].pbo2h);
                arg_BOpij[1] = bond[inxn].cBOp3 * pow(Double(dr2), bond[inxn].pbo4h);
                arg_BOpij[2] = bond[inxn].cBOp5 * pow(Double(dr2), bond[inxn].pbo6h);

                for (std::size_t k = 0; k < 3; ++k)
                {
                    _bo[k] = bond[inxn].swh[k] * exp(arg_BOpij[k]);
                }

                // Small modification exists in sigma-bond prime, see reac.f line 4444. sigma-bond prime is multiplied by
                // (1.d0 + 1.d-4) here.  Later in original reaxff code, sigma-bond prime is subtracted by
                // 0.01*vpar30 resulting in the following subtractions in abo(i1), bo(nbon), bos(nbon), bosi(nbon).
                // However, this modification is only applied to the energy calculation, not to the force calculation.
                // The subtraction by <cutoff_vpar30> is done after the derivative calculations so that the variables in the
                // force-calc routines use the original BOp values and the ones in the energy-calc routines are the modified value.

                _bo[0] *= (1.0 + cutoff_vpar30);

                // If the total <BOp> before the subtraction is greater than <cutoff_vpar30>,
                // get "final" <BOp> value and its derivatives.

                if (_bo.sum() > cutoff_vpar30)
                {
                    _dln_BOp[0] = bond[inxn].swh[0] * bond[inxn].pbo2 * arg_BOpij[0];
                    _dln_BOp[1] = bond[inxn].swh[1] * bond[inxn].pbo4 * arg_BOpij[1];
                    _dln_BOp[2] = bond[inxn].swh[2] * bond[inxn].pbo6 * arg_BOpij[2];
                    _dln_BOp = _dln_BOp / dr2;
                    arrayAppend(atomData[i].dln_BOp, InPlaceInit, _dln_BOp);
                    arrayAppend(atomData[j].dln_BOp, InPlaceInit, _dln_BOp);

                    dBOp = Vector3d(_bo * _dln_BOp).sum();
                    arrayAppend(atomData[i].dBOp, InPlaceInit, dBOp);
                    arrayAppend(atomData[j].dBOp, InPlaceInit, dBOp);
                    // After the derivative calculations are done, do the subtraction described above
                    // which results in the difference of bond-order  between the energy-calc and the force-calc.
                    _bo[0] -= cutoff_vpar30;
                    _bo_sum = _bo.sum();
                    arrayAppend(atomData[i].bo_sum, InPlaceInit, _bo_sum);
                    arrayAppend(atomData[j].bo_sum, InPlaceInit, _bo_sum);
                    arrayAppend(atomData[i].bo, InPlaceInit, _bo);
                    arrayAppend(atomData[j].bo, InPlaceInit, _bo);
                    atomData[i].deltap[0] += _bo_sum;
                    atomData[j].deltap[0] += _bo_sum;
                }
                else
                {
                    arrayAppend(atomData[i].dln_BOp, InPlaceInit, Vector3d{0.0});
                    arrayAppend(atomData[j].dln_BOp, InPlaceInit, Vector3d{0.0});
                    arrayAppend(atomData[i].dBOp, InPlaceInit, 0.0);
                    arrayAppend(atomData[j].dBOp, InPlaceInit, 0.0);
                    arrayAppend(atomData[i].bo, InPlaceInit, Vector3d{0.0});
                    arrayAppend(atomData[i].bo, InPlaceInit, Vector3d{0.0});
                    arrayAppend(atomData[i].bo_sum, InPlaceInit, 0.0);
                    arrayAppend(atomData[j].bo_sum, InPlaceInit, 0.0);
                }
            }
        }
    }

    void BOFULL()
    {
        std::size_t j, inxn;
        Double exppboc1i, exppboc2i, exppboc1j, exppboc2j;
        Double fn1, fn2, fn3, fn23, fn4, fn5, fn45, fn145, fn1145;
        Double BOp0, BOpsqr, _BO_sum, BOpij_2, _BO_nsum;
        Vector3d _BO;

        Double u1ij, u1ji, u1ij_inv2, u1ji_inv2, u45ij, u45ji;
        Double Cf1ij, Cf1ji, Cf1Aij, Cf1Bij, Cf45ij, Cf45ji, Cf1ij_div1, Cf1ji_div1;
        Double exph_45ij, exph_45ji, exp1, exp2, exp12, exp_delt22;
        Double pboc34, fn45_inv;

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            atomData[i].deltap[1] = atomData[i].deltap[0] + atom[atomData[i].type].Val - atom[atomData[i].type].Valval;
        }

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            exppboc1i = exp(-vpar1 * atomData[i].deltap[0]);
            exppboc2i = exp(-vpar2 * atomData[i].deltap[0]);

            fn2 = exppboc1i + exppboc1j;
            fn3 = (-1.0 / vpar2) * log(0.5 * (exppboc2i + exppboc2j));
            fn23 = fn2 + fn3;

            for (std::size_t b = 0; b < atomData[i].bonds.size(); ++b)
            {

                j = atomData[i].bonds[b];

                if (j < i)
                    continue;

                exppboc1j = exp(-vpar1 * atomData[j].deltap[0]);
                exppboc2j = exp(-vpar2 * atomData[j].deltap[0]);

                inxn = atom[atomData[i].type].inxn2[atomData[j].type] - 1;
                BOp0 = atomData[i].bo_sum[b];

                fn1 = 0.5 * ((atom[atomData[i].type].Val + fn2) / (atom[atomData[i].type].Val + fn23) + (atom[atomData[j].type].Val + fn2) / (atom[atomData[j].type].Val + fn23));
                // TODO ovc is either 1 or 0, so this doesn't make any sense. Probably some scaling i have missed or Fortran might not have booleans
                if (bond[inxn].ovc < 1.0e-3)
                    fn1 = 1.0;

                BOpsqr = BOp0 * BOp0;
                fn4 = 1.0 / (1.0 + exp(-bond[inxn].pboc3 * (bond[inxn].pboc4 * BOpsqr - atomData[i].deltap[1])) + bond[inxn].pboc5);
                fn5 = 1.0 / (1.0 + exp(-bond[inxn].pboc3 * (bond[inxn].pboc4 * BOpsqr - atomData[j].deltap[1])) + bond[inxn].pboc5);
                if (bond[inxn].v13cor < 1.0e-3)
                {
                    fn4 = 1.0;
                    fn5 = 1.0;
                }

                fn45 = fn4 * fn5;
                fn145 = fn1 * fn45;
                fn1145 = fn1 * fn145;

                // New bond order definition
                _BO_sum = BOp0 * fn145;
                _BO[1] = atomData[i].bo[b][1] * fn1145;
                _BO[2] = atomData[i].bo[b][2] * fn1145;
                if (_BO_sum < 1.0e-10)
                    _BO_sum = 0.0;
                if (_BO[1] < 1.0e-10)
                    _BO[1] = 0.0;
                if (_BO[2] < 1.0e-10)
                    _BO[2] = 0.0;

                // new sigma BO definition
                _BO[0] = _BO_sum - _BO[1] - _BO[2];
                arrayAppend(atomData[i].BO, InPlaceInit, _BO);
                arrayAppend(atomData[j].BO, InPlaceInit, _BO);

                arrayAppend(atomData[i].BO_sum, InPlaceInit, _BO_sum);
                arrayAppend(atomData[j].BO_sum, InPlaceInit, _BO_sum);

                // CALCULATION OF DERIVATIVE OF BOND ORDER
                // all following comes from Coding Methodology section:
                // part 1:

                u1ij = atom[atomData[i].type].Val + fn23;
                u1ji = atom[atomData[j].type].Val + fn23;

                // part 2:
                u1ij_inv2 = 1.0 / (u1ij * u1ij);
                u1ji_inv2 = 1.0 / (u1ji * u1ji);

                Cf1Aij = 0.5 * fn3 * (u1ij_inv2 + u1ji_inv2);
                Cf1Bij = -0.5 * ((u1ij - fn3) * u1ij_inv2 + (u1ji - fn3) * u1ji_inv2);

                // part 3:
                exp_delt22 = exppboc2i + exppboc2j;
                Cf1ij = (-Cf1Aij * bond[inxn].pboc1 * exppboc1i) + (Cf1Bij * exppboc2i) / (exp_delt22);
                Cf1ji = (-Cf1Aij * bond[inxn].pboc1 * exppboc1j) + (Cf1Bij * exppboc2j) / (exp_delt22);

                // part 4:
                pboc34 = bond[inxn].pboc3 * bond[inxn].pboc4;
                BOpij_2 = BOpsqr;

                u45ij = bond[inxn].pboc5 + bond[inxn].pboc3 * atomData[i].deltap[1] - pboc34 * BOpij_2;
                u45ji = bond[inxn].pboc5 + bond[inxn].pboc3 * atomData[j].deltap[1] - pboc34 * BOpij_2;

                // part 5:
                exph_45ij = exp(u45ij);
                exph_45ji = exp(u45ji);
                exp1 = 1.0 / (1.0 + exph_45ij);
                exp2 = 1.0 / (1.0 + exph_45ji);
                exp12 = exp1 * exp2;

                Cf45ij = -exph_45ij * exp12 * exp1;
                Cf45ji = -exph_45ji * exp12 * exp2;

                // if following conditions, <ovc> and <v13cor>, are satisfied, correction terms
                // fn1 and/or fn4 & fn5 are 1.d0 and their derivatives are zero

                if (bond[inxn].ovc < 1.0e-3)
                {
                    Cf1ij = 0.0;
                    Cf1ji = 0.0;
                }
                if (bond[inxn].v13cor < 1.0e-3)
                {
                    Cf45ij = 0.0;
                    Cf45ji = 0.0;
                }

                // part 6:
                fn45_inv = 1.0 / fn45;
                Cf1ij_div1 = Cf1ij / fn1;
                Cf1ji_div1 = Cf1ji / fn1;

                arrayAppend(atomData[i].A0, InPlaceInit, fn145);
                arrayAppend(atomData[i].A1, InPlaceInit, -2.0 * pboc34 * BOp0 * (Cf45ij + Cf45ji) * fn45_inv);
                arrayAppend(atomData[i].A2, InPlaceInit, Cf1ij_div1 + (bond[inxn].pboc3 * Cf45ij * fn45_inv));
                arrayAppend(atomData[i].A3, InPlaceInit, atomData[i].A2[atomData[i].A2.size() - 1] + Cf1ij_div1);

                arrayAppend(atomData[j].A0, InPlaceInit, atomData[i].A0[atomData[i].A0.size() - 1]);
                arrayAppend(atomData[j].A1, InPlaceInit, atomData[i].A1[atomData[i].A1.size() - 1]);
                arrayAppend(atomData[j].A2, InPlaceInit, Cf1ji_div1 + (bond[inxn].pboc3 * Cf45ji * fn45_inv));
                arrayAppend(atomData[j].A3, InPlaceInit, atomData[i].A2[atomData[i].A2.size() - 1] + Cf1ji_div1);
            }
        }

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            _BO_nsum = 0.0;
            for (std::size_t b = 0; b < atomData[i].BO_sum.size(); ++b)
            {
                _BO_nsum += atomData[i].BO_sum[b];
            }
            atomData[i].delta = -atom[atomData[i].type].Val - atom[atomData[i].type].Valval + _BO_nsum;
        }
    }

    void ENbond()
    {
        std::size_t j, itb, itb1, inxn, itype;
        Double drtb, qij;
        Double PEvdw, PEclmb;
        Double CEvdw, CEclmb;
        Vector3d dr, ff;

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            itype = atomData[i].type;
            PE[13] = PE[13] + CEchrge * (atom[itype].chi * atomData[i].q + pow(0.5 * atom[itype].eta * atomData[i].q, 2));

            for (std::size_t n = 0; n < atomData[i].neighbors.size(); ++n)
            {
                j = atomData[i].neighbors[n];

                if (j < i)
                    continue;

                itb = Int(atomData[i].dpq2[n] * Double(UDRi));
                itb1 = itb + 1;
                drtb = atomData[i].dpq2[n] - itb * Double(UDR);
                drtb = drtb * Double(UDRi);

                inxn = atom[atomData[i].type].inxn2[atomData[j].type] - 1;

                // van der Waals: TBL_Evdw[0] = TBL_Evdw_p ; TBL_Evdw[1] = TBL_Evdw_d
                PEvdw = (1.0 - drtb) * TBL_Evdw_p[itb][inxn] + drtb * TBL_Evdw_p[itb1][inxn];
                CEvdw = (1.0 - drtb) * TBL_Evdw_d[itb][inxn] + drtb * TBL_Evdw_d[itb1][inxn];

                // Coulomb:
                qij = atomData[i].q * atomData[j].q;
                PEclmb = (1.0 - drtb) * TBL_Eclmb_p[itb][inxn] + drtb * TBL_Eclmb_p[itb1][inxn];
                PEclmb *= qij;
                CEclmb = (1.0 - drtb) * TBL_Eclmb_d[itb][inxn] + drtb * TBL_Eclmb_d[itb1][inxn];
                CEclmb *= qij;

                PE[11] = PE[11] + PEvdw;
                PE[12] = PE[12] + PEclmb;

                dr = atomData[i].position - atomData[j].position;
                ff = (CEvdw + CEclmb) * dr;

                atomData[i].force = atomData[i].force - ff;
                atomData[j].force = atomData[j].force + ff;
            }
        }
    }

    void ForceBbo(const std::size_t &i, const std::size_t &j, const std::size_t &b, const Vector3d &coeff)
    {
        // * With the new bond-order definition, 1st term is the derivative of "full"-bond order,
        // * 2nd is for pi-bond order and 3rd is for pipi-bond order.
        Vector3d Cbond, dr, ff, cBO;

        Vector3d cf = Vector3d(coeff[0], coeff[1] - coeff[0], coeff[2] - coeff[0]);
        Cbond[0] = cf[0] * (atomData[i].A0[b] + atomData[i].BO_sum[b] * atomData[i].A1[b]) * atomData[i].dBOp[b] +
                   cf[1] * atomData[i].BO[b][1] * (atomData[i].dln_BOp[b][1] + atomData[i].A1[b] * atomData[i].dBOp[b]) +
                   cf[2] * atomData[i].BO[b][2] * (atomData[i].dln_BOp[b][2] + atomData[i].A1[b] * atomData[i].dBOp[b]);

        dr = atomData[i].position - atomData[j].position;
        ff = Cbond[0] * dr;

        atomData[i].force = atomData[i].force - ff;
        atomData[j].force = atomData[j].force + ff;

        // * 1st element is "full"-bond order.
        cBO = Vector3d(cf[0] * atomData[i].BO_sum[b], cf[1] * atomData[i].BO[b][1], cf[2] * atomData[i].BO[b][2]);

        Cbond[1] = cBO[0] * atomData[i].A2[b] + (cBO[1] + cBO[2]) * atomData[i].A3[b];

        for (std::size_t bj = 0; bj < atomData[j].bonds.size(); ++bj)
        {
            if (atomData[j].bonds[bj] != i)
                continue;

            Cbond[2] = cBO[0] * atomData[j].A2[bj] + (cBO[1] + cBO[2]) * atomData[j].A3[bj];
        }

        atomData[i].ccbnd = atomData[i].ccbnd + Cbond[1];
        atomData[j].ccbnd = atomData[j].ccbnd + Cbond[2];
    }

    void Ebond()
    {
        std::size_t j, inxn;
        Double exp_be12, PEbo, CEbo;
        Vector3d coeff;

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            for (std::size_t b = 0; b < atomData[i].bonds.size(); ++b)
            {
                j = atomData[i].bonds[b];
                // TODO might be wrong... maybe (probably only calculate bond order for nearest neighbor)
                if (j < i)
                    continue;

                inxn = atom[atomData[i].type].inxn2[atomData[j].type] - 1;

                exp_be12 = exp(bond[inxn].pbe1 * (1.0 - pow(atomData[i].BO[b][0], bond[inxn].pbe2)));
                PEbo = -bond[inxn].Desig * atomData[i].BO[b][0] * exp_be12 - bond[inxn].Depi * atomData[i].BO[b][1] - bond[inxn].Depipi * atomData[i].BO[b][2];

                PE[1] = PE[1] + PEbo;

                CEbo = -bond[inxn].Desig * exp_be12 * (1.0 - bond[inxn].pbe1 * bond[inxn].pbe2 * pow(atomData[i].BO[b][0], bond[inxn].pbe2));
                coeff = Vector3d(CEbo, -bond[inxn].Depi, -bond[inxn].Depipi);
                ForceBbo(i, j, b, coeff);
            }
        }
    }

    void Elnpr()
    {
        // TODO proper way of writing functions here, maybe initialize some array values as variables
        std::size_t j, inxn, itype, idEh;
        Vector3d coeff;

        // *Lone Pair Energy Terms
        Double Clp, CElp, PElp, dEh;
        Double explp1, expvd2, dElp, deltaE;

        // *Overcoordination Energy Terms
        Double sum_ovun1, sum_ovun2;
        Double deltalpcorr, PEover, DlpV_i;
        Double expovun2, expovun1;

        // *Undercoordination Energy Terms
        Double expovun2n, expovun6, expovun8;
        Double PEunder;

        Double CElp_b, CElp_d, CElp_bpp;

        Double div_expovun2, div_expovun2n, div_expovun1, div_expovun8;

        // preperation

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            itype = atomData[i].type;

            if (itype == 0)
                continue;

            deltaE = -atom[itype].Vale + atom[itype].Val + atomData[i].delta;

            dEh = deltaE * 0.5;
            // idEh = is_idEh * int(dEh)
            idEh = Int(dEh);
            explp1 = exp(-atom[itype].plp1 * pow(2.0 + deltaE - 2 * idEh, 2));
            Clp = 2.0 * atom[itype].plp1 * explp1 * (2.0 + deltaE - 2 * idEh);

            atomData[i].dDlp = Clp;

            atomData[i].nlp = explp1 - Double(idEh);
            atomData[i].deltalp = atom[itype].nlpopt - atomData[i].nlp;

            if (atom[itype].mass > 21.0)
                atomData[i].deltalp = 0.0;
        }

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            itype = atomData[i].type;

            sum_ovun1 = 0.0;
            sum_ovun2 = 0.0;

            for (std::size_t b = 0; b < atomData[i].bonds.size(); ++b)
            {
                j = atomData[i].bonds[b];
                inxn = atom[itype].inxn2[atomData[j].type] - 1;

                sum_ovun1 += bond[inxn].povun1 * bond[inxn].Desig * atomData[i].BO_sum[b];
                sum_ovun2 += (atomData[j].delta - atomData[j].deltalp) * (atomData[i].BO[b][1] + atomData[i].BO[b][2]);
            }

            // *Lone Pair
            expvd2 = exp(-75.0 * atomData[i].deltalp);
            dElp = atom[itype].plp2 * ((1.0 + expvd2) + 75.0 * atomData[i].deltalp * expvd2) / pow(1.0 + expvd2, 2);

            // *Over Coordinate + Common part with Under Coordinate
            expovun1 = atom[itype].povun3 * exp(atom[itype].povun4 * sum_ovun2);
            deltalpcorr = atomData[i].delta - atomData[i].deltalp / (1.0 + expovun1);
            expovun2 = exp(atom[itype].povun2 * deltalpcorr);

            // *if one atom flys away, this term becomes zero because the total bond-order becomes zero.
            // *Add a small value in the denominator to avoid it. See poten.f line 787,//
            // *hulpp=(1.0/(vov1+aval(ity1)+1e-8))
            DlpV_i = 1.0 / (deltalpcorr + atom[itype].Val + 1.0e-8);

            // *Under Coordinate
            expovun2n = 1.0 / expovun2;
            expovun6 = exp(atom[itype].povun6 * deltalpcorr);
            expovun8 = atom[itype].povun7 * exp(atom[itype].povun8 * sum_ovun2);

            div_expovun1 = 1.0 / (1.0 + expovun1);
            div_expovun2 = 1.0 / (1.0 + expovun2);
            div_expovun2n = 1.0 / (1.0 + expovun2n);
            div_expovun8 = 1.0 / (1.0 + expovun8);

            // *Energy Calculation
            PElp = atom[itype].plp2 * atomData[i].deltalp / (1.0 + expvd2);
            PEover = sum_ovun1 * DlpV_i * deltalpcorr * div_expovun2;
            PEunder = atom[itype].povun5 * (1.0 - expovun6) * div_expovun2n * div_expovun8;

            // *if the representitive atom is a resident ,sum their potential energies.
            PE[2] = PE[2] + PElp;
            PE[3] = PE[3] + PEover;
            PE[4] = PE[4] + PEunder;

            // *Coefficient Calculation
            CElp = dElp * atomData[i].dDlp;

            CEover[0] = deltalpcorr * DlpV_i * div_expovun2;
            CEover[1] = sum_ovun1 * DlpV_i * div_expovun2 * (1.0 - deltalpcorr * DlpV_i - atom[itype].povun2 * deltalpcorr * div_expovun2n);

            CEover[2] = CEover[1] * (1.0 - atomData[i].dDlp * div_expovun1);
            CEover[3] = CEover[1] * atomData[i].deltalp * atom[itype].povun4 * expovun1 * pow(div_expovun1, 2);

            CEunder[0] = (atom[itype].povun5 * atom[itype].povun6 * expovun6 * div_expovun8 + PEunder * atom[itype].povun2 * expovun2n) * div_expovun2n;
            CEunder[1] = -PEunder * atom[itype].povun8 * expovun8 * div_expovun8;
            CEunder[2] = CEunder[0] * (1.0 - atomData[i].dDlp * div_expovun1);
            CEunder[3] = CEunder[0] * atomData[i].deltalp * atom[itype].povun4 * expovun1 * pow(div_expovun1, 2) + CEunder[1];

            // *Force Calculation
            for (std::size_t b = 0; b < atomData[i].bonds.size(); ++b)
            {
                j = atomData[i].bonds[b];
                inxn = atom[itype].inxn2[atomData[j].type] - 1;

                CEover[4] = CEover[0] * bond[inxn].povun1 * bond[inxn].Desig;
                CEover[5] = CEover[3] * (1.0 - atomData[j].dDlp) * (atomData[i].BO[b][1] + atomData[i].BO[b][2]);
                CEover[6] = CEover[3] * (atomData[j].delta - atomData[j].deltalp);

                CEunder[4] = CEunder[3] * (1.0 - atomData[j].dDlp * (atomData[i].BO[b][1] + atomData[i].BO[b][2]));
                CEunder[5] = CEunder[3] * (atomData[j].delta - atomData[j].deltalp);

                CElp_b = CElp + CEover[2] + CEover[4] + CEunder[2];
                CElp_bpp = CEover[6] + CEunder[5];

                // coeff(1:3) = CElp_b + (/ 0.0, CElp_bpp, CElp_bpp/)
                coeff = Vector3d(CElp_b, CElp_bpp + CElp_b, CElp_bpp + CElp_b);
                ForceBbo(i, j, b, coeff);

                CElp_d = CEover[5] + CEunder[5];
                atomData[j].cdbnd = atomData[j].cdbnd + CElp_d;
            }
        }
    }

    void ForceB(const std::size_t &i, const std::size_t &j, const std::size_t &b, const Double &coeff)
    {
        Vector3d Cbond, dr, ff;

        Cbond[0] = coeff * (atomData[i].A0[b] + atomData[i].BO_sum[b] * atomData[i].A1[b]);
        dr = atomData[i].position - atomData[j].position;
        ff = Cbond[0] * atomData[i].dBOp[b] * dr;

        atomData[i].force = atomData[i].force - ff;
        atomData[i].force = atomData[j].force + ff;

        // A3 is not necessary anymore with the new BO def.
        Cbond[1] = coeff * atomData[i].BO_sum[b] * atomData[i].A2[b];
        Cbond[2] = coeff * atomData[i].BO_sum[b] * atomData[j].A2[b];

        atomData[i].ccbnd = atomData[i].ccbnd + Cbond[1];
        atomData[j].ccbnd = atomData[j].ccbnd + Cbond[2];
    }

    void ForceA3(const Double &coeff, const std::size_t &i, const std::size_t &j, const std::size_t &k, const Vector3d &da0, const Vector3d &da1, const Double &da0_0, const Double &da1_0)
    {
        Vector3d Ci, Ck;
        Vector3d fij, fjk, fijjk;
        Containers::StaticArray<2, Vector2d> Caa;
        Double CCisqr, coCC;

        Caa[1][1] = pow(da0_0, 2);
        Caa[1][0] = (da0 * da1).sum();
        Caa[0][1] = Caa[1][0];
        Caa[0][0] = pow(da1_0, 2);

        CCisqr = 1.0 / (da0_0 * da1_0);
        coCC = coeff * CCisqr;

        // Some of calculations are unnecessary due to the action-reaction relation.
        Ci[0] = -(Caa[1][0] / Caa[1][1]);
        Ci[1] = 1.0;

        Ck[0] = -1.0;
        Ck[1] = Caa[1][0] / Caa[0][0];

        fij = coCC * (Ci[0] * da0 + Ci[1] * da1);
        fjk = coCC * (Ck[0] * da0 + Ck[1] * da1);
        fijjk = -fij + fjk;

        atomData[i].force = atomData[i].force + fij;
        atomData[j].force = atomData[j].force + fijjk;
        atomData[k].force = atomData[k].force - fjk;
    }

    void Ehb()
    {
        std::size_t j, k, itype, jtype, ktype, inxnhb;
        Vector3d rjk, rij;
        Double rjk0, rij0;
        Double cos_ijk, theta_ijk, sin_ijk_half, sin_xhz4, cos_xhz1;
        Double exp_hb2, exp_hb3;
        Double PEhb;
        Vector3d CEhb, ff;

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            itype = atomData[i].type;
            for (std::size_t b = 0; b < atomData[i].bonds.size(); ++b)
            {
                j = atomData[i].bonds[b];
                jtype = atomData[j].type;

                if (jtype == 2 && atomData[i].BO_sum[b] > MINBO0)
                {
                    for (std::size_t ki; ki < atomData[i].neighbors.size(); ++ki)
                    {
                        k = atomData[i].neighbors[ki];
                        ktype = atomData[k].type;

                        inxnhb = atom[itype].inxn3hb[jtype][ktype];

                        if (j != k && i != k && inxnhb != 0)
                        {
                            inxnhb -= 1;

                            if (atomData[i].dpq2[k] >= Double(rctap2))
                                continue;

                            rjk = atomData[j].position - atomData[k].position;
                            rjk0 = sqrt((rjk * rjk).sum());

                            rij = atomData[i].position - atomData[j].position;
                            rij0 = sqrt((rij * rij).sum());

                            cos_ijk = -(rij * rjk).sum() / (rij0 * rjk0);
                            if (cos_ijk > MAXANGLE)
                                cos_ijk = MAXANGLE;
                            if (cos_ijk < MINANGLE)
                                cos_ijk = MINANGLE;

                            theta_ijk = acos(cos_ijk);

                            sin_ijk_half = sin(0.5 * theta_ijk);
                            sin_xhz4 = pow(sin_ijk_half, 4);
                            cos_xhz1 = (1.0 - cos_ijk);

                            exp_hb2 = exp(-h_bond[inxnhb].phb2 * atomData[i].BO_sum[b]);
                            exp_hb3 = exp(-h_bond[inxnhb].phb3 * (h_bond[inxnhb].r0hb / rjk0 + rjk0 / h_bond[inxnhb].r0hb - 2.0));

                            PEhb = h_bond[inxnhb].phb1 * (1.0 - exp_hb2) * exp_hb3 * sin_xhz4;

                            PE[10] = PE[10] + PEhb;

                            CEhb[0] = h_bond[inxnhb].phb1 * h_bond[inxnhb].phb2 * exp_hb2 * exp_hb3 * sin_xhz4;
                            CEhb[1] = -0.5 * h_bond[inxnhb].phb1 * (1.0 - exp_hb2) * exp_hb3 * cos_xhz1;
                            CEhb[2] = -PEhb * h_bond[inxnhb].phb3 * (-h_bond[inxnhb].r0hb / pow(rjk0, 2) + 1.0 / h_bond[inxnhb].r0hb) * (1.0 / rjk0);

                            ForceB(i, j, b, CEhb[0]);
                            ForceA3(CEhb[1], i, j, k, rij, rjk, rij0, rjk0);

                            ff = CEhb[2] * rjk;

                            atomData[j].force = atomData[k].force - ff;
                            atomData[k].force = atomData[k].force + ff;
                        }
                    }
                }
            }
        }
    }

    void E3b()
    {
        std::size_t i, k, itype, jtype, ktype, inxn;
        // Valency Energy Calculation Terms:
        Double PEval, fn7ij, fn7jk, fn8j, delta_ang;
        Vector3d rij, rjk;
        Double rij0, rjk0, SBO2, SBO;
        Double sum_BO8, theta_ijk, theta0, theta_diff, exp2;
        Double sin_ijk, cos_ijk;
        Double Cf7ij, Cf7jk, Cf8j, Ctheta0, CSBO2; // ! [unused] Ctheta_diff
        Vector2d dsBO;
        Double BOij_p4, BOjk_p4, exp3ij, exp3jk, exp6, exp7, trm8;
        Double sum_SBO1, prod_SBO;

        // Penalty Energy Calculation Terms:
        Double fn9, PEpen;
        Double exp_pen3, exp_pen4, exp_pen2ij, exp_pen2jk;
        Double trm_pen34, Cf9j;

        // Conjugation Energy (3body) Calculation Terms:
        Double delta_val, PEcoa;
        Double sum_BOi, sum_BOk, exp_coa3i, exp_coa3k, exp_coa4i, exp_coa4k;
        Double exp_coa2;

        Containers::StaticArray<9, Double> CEval;
        Containers::StaticArray<5, Double> CEcoa;
        Vector3d CEpen, CE3body_d, CE3body_b, coeff;
        Double BOij, BOjk, CE3body_a, coeff_val;

        for (std::size_t j = 0; j < NATOMS; ++j)
        {
            jtype = atomData[j].type;
            sum_BO8 = 0.0;
            sum_SBO1 = 0.0;

            for (std::size_t b = 0; b < atomData[j].bonds.size(); ++b)
            {
                sum_BO8 = sum_BO8 - pow(atomData[j].BO_sum[b], 8);
                sum_SBO1 = sum_SBO1 + atomData[j].BO[b][1] + atomData[j].BO[b][2];
            }

            prod_SBO = exp(sum_BO8);
            delta_ang = atomData[j].delta + atom[jtype].Val - atom[jtype].Valangle;

            for (std::size_t b0 = 0; b0 + 1 < atomData[j].bonds.size(); ++b0)
            {
                // cutof2_esub is used as the BO cutoff in the original ReaxFF code.
                BOij = atomData[j].BO_sum[b0] - cutof2_esub;
                if (BOij > 0.0)
                {
                    i = atomData[j].bonds[b0];
                    itype = atomData[i].type;

                    rij = atomData[i].position - atomData[j].position;
                    rij0 = sqrt((rij * rij).sum());

                    for (std::size_t b1 = 1; b1 < atomData[j].bonds.size(); ++b1)
                    {
                        BOjk = atomData[j].BO_sum[b1] - cutof2_esub;

                        if (BOjk <= 0.0)
                            continue;

                        if (atomData[j].BO_sum[b0] * atomData[j].BO_sum[b1] <= cutof2_esub)
                            continue;

                        k = atomData[j].bonds[b1];
                        ktype = atomData[k].type;

                        rjk = atomData[j].position - atomData[k].position;
                        rjk0 = sqrt((rjk * rjk).sum());

                        cos_ijk = -(rij * rjk).sum() / (rij0 * rjk0);
                        if (cos_ijk > MAXANGLE)
                            cos_ijk = MAXANGLE;
                        if (cos_ijk < MINANGLE)
                            cos_ijk = MINANGLE;

                        theta_ijk = acos(cos_ijk);
                        sin_ijk = sin(theta_ijk);

                        // Check the type of 3atoms combination
                        inxn = atom[itype].inxn3[jtype][ktype];
                        if (inxn == 0)
                            continue;

                        inxn -= 1;

                        // PEval part:
                        BOij_p4 = pow(BOij, angle[inxn].pval4);
                        exp3ij = exp(-atom[jtype].pval3 * BOij_p4);
                        fn7ij = 1.0 - exp3ij;
                        BOjk_p4 = pow(BOjk, angle[inxn].pval4);
                        exp3jk = exp(-atom[jtype].pval3 * BOjk_p4);
                        fn7jk = 1.0 - exp3jk;

                        exp6 = exp(angle[inxn].pval6 * delta_ang);
                        exp7 = exp(-angle[inxn].pval7 * delta_ang);
                        trm8 = 1.0 + exp6 + exp7;

                        fn8j = atom[jtype].pval5 - (atom[jtype].pval5 - 1.0) * (2.0 + exp6) / trm8;

                        SBO = sum_SBO1 + (1.0 - prod_SBO) * (-delta_ang - angle[inxn].pval8 * atomData[j].nlp);
                        if (SBO <= 0)
                        {
                            SBO2 = 0.0;
                        }
                        else if (SBO <= 1)
                        {
                            SBO2 = pow(SBO, angle[inxn].pval9);
                        }
                        else if (SBO <= 2)
                        {
                            SBO2 = 2.0 - pow(2.0 - SBO, angle[inxn].pval9);
                        }
                        else
                        {
                            SBO2 = 2.0;
                        }

                        theta0 = Constantsd::pi() - angle[inxn].theta00 * (1.0 - exp(-angle[inxn].pval10 * (2.0 - SBO2)));
                        theta_diff = theta0 - theta_ijk;
                        exp2 = exp(-angle[inxn].pval2 * theta_diff * theta_diff);

                        PEval = fn7ij * fn7jk * fn8j * (angle[inxn].pval1 - angle[inxn].pval1 * exp2);

                        // PEval derivative part:
                        Cf7ij = atom[jtype].pval3 * angle[inxn].pval4 * pow(BOij, angle[inxn].pval4 - 1.0) * exp3ij;
                        Cf7jk = atom[jtype].pval3 * angle[inxn].pval4 * pow(BOjk, angle[inxn].pval4 - 1.0) * exp3jk;

                        Cf8j = (1.0 - atom[jtype].pval5) / (trm8 * trm8) * (angle[inxn].pval6 * exp6 * trm8 - (2.0 + exp6) * (angle[inxn].pval6 * exp6 - angle[inxn].pval7 * exp7));
                        // ! [unused] Ctheta_diff = 2.0 * angle[inxn].pval2 * theta_diff * exp2 / (1.0 - exp2);
                        Ctheta0 = angle[inxn].pval10 * angle[inxn].theta00 * exp(-angle[inxn].pval10 * (2.0 - SBO2));

                        if (SBO <= 0 || SBO > 2)
                        {
                            CSBO2 = 0.0;
                        }
                        else if (SBO > 0 && SBO <= 1)
                        {
                            CSBO2 = angle[inxn].pval9 * pow(SBO, angle[inxn].pval9 - 1.0);
                        }
                        else if (SBO > 1 && SBO <= 2)
                        {
                            CSBO2 = angle[inxn].pval9 * pow(2.0 - SBO, angle[inxn].pval9 - 1.0);
                        }

                        dsBO[0] = -8.0 * prod_SBO * (delta_ang + angle[inxn].pval8 * atomData[j].nlp);
                        dsBO[1] = (prod_SBO - 1.0) * (1.0 - angle[inxn].pval8 * atomData[j].dDlp);

                        CEval[0] = Cf7ij * fn7jk * fn8j * angle[inxn].pval1 * (1.0 - exp2);
                        CEval[1] = fn7ij * Cf7jk * fn8j * angle[inxn].pval1 * (1.0 - exp2);
                        CEval[2] = fn7ij * fn7jk * Cf8j * angle[inxn].pval1 * (1.0 - exp2);

                        CEval[3] = 2.0 * angle[inxn].pval1 * angle[inxn].pval2 * fn7ij * fn7jk * fn8j * exp2 * theta_diff;
                        CEval[4] = CEval[3] * Ctheta0 * CSBO2;
                        CEval[5] = CEval[4] * dsBO[0];
                        CEval[6] = CEval[4] * dsBO[1];
                        CEval[7] = CEval[3] / sin_ijk;

                        // PEpen part:
                        exp_pen3 = exp(-angle[inxn].ppen3 * atomData[j].delta);
                        exp_pen4 = exp(angle[inxn].ppen4 * atomData[j].delta);

                        fn9 = (2.0 + exp_pen3) / (1.0 + exp_pen3 + exp_pen4);
                        exp_pen2ij = exp(-angle[inxn].ppen2 * (BOij - 2.0) * (BOij - 2.0));
                        exp_pen2jk = exp(-angle[inxn].ppen2 * (BOjk - 2.0) * (BOjk - 2.0));

                        PEpen = angle[inxn].ppen1 * fn9 * exp_pen2ij * exp_pen2jk;

                        // PEpen derivative part:
                        trm_pen34 = 1.0 * exp_pen3 + exp_pen4;
                        Cf9j = (-angle[inxn].ppen3 * exp_pen3 * trm_pen34 - (2.0 + exp_pen3) * (-angle[inxn].ppen3 * exp_pen3 + angle[inxn].ppen4 * exp_pen4)) / (trm_pen34 * trm_pen34);

                        CEpen[0] = Cf9j / fn9;
                        CEpen[1] = -2.0 * angle[inxn].ppen2 * (BOij - 2.0);
                        CEpen[2] = -2.0 * angle[inxn].ppen2 * (BOjk - 2.0);

                        CEpen = CEpen * PEpen;

                        // PEcoa part:
                        sum_BOi = atomData[i].delta + atom[itype].Val;
                        sum_BOk = atomData[k].delta + atom[ktype].Val;
                        delta_val = atomData[j].delta + atom[jtype].Val - atom[jtype].Valval;

                        exp_coa2 = exp(angle[inxn].pcoa2 * delta_val);
                        exp_coa3i = exp(-angle[inxn].pcoa3 * pow(-BOij + sum_BOi, 2));
                        exp_coa3k = exp(-angle[inxn].pcoa3 * pow(-BOjk + sum_BOk, 2));
                        exp_coa4i = exp(-angle[inxn].pcoa4 * pow(BOij - 1.5, 2));
                        exp_coa4k = exp(-angle[inxn].pcoa4 * pow(BOjk - 1.5, 2));

                        PEcoa = angle[inxn].pcoa1 / (1.0 + exp_coa2) * exp_coa3i * exp_coa3k * exp_coa4i * exp_coa4k;

                        // PEcoa derivative part:
                        CEcoa[0] = -2.0 * angle[inxn].pcoa4 * (BOij - 1.5);          // dBOij
                        CEcoa[1] = -2.0 * angle[inxn].pcoa4 * (BOjk - 1.5);          // dBOjk
                        CEcoa[2] = -angle[inxn].pcoa2 * exp_coa2 / (1.0 + exp_coa2); // dDj
                        CEcoa[3] = -2.0 * angle[inxn].pcoa3 * (-BOij + sum_BOi);
                        CEcoa[4] = -2.0 * angle[inxn].pcoa3 * (-BOjk + sum_BOk);
                        for (std::size_t x = 0; x < 5; ++x)
                        {
                            CEcoa[x] = CEcoa[x] * PEcoa;
                        }

                        // if the j-atom is a resident count the potential energies.
                        PE[5] = PE[5] + PEval;
                        PE[6] = PE[6] + PEpen;
                        PE[7] = PE[7] + PEcoa;

                        CE3body_b[0] = CEpen[1] + CEcoa[0] - CEcoa[3] + CEval[0]; // BO_ij
                        CE3body_b[1] = CEpen[2] + CEcoa[1] - CEcoa[4] + CEval[1]; // BO_jk

                        CE3body_d[0] = CEpen[0] + CEcoa[2] + CEval[2] + CEval[6]; // delta_j
                        CE3body_d[1] = CEcoa[3];                                  // delta_i
                        CE3body_d[2] = CEcoa[4];                                  // delta_k

                        CE3body_a = CEval[7];

                        // Force calculation
                        for (std::size_t ib = 0; ib < atomData[i].bonds.size(); ++ib)
                        {
                            if (atomData[i].bonds[ib] == j)
                            {
                                ForceB(i, j, ib, CE3body_b[0]);
                                break; // BO_ij
                            }
                        }

                        for (std::size_t kb = 0; kb < atomData[k].bonds.size(); ++kb)
                        {
                            if (atomData[k].bonds[kb] == j)
                            {
                                ForceB(j, k, kb, CE3body_b[1]);
                                break; // BO_jk
                            }
                        }

                        for (std::size_t jb = 0; jb < atomData[j].bonds.size(); ++jb)
                        {
                            coeff_val = CE3body_d[0] + CEval[5] * pow(atomData[j].BO_sum[jb], 7);
                            coeff = Vector3d(coeff_val, CEval[4] + coeff_val, CEval[4] + coeff_val);

                            const std::size_t n = atomData[j].bonds[jb];
                            ForceBbo(j, n, jb, coeff);
                        }

                        atomData[i].cdbnd = atomData[i].cdbnd + CE3body_d[1];
                        atomData[k].cdbnd = atomData[k].cdbnd + CE3body_d[2];

                        ForceA3(CE3body_a, i, j, k, rij, rjk, rij0, rjk0);
                    }
                }
            }
        }
    }

    void ForceA4(const Double &coeff, const std::size_t &i, const std::size_t &j, const std::size_t &k, const std::size_t &l, const Vector3d &da0, const Vector3d &da1, const Vector3d &da2, const Double &da0_0, const Double &da1_0, const Double &da2_0)
    {
        Vector3d rij, rjk, rkl;
        Vector3d fij, fjk, fkl, fijjk, fjkkl;
        Vector3d Cwi, Cwj, Cwl;
        Containers::StaticArray<3, Vector3d> Caa;
        Vector2d Daa;
        Double DDisqr, coDD, com;

        Caa[0][0] = pow(da0_0, 2);
        Caa[0][1] = (da0 * da1).sum();
        Caa[0][2] = (da0 * da2).sum();

        Caa[1][0] = Caa[0][1];
        Caa[1][1] = pow(da1_0, 2);
        Caa[1][2] = (da1 * da2).sum();
        
        Caa[2][0] = Caa[0][2];
        Caa[2][1] = Caa[1][2];
        Caa[2][2] = pow(da2_0, 2);

        Daa[0] = Caa[0][0] * Caa[1][1] - Caa[0][1] * Caa[0][1];
        Daa[1] = Caa[1][1] * Caa[2][2] - Caa[1][2] * Caa[1][2];

        DDisqr = 1.0 / sqrt(Daa[0] * Daa[1]);
        coDD = coeff * DDisqr;

        com = Caa[1][0] * Caa[1][2] - Caa[0][2] * Caa[1][1];

        // Some of calculations are unnecessary due to the action-reaction relation.
        Cwi[0] = Caa[1][1] / Daa[0] * com;
        Cwi[1] = -(Caa[1][2] + Caa[0][1] / Daa[0] * com);
        Cwi[2] = Caa[1][1];

        Cwj[0] = -(Caa[1][2] + (Caa[1][1] + Caa[1][0]) / Daa[0] * com);
        Cwj[1] = -(Caa[1][2] - 2 * Caa[0][2] - Caa[2][2] / Daa[1] * com - (Caa[0][0] + Caa[1][0]) / Daa[0] * com);
        Cwj[2] = -(Caa[1][0] + Caa[1][1] + Caa[1][2] / Daa[1] * com);

        Cwl[0] = -Caa[1][1];
        Cwl[1] = (Caa[1][0] + Caa[2][1] / Daa[1] * com);
        Cwl[2] = -(Caa[1][1] / Daa[1] * com);

        rij = da0;
        rjk = da1;
        rkl = da2;

        fij = coDD * (Cwi[0] * rij + Cwi[1] * rjk + Cwi[2] * rkl);
        fjk = coDD * ((Cwj[0] + Cwi[0]) * rij + (Cwj[1] + Cwi[1]) * rjk + (Cwj[2] + Cwi[2]) * rkl);
        fkl = -coDD * (Cwl[0] * rij + Cwl[1] * rjk + Cwl[2] * rkl);

        fijjk = -fij + fjk;
        fjkkl = -fjk + fkl;

        atomData[i].force = atomData[i].force + fij;
        atomData[j].force = atomData[j].force + fijjk;
        atomData[k].force = atomData[k].force + fjkkl;
        atomData[l].force = atomData[l].force - fkl;
    }

    void E4b()
    {
        std::size_t itype, jtype, ktype, ltype;
        std::size_t i, k, l;
        std::size_t inxn;

        // angles
        Vector3d cos_ijkl;
        Double cos_ijkl_sqr, cos_2ijkl, sin_ijk, sin_jkl, tan_ijk_i, tan_jkl_i; // ! [unused] sin_ijkl
        Double cos_ijk, cos_jkl, theta_ijk, theta_jkl, omega_ijkl;

        // vectors
        Vector3d rjk, rij, rkl, crs_ijk, crs_jkl;
        Double rjk0, rij0, rkl0, crs_ijk0, crs_jkl0;

        Double delta_ang_j, delta_ang_k, delta_ang_jk;
        Double exp_tor1, exp_tor3, exp_tor4, exp_tor34_i, fn10, fn11, dfn11, fn12, PEtors, PEconj, cmn;
        Vector3d exp_tor2;

        // coefficents
        Containers::StaticArray<9, Double> CEtors;
        Containers::StaticArray<6, Double> CEconj;
        Vector3d C4body_a, C4body_b, C4body_b_jk;

        Double Cconj, BOij, BOjk, BOkl, btb2;

        for (std::size_t j = 0; j < NATOMS; ++j)
        {
            jtype = atomData[j].type;
            delta_ang_j = atomData[j].delta + atom[jtype].Val - atom[jtype].Valangle;

            for (std::size_t k0; k0 < atomData[j].bonds.size(); ++k0)
            {
                if (atomData[j].BO_sum[k0] <= cutof2_esub)
                {
                    BOjk = atomData[j].BO_sum[k0] - cutof2_esub;
                    k = atomData[j].bonds[k0];

                    if (j < k)
                        continue;

                    ktype = atomData[k].type;
                    delta_ang_k = atomData[k].delta + atom[ktype].Val - atom[ktype].Valangle;
                    delta_ang_jk = delta_ang_j + delta_ang_k;

                    rjk = atomData[j].position - atomData[k].position;
                    rjk0 = sqrt(rjk.dot());

                    for (std::size_t i0; i0 < atomData[j].bonds.size(); ++i0)
                    {
                        if (atomData[j].BO_sum[i0] > cutof2_esub && (atomData[j].BO_sum[i0] * atomData[j].BO_sum[i0]) > cutof2_esub)
                        {
                            BOij = atomData[j].BO_sum[i0] - cutof2_esub;
                            i = atomData[j].bonds[i0];

                            if (i == k)
                                continue;

                            itype = atomData[i].type;
                            rij = atomData[i].position - atomData[j].position;
                            rij0 = sqrt(rij.dot());

                            // Calculate the angle i-j-k
                            cos_ijk = -(rij * rjk).sum() / (rij0 * rjk0);
                            if (cos_ijk > MAXANGLE)
                                cos_ijk = MAXANGLE;
                            if (cos_ijk < MINANGLE)
                                cos_ijk = MINANGLE;

                            theta_ijk = acos(cos_ijk);
                            sin_ijk = sin(theta_ijk);
                            tan_ijk_i = 1.0 / tan(theta_ijk);

                            crs_ijk = Math::cross(rij, rjk);
                            crs_ijk0 = sqrt(crs_ijk.dot());
                            if (crs_ijk0 < NSMALL)
                                crs_ijk0 = NSMALL;

                            for (std::size_t l0; l0 < atomData[k].bonds.size(); ++l0)
                            {
                                if (atomData[k].BO_sum[l0] > cutof2_esub && (atomData[k].BO_sum[l0] * atomData[k].BO_sum[l0]) > cutof2_esub)
                                {
                                    BOkl = atomData[k].BO_sum[l0] - cutof2_esub;
                                    l = atomData[k].bonds[l0];
                                    ltype = atomData[l].type;
                                    inxn = atom[itype].inxn4[jtype][ktype][ltype];

                                    if (inxn != 0 && i != l && j != l)
                                    {
                                        inxn -= 1;
                                        // cutoff condition to ignore bonding
                                        if (pow(atomData[j].BO_sum[i0] * atomData[j].BO_sum[k0], 2) * atomData[k].BO_sum[l0] > MINBO0)
                                        {
                                            rkl = atomData[k].position - atomData[l].position;
                                            rkl0 = sqrt(rkl.dot());

                                            exp_tor2[0] = exp(-torsion[inxn].ptor2 * BOij); // i-j
                                            exp_tor2[1] = exp(-torsion[inxn].ptor2 * BOjk); // j-k
                                            exp_tor2[2] = exp(-torsion[inxn].ptor2 * BOkl); // k-l

                                            exp_tor3 = exp(-torsion[inxn].ptor3 * delta_ang_jk);
                                            exp_tor4 = exp(torsion[inxn].ptor4 * delta_ang_jk);
                                            exp_tor34_i = 1.0 / (1.0 + exp_tor3 + exp_tor4);

                                            fn10 = (1.0 - exp_tor2[0]) * (1.0 - exp_tor2[1]) * (1.0 - exp_tor2[2]);
                                            fn11 = (2.0 + exp_tor3) / (1.0 + exp_tor3 + exp_tor4);

                                            fn12 = exp(-torsion[inxn].pcot2 * (pow(BOij - 1.5, 2) + pow(BOjk - 1.5, 2) + pow(BOkl - 1.5, 2)));
                                            // pi-bond value used here is not the subtracted one but the original value
                                            btb2 = 2.0 - atomData[j].BO[k0][1] - fn11;
                                            exp_tor1 = exp(pow(torsion[inxn].ptor1 * btb2, 2));

                                            // Get angle variables i-j-k, j-k-l, i-j-k-l
                                            cos_jkl = -(rjk * rkl).sum() / (rjk0 * rkl0);
                                            if (cos_jkl > MAXANGLE)
                                                cos_jkl = MAXANGLE;
                                            if (cos_jkl < MINANGLE)
                                                cos_jkl = MINANGLE;

                                            theta_jkl = acos(cos_jkl);
                                            sin_jkl = sin(theta_jkl);
                                            tan_jkl_i = 1.0 / tan(theta_jkl);

                                            crs_jkl = Math::cross(rjk, rkl);
                                            crs_jkl0 = sqrt(crs_jkl.dot());
                                            if (crs_jkl0 < NSMALL)
                                                crs_jkl0 = NSMALL;

                                            cos_ijkl[0] = (crs_ijk * crs_jkl).sum() / (crs_ijk0 * crs_jkl0);
                                            if (cos_ijkl[0] > MAXANGLE)
                                                cos_ijkl[0] = MAXANGLE;
                                            if (cos_ijkl[0] < MINANGLE)
                                                cos_ijkl[0] = MINANGLE;

                                            omega_ijkl = acos(cos_ijkl[0]);
                                            cos_ijkl_sqr = cos_ijkl[0] * cos_ijkl[0];
                                            cos_2ijkl = cos(2.0 * omega_ijkl);
                                            cos_ijkl[1] = 1.0 - cos_2ijkl;
                                            cos_ijkl[2] = 1.0 + cos(3.0 * omega_ijkl);
                                            // ! [unused] sin_ijkl = sin(omega_ijkl);

                                            PEtors = 0.5 * fn10 * sin_ijk * sin_jkl * (torsion[inxn].V1 * (1.0 + cos_ijkl[0]) + torsion[inxn].V2 * exp_tor1 * cos_ijkl[1] + torsion[inxn].V3 * cos_ijkl[2]);
                                            PEconj = torsion[inxn].pcot1 * fn12 * (1.0 + (cos_ijkl_sqr - 1.0) * sin_ijk * sin_jkl);

                                            PE[8] = PE[8] + PEtors;
                                            PE[9] = PE[9] + PEconj;
                                            // Force coefficient calculation
                                            // TOrsional term
                                            CEtors[0] = 0.5 * sin_ijk * sin_jkl * (torsion[inxn].V1 * (1.0 + cos_ijkl[0]) + torsion[inxn].V2 * exp_tor1 * cos_ijkl[1] + torsion[inxn].V3 * cos_ijkl[2]);
                                            CEtors[1] = -torsion[inxn].ptor2 * fn10 * sin_ijk * sin_jkl * torsion[inxn].V2 * exp_tor1 * btb2 * cos_ijkl[1];

                                            dfn11 = (-torsion[inxn].ptor3 * exp_tor3 + (torsion[inxn].ptor3 * exp_tor3 - torsion[inxn].ptor4 * exp_tor4) * (2.0 + exp_tor3) * exp_tor34_i) * exp_tor34_i;

                                            CEtors[2] = CEtors[1] * dfn11;
                                            CEtors[3] = CEtors[0] * torsion[inxn].ptor2 * exp_tor2[0] * (1.0 - exp_tor2[1]) * (1.0 - exp_tor2[2]);
                                            CEtors[4] = CEtors[0] * torsion[inxn].ptor2 * (1.0 - exp_tor2[0]) * exp_tor2[1] * (1.0 - exp_tor2[2]);
                                            CEtors[5] = CEtors[0] * torsion[inxn].ptor2 * (1.0 - exp_tor2[0]) * (1.0 - exp_tor2[1]) * exp_tor2[2];

                                            cmn = -0.5 * fn10 * (torsion[inxn].V1 * (1.0 + cos_ijkl[0]) + torsion[inxn].V2 * exp_tor1 * cos_ijkl[1] + torsion[inxn].V3 * cos_ijkl[2]);

                                            CEtors[6] = cmn * sin_jkl * tan_ijk_i;
                                            CEtors[7] = cmn * sin_ijk * tan_jkl_i;
                                            CEtors[8] = fn10 * sin_ijk * sin_jkl * (0.5 * torsion[inxn].V1 - 2.0 * torsion[inxn].V2 * exp_tor1 * cos_ijkl[0] + 1.5 * torsion[inxn].V3 * (cos_2ijkl + 2.0 * cos_ijkl_sqr));

                                            // Conjugation Energy
                                            Cconj = -2.0 * torsion[inxn].pcot2 * PEconj;
                                            CEconj[0] = Cconj * (BOij - 1.5);
                                            CEconj[1] = Cconj * (BOjk - 1.5);
                                            CEconj[2] = Cconj * (BOkl - 1.5);

                                            CEconj[3] = -torsion[inxn].pcot1 * fn12 * (cos_ijkl_sqr - 1.0) * tan_ijk_i * sin_jkl;
                                            CEconj[4] = -torsion[inxn].pcot1 * fn12 * (cos_ijkl_sqr - 1.0) * sin_ijk * tan_jkl_i;
                                            CEconj[5] = 2.0 * torsion[inxn].pcot1 * fn12 * cos_ijkl[0] * sin_ijk * sin_jkl;

                                            C4body_b = Vector3d(CEconj[0], CEconj[1], CEconj[2]) + Vector3d(CEtors[3], CEtors[4], CEtors[5]); // dBOij, dBOjk, dBOkl
                                            C4body_a = Vector3d(CEconj[3], CEconj[4], CEconj[5]) + Vector3d(CEtors[6], CEtors[7], CEtors[8]); // ijk, jkl, ijkl

                                            atomData[j].cdbnd = atomData[j].cdbnd + CEtors[2];
                                            atomData[k].cdbnd = atomData[k].cdbnd + CEtors[2];

                                            for (std::size_t ib = 0; ib < atomData[i].bonds.size(); ++ib)
                                            {
                                                if (atomData[i].bonds[ib] == j)
                                                {
                                                    ForceB(i, j, ib, C4body_b[0]);
                                                    break;
                                                }
                                            }

                                            // To take care of the derivative of BOpi(j,k), add <Ctors(2)> to 
                                            // the full BOjk derivative coefficient <C4body_b(2)>, but only pi-bond component
                                            C4body_b_jk = Vector3d(C4body_b[1], CEtors[1] + C4body_b[1], C4body_b[1]);

                                            for (std::size_t kb = 0; kb < atomData[k].bonds.size(); ++kb)
                                            {
                                                if (atomData[k].bonds[kb] == j)
                                                {
                                                    ForceBbo(j, k, kb, C4body_b_jk);
                                                    break;
                                                }
                                            }

                                            for (std::size_t lb = 0; lb < atomData[l].bonds.size(); ++lb)
                                            {
                                                if (atomData[l].bonds[lb] == k)
                                                {
                                                    ForceB(k, l, lb, C4body_b[2]);
                                                    break;
                                                }
                                            }

                                            ForceA3(C4body_a[0], i, j, k, rij, rjk, rij0, rjk0);
                                            ForceA3(C4body_a[1], j, k, l, rjk, rkl, rjk0, rkl0);
                                            ForceA4(C4body_a[2], i, j, k, l, rij, rjk, rkl, rij0, rjk0, rkl0);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void ForceD(const std::size_t &i, const Double &coeff)
    {
        std::size_t j;
        Vector3d Cbond, dr, ff;

        for (std::size_t b = 0; b < atomData[i].bonds.size(); ++b)
        {
            j = atomData[i].bonds[b];
            if (j < i)
                continue;
            
            Cbond[0] = coeff * (atomData[i].A0[b] + atomData[i].BO_sum[b] * atomData[i].A1[b]);
            dr = atomData[i].position - atomData[j].position;
            ff = Cbond[0] * atomData[i].dBOp[b] * dr;

            atomData[i].force = atomData[i].force - ff;
            atomData[j].force = atomData[j].force + ff;
            
            Cbond[1] = coeff * atomData[i].BO_sum[b] * atomData[i].A2[b];
            
            for (std::size_t ib = 0; ib < atomData[j].bonds.size(); ++ib)
            {
                if (atomData[j].bonds[ib] == i)
                {
                    Cbond[2] = coeff * atomData[i].BO_sum[b] * atomData[j].A2[ib];
                    break;
                }
            }

            atomData[i].ccbnd = atomData[i].ccbnd + Cbond[1];
            atomData[j].ccbnd = atomData[i].ccbnd + Cbond[2];
        }
    }

    void ForceBondedTerms()
    {
        std::size_t j;
        Vector3d dr, ff;

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            ForceD(i, atomData[i].cdbnd);
            for (std::size_t b = 0; b < atomData[i].bonds.size(); ++b)
            {
                j = atomData[i].bonds[b];
                if (j < i)
                    continue;
                dr = atomData[i].position - atomData[j].position;
                ff = atomData[i].ccbnd * atomData[i].dBOp[b] * dr;
                atomData[i].force = atomData[i].force - ff;
                atomData[j].force = atomData[j].force + ff;
            }
            atomData[i].ccbnd = 0.0;
        }
    }

    void FORCE()
    {
        // calculate BO prime
        BOPRIM();
        // calculate full BO
        BOFULL();
        ENbond();
        Ebond();
        Elnpr();
        Ehb();
        E3b();
        E4b();
        ForceBondedTerms();

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            astr[0] = astr[0] + atomData[i].position[0] * atomData[i].force[0];
            astr[1] = astr[1] + atomData[i].position[1] * atomData[i].force[1];
            astr[2] = astr[2] + atomData[i].position[2] * atomData[i].force[2];
            astr[3] = astr[3] + atomData[i].position[1] * atomData[i].force[2];
            astr[4] = astr[4] + atomData[i].position[2] * atomData[i].force[0];
            astr[5] = astr[5] + atomData[i].position[0] * atomData[i].force[1];
        }
    }

    void vkick(const Double &dtf)
    {
        std::size_t itype;

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            itype = atomData[i].type;
            atomData[i].velocity = atomData[i].velocity + dtf * atom[itype].dthm * atomData[i].force;
        }
    }
}
