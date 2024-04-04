/* A B C D E F G H I J K L M N O P Q R S T U V W X Y Z */

#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>

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
        // TODO declare variables outside loop for performance
        Double gssum = 0.0;
        Double gtsum = 0.0;

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
                const std::size_t j = atomData[i].neighbors[n];

                // ? qeq_initialize()
                if (j < i)
                    break;

                const std::size_t itb = Int(atomData[i].dpq2[n] * Double(UDRi));

                const std::size_t itb1 = itb + 1;
                Double drtb = atomData[i].dpq2[n] - itb * Double(UDR);

                drtb = drtb * Double(UDRi);

                const std::size_t inxn = atom[atomData[i].type].inxn2[atomData[j].type] - 1;
                arrayAppend(atomData[i].hessian, InPlaceInit, (1.0 - drtb) * TBL_Eclmb_QEq[itb][inxn] + drtb * TBL_Eclmb_QEq[itb1][inxn]);
                arrayAppend(atomData[j].hessian, InPlaceInit, (1.0 - drtb) * TBL_Eclmb_QEq[itb][inxn] + drtb * TBL_Eclmb_QEq[itb1][inxn]);

                // ? get_gradient [Gnew]

                gssum += atomData[i].hessian[n] * atomData[j].qs;
                gtsum += atomData[i].hessian[n] * atomData[j].qt;

                // TODO might need to do atomData[j] - idk
            }
            const std::size_t _type = atomData[i].type;
            const Double _eta = atom[_type].eta;
            atomData[i].gs = atom[_type].chi - _eta * atomData[i].qs - gssum;
            atomData[i].gt = 1.0 - _eta * atomData[i].qt - gtsum;

            atomData[i].hs = atomData[i].gs;
            atomData[i].ht = atomData[i].gt;
        }

        Gnew[0] = dot_product(1);
        Gnew[1] = dot_product(2);

        Double GEst2 = 1.0e99;

        for (std::size_t nstep_qeq = 0; nstep_qeq < NMAXQEq; ++nstep_qeq)
        {
            // ? get_hsh()
            Double Est = 0.0, hshs_sum = 0.0, hsht_sum = 0.0;

            for (std::size_t i = 0; i < NATOMS; ++i)
            {
                const std::size_t _type = atomData[i].type;
                const Double _eta = atom[_type].eta;
                Double _hshs = _eta * atomData[i].hs;
                Double _hsht = _eta * atomData[i].ht;
                Est += atom[_type].chi * atomData[i].q + 0.5 * _eta * atomData[i].q * atomData[i].q;

                for (std::size_t n = 0; n < atomData[i].neighbors.size(); ++n)
                {
                    const std::size_t j = atomData[i].neighbors[n];

                    if (j < i)
                        break;

                    _hshs += atomData[i].hessian[n] * atomData[j].hs;
                    _hsht += atomData[i].hessian[n] * atomData[j].ht;
                    // get half of potential energy, then sum it up if atoms are resident
                    // Double Est1 = 0.5 * atomData[i].hessian[n] * atomData[i].q * atomData[j].q;
                    Double Est1 = atomData[i].hessian[n] * atomData[i].q * atomData[j].q;
                    Est += Est1;
                    // Est += Est1;
                }

                hshs_sum += _hshs * atomData[i].hs;
                hsht_sum += _hsht * atomData[i].ht;
            }

            // ? QEq
            Double GEst1 = Est;
            if (0.5 * (abs(GEst2) + abs(GEst1)) < QEq_tol)
                break;
            if (abs(GEst2) > 0.0 && (abs(GEst1 / GEst2 - 1.0) < QEq_tol))
                break;
            GEst2 = GEst1;
            // line minimization factor of <s> vector
            Vector2d g_h, h_hsh;
            g_h[0] = dot_product(3);
            h_hsh[0] = hshs_sum;

            // line minimization factor of <t> vector
            g_h[1] = dot_product(4);
            h_hsh[1] = hsht_sum;

            // TODO MPI_ALLREDUCE does some buffer stuff, might be important

            Vector2d lmin = g_h / h_hsh;
            Double ssum = 0.0, tsum = 0.0;

            for (std::size_t i = 0; i < NATOMS; ++i)
            {
                atomData[i].qs += lmin[0] * atomData[i].hs;
                atomData[i].qt += lmin[1] * atomData[i].ht;
                ssum += atomData[i].qs;
                tsum += atomData[i].qt;
            }

            const Double mu = ssum / tsum;

            gssum = 0.0;
            gtsum = 0.0;

            for (std::size_t i = 0; i < NATOMS; ++i)
            {
                atomData[i].q = atomData[i].qs - mu * atomData[i].qt;
                for (std::size_t n = 0; n < atomData[i].neighbors.size(); ++n)
                {
                    const std::size_t j = atomData[i].neighbors[n];

                    // ? get_gradient [Gnew]
                    if (j < i)
                        break;

                    gssum += atomData[i].hessian[n] * atomData[j].qs;
                    gtsum += atomData[i].hessian[n] * atomData[j].qt;
                }
                const std::size_t _type = atomData[i].type;
                const Double _eta = atom[_type].eta;
                atomData[i].gs = atom[_type].chi - _eta * atomData[i].qs - gssum;
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
        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            atomData[i].deltap[0] = -atom[atomData[i].type].Val;
        }

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            for (std::size_t n = 0; n < atomData[i].neighbors.size(); ++n)
            {
                const std::size_t j = atomData[i].neighbors[n];

                if (j < i)
                    break;

                const std::size_t inxn = atom[atomData[i].type].inxn2[atomData[j].type] - 1;
                const Float dpq2 = atomData[i].dpq2[n];

                Vector3d _bo, _dln_BOp;
                Double _bo_sum;

                for (std::size_t b = 0; b < atomData[i].bonds.size(); ++b)
                {
                    if (atomData[i].bonds[b] != j)
                        continue;

                    arg_BOpij[0] = bond[inxn].cBOp1 * pow(Double(dpq2), bond[inxn].pbo2h);
                    arg_BOpij[1] = bond[inxn].cBOp3 * pow(Double(dpq2), bond[inxn].pbo4h);
                    arg_BOpij[2] = bond[inxn].cBOp5 * pow(Double(dpq2), bond[inxn].pbo6h);

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
                        _dln_BOp = _dln_BOp / dpq2;
                        arrayAppend(atomData[i].dln_BOp, InPlaceInit, _dln_BOp);
                        arrayAppend(atomData[j].dln_BOp, InPlaceInit, _dln_BOp);

                        Double dBOp = Vector3d(_bo * _dln_BOp).sum();
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
    }

    void BOFULL()
    {
        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            atomData[i].deltap[1] = atomData[i].deltap[0] + atom[atomData[i].type].Val - atom[atomData[i].type].Valval;
        }

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            const Double exppboc1i = exp(-vpar1 * atomData[i].deltap[0]);
            const Double exppboc2i = exp(-vpar2 * atomData[i].deltap[0]);

            for (std::size_t n = 0; n < atomData[i].neighbors.size(); ++n)
            {
                const std::size_t j = atomData[i].neighbors[n];

                if (j < i)
                    break;

                const Double exppboc1j = exp(-vpar1 * atomData[j].deltap[0]);
                const Double exppboc2j = exp(-vpar2 * atomData[j].deltap[0]);

                const std::size_t inxn = atom[atomData[i].type].inxn2[atomData[j].type] - 1;

                const Double fn2 = exppboc1i + exppboc1j;
                const Double fn3 = (-1.0 / vpar2) * log(0.5 * (exppboc2i + exppboc2j));
                const Double fn23 = fn2 + fn3;

                for (std::size_t b = 0; b < atomData[i].bonds.size(); ++b)
                {

                    if (atomData[i].bonds[b] != j)
                        continue;

                    const Double BOp0 = atomData[i].bo_sum[b];

                    Double fn1 = 0.5 * ((atom[atomData[i].type].Val + fn2) / (atom[atomData[i].type].Val + fn23) + (atom[atomData[j].type].Val + fn2) / (atom[atomData[j].type].Val + fn23));
                    // TODO ovc is either 1 or 0, so this doesn't make any sense. Probably some scaling i have missed or Fortran might not have booleans
                    if (bond[inxn].ovc < 1.0e-3)
                        fn1 = 1.0;

                    const Double BOpsqr = BOp0 * BOp0;
                    Double fn4 = 1.0 / (1.0 + exp(-bond[inxn].pboc3 * (bond[inxn].pboc4 * BOpsqr - atomData[i].deltap[1])) + bond[inxn].pboc5);
                    Double fn5 = 1.0 / (1.0 + exp(-bond[inxn].pboc3 * (bond[inxn].pboc4 * BOpsqr - atomData[j].deltap[1])) + bond[inxn].pboc5);
                    if (bond[inxn].v13cor < 1.0e-3)
                    {
                        fn4 = 1.0;
                        fn5 = 1.0;
                    }

                    const Double fn45 = fn4 * fn5;
                    const Double fn145 = fn1 * fn45;
                    const Double fn1145 = fn1 * fn145;

                    // New bond order definition
                    Vector3d _BO;
                    Double _BO_sum;
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

                    const Double u1ij = atom[atomData[i].type].Val + fn23;
                    const Double u1ji = atom[atomData[j].type].Val + fn23;

                    // part 2:
                    const Double u1ij_inv2 = 1.0 / (u1ij * u1ij);
                    const Double u1ji_inv2 = 1.0 / (u1ji * u1ji);

                    const Double Cf1Aij = 0.5 * fn3 * (u1ij_inv2 + u1ji_inv2);
                    const Double Cf1Bij = -0.5 * ((u1ij - fn3) * u1ij_inv2 + (u1ji - fn3) * u1ji_inv2);

                    // part 3:
                    const Double exp_delt22 = exppboc2i + exppboc2j;
                    Double Cf1ij = (-Cf1Aij * bond[inxn].pboc1 * exppboc1i) + (Cf1Bij * exppboc2i) / (exp_delt22);
                    Double Cf1ji = (-Cf1Aij * bond[inxn].pboc1 * exppboc1j) + (Cf1Bij * exppboc2j) / (exp_delt22);

                    // part 4:
                    const Double pboc34 = bond[inxn].pboc3 * bond[inxn].pboc4;
                    const Double BOpij_2 = BOpsqr;

                    const Double u45ij = bond[inxn].pboc5 + bond[inxn].pboc3 * atomData[i].deltap[1] - pboc34 * BOpij_2;
                    const Double u45ji = bond[inxn].pboc5 + bond[inxn].pboc3 * atomData[j].deltap[1] - pboc34 * BOpij_2;

                    // part 5:
                    const Double exph_45ij = exp(u45ij);
                    const Double exph_45ji = exp(u45ji);
                    const Double exp1 = 1.0 / (1.0 + exph_45ij);
                    const Double exp2 = 1.0 / (1.0 + exph_45ji);
                    const Double exp12 = exp1 * exp2;

                    Double Cf45ij = -exph_45ij * exp12 * exp1;
                    Double Cf45ji = -exph_45ji * exp12 * exp2;

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
                    const Double fn45_inv = 1.0 / fn45;
                    const Double Cf1ij_div1 = Cf1ij / fn1;
                    const Double Cf1ji_div1 = Cf1ji / fn1;

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
        }

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            Double _BO_nsum = 0.0;
            for (std::size_t b = 0; b < atomData[i].BO_sum.size(); ++b)
            {
                _BO_nsum += atomData[i].BO_sum[b];
            }
            atomData[i].delta = -atom[atomData[i].type].Val - atom[atomData[i].type].Valval + _BO_nsum;
        }
    }

    void ENbond()
    {
        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            PE[13] = PE[13] + CEchrge * (atom[atomData[i].type].chi * atomData[i].q + pow(0.5 * atom[atomData[i].type].eta * atomData[i].q, 2));
            for (std::size_t n = 0; n < atomData[i].neighbors.size(); ++n)
            {
                const std::size_t j = atomData[i].neighbors[n];

                if (j < i)
                    break;

                const std::size_t itb = Int(atomData[i].dpq2[n] * Double(UDRi));
                const std::size_t itb1 = itb + 1;
                Double drtb = atomData[i].dpq2[n] - itb * Double(UDR);
                drtb = drtb * Double(UDRi);

                const std::size_t inxn = atom[atomData[i].type].inxn2[atomData[j].type] - 1;

                // van der Waals: TBL_Evdw[0] = TBL_Evdw_p ; TBL_Evdw[1] = TBL_Evdw_d
                const Double PEvdw = (1.0 - drtb) * TBL_Evdw_p[itb][inxn] + drtb * TBL_Evdw_p[itb1][inxn];
                const Double CEvdw = (1.0 - drtb) * TBL_Evdw_d[itb][inxn] + drtb * TBL_Evdw_d[itb1][inxn];

                // Coulomb:
                const Double qij = atomData[i].q * atomData[j].q;
                Double PEclmb = (1.0 - drtb) * TBL_Eclmb_p[itb][inxn] + drtb * TBL_Eclmb_p[itb1][inxn];
                PEclmb *= qij;
                Double CEclmb = (1.0 - drtb) * TBL_Eclmb_d[itb][inxn] + drtb * TBL_Eclmb_d[itb1][inxn];
                CEclmb *= qij;

                PE[11] = PE[11] + PEvdw;
                PE[12] = PE[12] + PEclmb;

                Vector3d dr = atomData[i].position - atomData[j].position;
                Vector3d ff = (CEvdw + CEclmb) * dr;

                atomData[i].force = atomData[i].force - ff;
                atomData[j].force = atomData[j].force + ff;
            }
        }
    }

    void ForceBbo(const std::size_t &i, const std::size_t &j, const std::size_t &b, const Vector3d &coeff)
    {
        // * With the new bond-order definition, 1st term is the derivative of "full"-bond order,
        // * 2nd is for pi-bond order and 3rd is for pipi-bond order.
        Vector3d cf = Vector3d(coeff[0], coeff[1] - coeff[0], coeff[2] - coeff[0]);
        Vector3d Cbond;
        Cbond[0] = cf[0] * (atomData[i].A0[b] + atomData[i].BO_sum[b] * atomData[i].A1[b]) * atomData[i].dBOp[b] +
                   cf[1] * atomData[i].BO[b][1] * (atomData[i].dln_BOp[b][1] + atomData[i].A1[b] * atomData[i].dBOp[b]) +
                   cf[2] * atomData[i].BO[b][2] * (atomData[i].dln_BOp[b][2] + atomData[i].A1[b] * atomData[i].dBOp[b]);
        
        Vector3d dr = atomData[i].position - atomData[j].position;
        Vector3d ff = Cbond[0] * dr;

        atomData[i].force = atomData[i].force - ff;
        atomData[j].force = atomData[j].force + ff;

        // * 1st element is "full"-bond order.
        Vector3d cBO = Vector3d(cf[0] * atomData[i].BO_sum[b], cf[1] * atomData[i].BO[b][1], cf[2] * atomData[i].BO[b][2]);

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
        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            for (std::size_t n = 0; n < atomData[i].neighbors.size(); ++n)
            {
                const std::size_t j = atomData[i].neighbors[n];

                if (j < i)
                    break;

                const std::size_t inxn = atom[atomData[i].type].inxn2[atomData[j].type] - 1;
                for (std::size_t b = 0; b < atomData[i].bonds.size(); ++b)
                {
                    if (atomData[i].bonds[b] != j)
                        continue;
                    // TODO might be wrong... maybe (probably only calculate bond order for nearest neighbor)

                    const Double exp_be12 = exp(bond[inxn].pbe1 * (1.0 - pow(atomData[i].BO[b][0], bond[inxn].pbe2)));
                    const Double PEbo = -bond[inxn].Desig * atomData[i].BO[b][0] * exp_be12 - bond[inxn].Depi * atomData[i].BO[b][1] - bond[inxn].Depipi * atomData[i].BO[b][2];

                    PE[1] = PE[1] + PEbo;

                    const Double CEbo = -bond[inxn].Desig * exp_be12 * (1.0 - bond[inxn].pbe1 * bond[inxn].pbe2 * pow(atomData[i].BO[b][0], bond[inxn].pbe2));
                    Vector3d coeff = Vector3d(CEbo, -bond[inxn].Depi, -bond[inxn].Depipi);

                    ForceBbo(i, j, b, coeff);
                }
            }
        }
    }

    void Elnpr()
    {
        // TODO proper way of writing these functions here
        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            std::size_t _type = atomData[i].type;

            if (_type == 0)
                continue;
            
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
    }
}
