/* A B C D E F G H I J K L M N O P Q R S T U V W X Y Z */

#include <Magnum/Math/Functions.h>

#include <Magnum/MeshTools/Compile.h>

#include <Magnum/Primitives/Cube.h>
#include <Magnum/Primitives/Icosphere.h>

#include <Magnum/Trade/MeshData.h>

#include "../components/scene/Scene.h"
#include "Functions.h"
#include "Simulation.h"

namespace Magnum
{
    Simulation::Simulation(Scene3D &_scene, SceneGraph::DrawableGroup3D &_drawables)
    {
        Debug{} << "Simulation has loaded";
        atomInstanceData = Containers::Array<AtomInstanceData>{};
        atomData = Containers::Array<AtomData>{};

        atomShader = Shaders::PhongGL{Shaders::PhongGL::Configuration{}
                                          .setFlags(Shaders::PhongGL::Flag::VertexColor |
                                                    Shaders::PhongGL::Flag::InstancedTransformation)};
        atomShader.setShininess(5.0f);
        atomShader.setSpecularColor(Color3({0.1f}));
        atomInstanceBuffer = GL::Buffer{};
        atomMesh = MeshTools::compile(Primitives::icosphereSolid(2));
        atomMesh.addVertexBufferInstanced(atomInstanceBuffer, 1, 0,
                                          Shaders::PhongGL::TransformationMatrix{},
                                          Shaders::PhongGL::NormalMatrix{},
                                          Shaders::PhongGL::Color3{});
        atomInstanceBuffer.setData(atomInstanceData, GL::BufferUsage::DynamicDraw);
        atomMesh.setInstanceCount(atomInstanceData.size());
        auto atomObject = new Object3D{&_scene};
        new AtomDrawable{*atomObject, atomShader, atomMesh, _drawables};

        //? OCTREE
        octree.emplace(Vector3{0}, 100.0f, rctap0);
        octreeShader = Shaders::FlatGL3D{Shaders::FlatGL3D::Configuration{}
                                             .setFlags(Shaders::FlatGL3D::Flag::VertexColor |
                                                       Shaders::FlatGL3D::Flag::InstancedTransformation)};
        octreeInstanceBuffer = GL::Buffer{};
        octreeMesh = MeshTools::compile(Primitives::cubeWireframe());
        octreeMesh.addVertexBufferInstanced(octreeInstanceBuffer, 1, 0,
                                            Shaders::FlatGL3D::TransformationMatrix{},
                                            Shaders::FlatGL3D::Color3{});
        octreeInstanceBuffer.setData(octreeInstanceData, GL::BufferUsage::DynamicDraw);
        octreeMesh.setInstanceCount(octreeInstanceData.size());
        auto octreeObject = new Object3D{&_scene};
        new FlatGLDrawable{*octreeObject, octreeShader, octreeMesh, _drawables};
    }

    void Simulation::RUN(const SimulationParameters &parameters)
    {
        //? Parameters
        NATOMS = parameters.atomCount;
        atomRadius = parameters.atomRadius;
        randomVelocity = parameters.randomVelocity;

        Debug{} << "Simulation running" << NATOMS;
        arrayResize(atomInstanceData, NATOMS);
        arrayResize(atomData, NATOMS);
        arrayResize(atomFloatPositions, NATOMS);

        GETPARAMS();
        INITSYSTEM();

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            const Vector3d tempPosition = Vector3d(std::rand(), std::rand(), std::rand()) / Double(RAND_MAX);
            const Vector3d tempVelocity = Vector3d(std::rand(), std::rand(), std::rand()) / Double(RAND_MAX);
            const std::size_t tempType = std::rand() / ((RAND_MAX + 1u) / 2); // 0 or 1

            atomData[i].position = tempPosition * 200.0 - Vector3d(100.0);
            atomData[i].position.y() *= 0.5;
            atomData[i].velocity = (tempVelocity * 200.0 - Vector3d{100.0}).resized(randomVelocity);
            atomData[i].type = tempType;
            atomFloatPositions[i] = Vector3(atomData[i].position);

            atomInstanceData[i].transformationMatrix = Matrix4::translation(atomFloatPositions[i]) * Matrix4::scaling(Vector3{atomRadius});
            atomInstanceData[i].normalMatrix = atomInstanceData[i].transformationMatrix.normalMatrix();
            atomInstanceData[i].color = Color3(1.0f, 0.0f, 0.0f);
        }

        atomMesh.setInstanceCount(atomInstanceData.size());
        octree->setPoints(atomFloatPositions);
        octree->build();
        Debug{} << " [Octree] Allocated nodes:" << octree->numAllocatedNodes();
        Debug{} << " [Octree] Max number of points per node:" << octree->maxNumPointInNodes();

        // QeQ
        // arrayResize(qsfp, NATOMS);
        // arrayResize(qsfv, NATOMS);
        // arrayResize(q, NATOMS); //! Should initialize per atom

        // Atomtype (should be given)
    }

    void Simulation::GETPARAMS()
    {
        // arrayResize(TBL_Eclmb_QEq, 4); // +1 from last inxn2
        for (std::size_t i = 0; i < NTABLE; ++i)
        {
            arrayResize(TBL_Eclmb_QEq[i], nso);
        }
        for (std::size_t i = 0; i < nso; ++i)
        {
            arrayResize(atom[i].gamW, nso);
            arrayResize(atom[i].gamij, nso);
            arrayResize(atom[i].inxn2, nso);
            arrayResize(atom[i].r0s, nso);
            arrayResize(atom[i].r0p, nso);
            arrayResize(atom[i].r0pp, nso);
        }

        atom[0].inxn2[0] = 1;
        atom[1].inxn2[1] = 2;
        atom[0].inxn2[1] = 3;
        atom[1].inxn2[0] = 3;

        // atom
        Containers::StaticArray<3, std::string> _name{"H", "O", "X"};
        Containers::StaticArray<3, Double> _vop{33.2894, 11.7301, 2.5};
        Containers::StaticArray<3, Double> _gam{0.82, 1.095, 1.0};
        Containers::StaticArray<3, Double> _eta{9.6093, 1.0548, -0.1};
        Containers::StaticArray<3, Double> _chi{3.7248, 8.5, 8.50};
        Containers::StaticArray<3, Double> _rat{0.893, 1.245, -0.1};
        Containers::StaticArray<3, Double> _rapt{-0.1, 1.0548, -0.1};
        Containers::StaticArray<3, Double> _vnq{-0.1, 0.9049, -0.1};
        Containers::StaticArray<3, Double> _Val{1.0, 2.0, 2.0};
        Containers::StaticArray<3, Double> _Valval{1.0, 4.0, 2.0};
        Containers::StaticArray<3, Double> _bo131{3.0408, 3.5357, 8.741};
        Containers::StaticArray<3, Double> _bo132{2.4197, 0.6653, 13.364};
        Containers::StaticArray<3, Double> _bo133{0.0003, 0.0021, 0.669};

        // bond
        Containers::StaticArray<3, Double> _pbo1{-0.079, -0.1225, -0.0924};
        Containers::StaticArray<3, Double> _pbo2{6.0552, 5.50000, 4.27780};
        Containers::StaticArray<3, Double> _pbo3{1.0000, -0.1055, 1.00000};
        Containers::StaticArray<3, Double> _pbo4{0.0000, 9.00000, 0.00000};
        Containers::StaticArray<3, Double> _pbo5{1.0000, 1.00000, 1.00000};
        Containers::StaticArray<3, Double> _pbo6{6.0000, 29.7503, 6.00000};
        Containers::StaticArray<3, Double> _ovc{0.00000, 1.00000, 0.00000};
        Containers::StaticArray<3, Double> _v13cor{1.0000, 1.0000, 1.0000};

        // ? change atom type values
        for (std::size_t i = 0; i < nso; ++i)
        {
            atom[i].name = _name[i];
            atom[i].vop = _vop[i];
            atom[i].gam = _gam[i];
            atom[i].eta = _eta[i];
            atom[i].chi = _chi[i];
            atom[i].rat = _rat[i];
            atom[i].rapt = _rapt[i];
            atom[i].vnq = _vnq[i];
            atom[i].Val = _Val[i];
            atom[i].Valval = _Valval[i];
            atom[i].bo131 = _bo131[i];
            atom[i].bo132 = _bo132[i];
            atom[i].bo133 = _bo133[i];
        }

        for (std::size_t i = 0; i < nboty; ++i)
        {
            bond[i].pbo1 = _pbo1[i];
            bond[i].pbo2 = _pbo2[i];
            bond[i].pbo3 = _pbo3[i];
            bond[i].pbo4 = _pbo4[i];
            bond[i].pbo5 = _pbo5[i];
            bond[i].pbo6 = _pbo6[i];
            bond[i].ovc = _ovc[i];
            bond[i].v13cor = _v13cor[i];
        }

        for (std::size_t i = 0; i < nso; ++i)
        {
            for (std::size_t j = 0; j < nso; ++j)
            {
                // Terms for the Bond Order Calculation:
                atom[i].r0s[j] = 0.5 * (atom[i].rat + atom[j].rat);
                atom[i].r0p[j] = 0.5 * (atom[i].rapt + atom[j].rapt);
                atom[i].r0pp[j] = 0.5 * (atom[i].vnq + atom[j].vnq);

                // Terms used in van der Waals calc:
                atom[i].gamW[j] = sqrt(atom[i].vop * atom[j].vop);
                atom[i].gamij[j] = pow((atom[i].gam * atom[j].gam), -1.5);

                // ? Bonds
                std::size_t inxn = atom[i].inxn2[j]; // -1, because values start from 1 (rxmd)
                if (inxn == 0)
                    continue;
                inxn -= 1;

                // ! In BOp calculation, <swh> will be multiplied to <BOp> to remove
                // ! BOpi and BOpipi for bonding interaction of atoms with a hydrogen.
                if (atom[i].rat > 0.0 && atom[j].rat > 0.0)
                    bond[inxn].swh[0] = 1.0;
                if (atom[i].rapt > 0.0 && atom[j].rapt > 0.0)
                    bond[inxn].swh[1] = 1.0;
                if (atom[i].vnq > 0.0 && atom[j].vnq > 0.0)
                    bond[inxn].swh[2] = 1.0;

                if (atom[i].r0s[j] <= 0.0)
                    bond[inxn].cBOp1 = 0.0;
                else
                    bond[inxn].cBOp1 = bond[inxn].pbo1 / pow((atom[i].r0s[j]), bond[inxn].pbo2);
                if (atom[i].r0p[j] <= 0.0)
                    bond[inxn].cBOp3 = 0.0;
                else
                    bond[inxn].cBOp3 = bond[inxn].pbo3 / pow((atom[i].r0p[j]), bond[inxn].pbo4);
                if (atom[i].r0pp[j] <= 0.0)
                    bond[inxn].cBOp5 = 0.0;
                else
                    bond[inxn].cBOp5 = bond[inxn].pbo5 / pow((atom[i].r0pp[j]), bond[inxn].pbo6);

                bond[inxn].pbo2h = 0.5 * bond[inxn].pbo2;
                bond[inxn].pbo4h = 0.5 * bond[inxn].pbo4;
                bond[inxn].pbo6h = 0.5 * bond[inxn].pbo6;

                bond[inxn].pboc1 = 50.0; // ! from ffield (same as vpar1 and vpar2)
                bond[inxn].pboc2 = 9.5469;
                bond[inxn].pboc3 = sqrt(atom[i].bo132 * atom[j].bo132);
                bond[inxn].pboc4 = sqrt(atom[i].bo131 * atom[j].bo131);
                bond[inxn].pboc5 = sqrt(atom[i].bo133 * atom[j].bo133);
            }
        }

        // for (std::size_t i; i < nboty; i++)
        //{
        //     read type a, read type b (1 1, 2 2, 1 2)
        //     inxn2[type a][type b] = i;
        // }
    }

    void Simulation::INITSYSTEM()
    {
        // ? CUTOFFLENGTH()
        cutoff_vpar30 = cutof2_bo * vpar30;

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

        // unit distance in r^2 scale
        UDR = rctap2 / NTABLE;
        UDRi = 1.0f / UDR;
        Debug{} << "UDR:" << UDR << UDRi;

        for (std::size_t i = 0; i < nso; ++i)
        {
            for (std::size_t j = 0; j < nso; ++j)
            {
                std::size_t inxn = atom[i].inxn2[j];
                // inxn2 doesnt use 0, thats the "null" value, any table using it should have +1 value
                if (inxn == 0)
                    continue;
                // Double gamWij = gamW[i][j];
                inxn -= 1;

                for (std::size_t k = 0; k < NTABLE; ++k)
                {
                    Double dr2 = UDR * k;
                    Double dr1 = sqrt(dr2);

                    Double dr3 = dr1 * dr2;
                    Double dr4 = dr2 * dr2;
                    Double dr5 = dr1 * dr2 * dr2;
                    Double dr6 = dr2 * dr2 * dr2;
                    Double dr7 = dr1 * dr2 * dr2 * dr2;

                    Double Tap = CTap[7] * dr7 + CTap[6] * dr6 + CTap[5] * dr5 + CTap[4] * dr4 + CTap[0];
                    Double dr3gamij = pow((dr3 + atom[i].gamij[j]), -1.0 / 3.0);

                    TBL_Eclmb_QEq[k][inxn] = Tap * Cclmb0_qeq * dr3gamij;
                }

                // ? CUTOFFLENGTH()
                Double dr = 1.0;
                Double BOsig = 1.0;
                while (BOsig > MINBOSIG)
                {
                    dr += 0.01;
                    BOsig = exp(pow(bond[inxn].pbo1 * (dr / atom[i].r0s[j]), bond[inxn].pbo2));
                }
                bond[inxn].rc = dr;
                bond[inxn].rc2 = dr * dr;
            }
        }
    }

    void Simulation::UPDATE_ATOMS()
    {
        QEq();
        // FORCE();
        constexpr Float dt = 1.0 / 1.2;

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            Vector3d pos = atomData[i].position + atomData[i].velocity * dt;
            for (std::size_t j = 0; j < 3; ++j)
            {
                if (pos[j] < -100.0 || pos[j] > 100.0)
                    atomData[i].velocity[j] = -atomData[i].velocity[j];
                pos[j] = Math::clamp(pos[j], -100.0, 100.0);
            }

            atomData[i].position = pos;
            atomFloatPositions[i] = Vector3(pos);
            atomInstanceData[i].transformationMatrix.translation() = atomFloatPositions[i];
        }
        atomInstanceBuffer.setData(atomInstanceData, GL::BufferUsage::DynamicDraw);
    }

    void Simulation::UPDATE_OCTREE()
    {
        for (std::size_t i = 0; i < NATOMS; i++)
        {
            atomInstanceData[i].color = Color3(1.0, 0.0, 0.0);
            arrayResize(atomData[i].hessian, 0);
            arrayResize(atomData[i].neighbors, 0);
            arrayResize(atomData[i].bonds, 0);
            arrayResize(atomData[i].bo, 0);
            arrayResize(atomData[i].dln_BOp, 0);
            arrayResize(atomData[i].dBOp, 0);
            arrayResize(atomData[i].bo_sum, 0);
            arrayResize(atomData[i].A0, 0);
            arrayResize(atomData[i].A1, 0);
            arrayResize(atomData[i].A2, 0);
            arrayResize(atomData[i].A3, 0);
        }
        //
        octreeCollisionDetection();
        octree->update();
        arrayResize(octreeInstanceData, 0);
        arrayAppend(octreeInstanceData, InPlaceInit,
                    Matrix4::translation(octree->center()) *
                        Matrix4::scaling(Vector3{octree->halfWidth()}),
                    0x00ffff_rgbf);

        if (drawOctreeBounds)
        {

            const auto &activeTreeNodeBlocks = octree->activeTreeNodeBlocks();
            for (OctreeNodeBlock *const pNodeBlock : activeTreeNodeBlocks)
            {
                for (std::size_t childIdx = 0; childIdx < 8; ++childIdx)
                {
                    const OctreeNode &pNode = pNodeBlock->_nodes[childIdx];

                    /* Non-empty node */
                    if (!pNode.isLeaf() || pNode.pointCount() > 0)
                    {
                        /*INFO missing _arcballCamera->viewMatrix()*/
                        const Matrix4 t = Matrix4::translation(pNode.center()) *
                                          Matrix4::scaling(Vector3{pNode.halfWidth()});
                        arrayAppend(octreeInstanceData, InPlaceInit, t, 0x197f99_rgbf);
                    }
                }
            }
        }
        octreeInstanceBuffer.setData(octreeInstanceData, GL::BufferUsage::DynamicDraw);
        octreeMesh.setInstanceCount(octreeInstanceData.size());
    }

    void Simulation::octreeCollisionDetection()
    {
        const OctreeNode &rootNode = octree->rootNode();
        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            checkCollisionWithSubTree(rootNode, i,
                                      atomData[i].position, atomData[i].velocity,
                                      Range3D::fromCenter(atomFloatPositions[i], Vector3{rctap}));
        }
    }

    void Simulation::checkCollisionWithSubTree(
        const OctreeNode &node,
        std::size_t i,
        const Vector3d &ppos,
        const Vector3d &pvel,
        const Range3D &bounds)
    {
        if (!node.looselyOverlaps(bounds))
            return;

        if (!node.isLeaf())
        {
            for (std::size_t childIdx = 0; childIdx < 8; ++childIdx)
            {
                const OctreeNode &child = node.childNode(childIdx);
                checkCollisionWithSubTree(child, i, ppos, pvel, bounds);
            }
        }

        for (const OctreePoint *const point : node.pointList())
        {
            const std::size_t j = point->idx();
            if (j > i)
            {
                // const Vector3d qpos = atomData[j].position;
                //  const Vector3 qvel = _atomVelocities[j];
                //  const Vector3 velpq = pvel - qvel;
                const Vector3d pospq = ppos - atomData[j].position; // qpos
                // const Float vp = Math::dot(velpq, pospq);
                /* INFO Velocity vector, to stop attraction*/
                // if (vp < 0.0f)
                // const Float dpq = pospq.length();
                const Float dpq2 = (pospq * 0.5).dot();

                // TODO copyptr(6) probably uses all atoms, so rewrite everything :)
                if (dpq2 < rctap2)
                {
                    atomInstanceData[i].color = Color3(1.0f, 0.5f + atomInstanceData[i].color.y(), 1.0f);
                    atomInstanceData[j].color = Color3(1.0f, 0.5f + atomInstanceData[j].color.y(), 1.0f);

                    // TODO
                    // array with atype[atom] + inxn(atype, atype);
                    // ity (i type), jty (j type)

                    // get table value from distance
                    std::size_t itb = Int(dpq2 * UDRi);
                    Double drtb = dpq2 - itb * UDR;
                    drtb = drtb * Double(UDRi);

                    std::size_t inxn = atom[atomData[i].type].inxn2[atomData[j].type] - 1;

                    arrayAppend(atomData[i].hessian, InPlaceInit, (1.0 - drtb) * TBL_Eclmb_QEq[itb][inxn] + drtb * TBL_Eclmb_QEq[itb][inxn]); // TODO could be itb + 1
                    arrayAppend(atomData[i].neighbors, InPlaceInit, j);

                    arrayAppend(atomData[j].hessian, InPlaceInit, (1.0 - drtb) * TBL_Eclmb_QEq[itb][inxn] + drtb * TBL_Eclmb_QEq[itb][inxn]);
                    arrayAppend(atomData[j].neighbors, InPlaceInit, i);

                    //? BOPRIM
                    atomData[i].deltap[0] = -atom[atomData[i].type].Val;
                    atomData[j].deltap[0] = -atom[atomData[j].type].Val;

                    Vector3d _bo, _dln_BOp;
                    Double _bo_sum;

                    if (dpq2 < Float(bond[inxn].rc2))
                    {
                        arrayAppend(atomData[i].bonds, InPlaceInit, j);
                        arrayAppend(atomData[j].bonds, InPlaceInit, i);

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
                            arrayAppend(atomData[i].dBOp, InPlaceInit, 0.0);
                            arrayAppend(atomData[j].dBOp, InPlaceInit, 0.0);
                            arrayAppend(atomData[i].bo_sum, InPlaceInit, 0.0);
                            arrayAppend(atomData[j].bo_sum, InPlaceInit, 0.0);
                        }
                    }

                    //? BOFULL
                    atomData[i].deltap[1] = atomData[i].deltap[0] + atom[atomData[i].type].Val - atom[atomData[i].type].Valval;
                    atomData[j].deltap[1] = atomData[j].deltap[0] + atom[atomData[j].type].Val - atom[atomData[j].type].Valval;

                    Double exppboc1i = exp(-vpar1 * atomData[i].deltap[0]);
                    Double exppboc2i = exp(-vpar2 * atomData[i].deltap[0]);
                    Double exppboc1j = exp(-vpar1 * atomData[j].deltap[0]);
                    Double exppboc2j = exp(-vpar2 * atomData[j].deltap[0]);

                    Double fn2 = exppboc1i + exppboc1j;
                    Double fn3 = (-1.0 / vpar2) * log(0.5 * (exppboc2i + exppboc2j));

                    Double fn23 = fn2 + fn3;
                    // TODO Double BOp0 = atomData[i].bo_sum[some index idk];
                    Double BOp0 = 1.0;

                    Double fn1 = 0.5 * ((atom[atomData[i].type].Val + fn2) / (atom[atomData[i].type].Val + fn23) + (atom[atomData[j].type].Val + fn2) / (atom[atomData[j].type].Val + fn23));

                    // TODO ovc is either 1 or 0, so this doesn't make any sense. Probably some scaling i have missed;
                    if (bond[inxn].ovc < 1.0e-3)
                        fn1 = 1.0;

                    Double BOpsqr = BOp0 * BOp0;
                    Double fn4 = 1.0 / (1.0 + exp(-bond[inxn].pboc3 * (bond[inxn].pboc4 * BOpsqr - atomData[i].deltap[1])) + bond[inxn].pboc5);
                    Double fn5 = 1.0 / (1.0 + exp(-bond[inxn].pboc3 * (bond[inxn].pboc4 * BOpsqr - atomData[j].deltap[1])) + bond[inxn].pboc5);

                    if (bond[inxn].v13cor < 1.0e-3)
                    {
                        fn4 = 1.0;
                        fn5 = 1.0;
                    }

                    Double fn45 = fn4 * fn5;
                    Double fn145 = fn1 * fn45;
                    Double fn1145 = fn1 * fn145;

                    // New bond order definition
                    Vector3d _BO;
                    Double _BO_sum;
                    _BO_sum = _bo_sum * fn145;
                    _BO[1] = _bo[1] * fn1145;
                    _BO[2] = _bo[2] * fn1145;
                    if (_BO_sum < 1.0e-10)
                        _BO_sum = 0.0;
                    if (_BO[1] < 1.0e-10)
                        _BO[1] = 0.0;
                    if (_BO[2] < 1.0e-10)
                        _BO[2] = 0.0;

                    // new sigma BO definition
                    _BO[0] = _bo[0] - _BO[1] - _BO[2];
                    // TODO assign BO to atom;
                    // _BO(i) = _BO(j)

                    // CALCULATION OF DERIVATIVE OF BOND ORDER
                    // all following comes from Coding Methodology section:
                    // part 1:

                    Double u1ij = atom[atomData[i].type].Val + fn23;
                    Double u1ji = atom[atomData[j].type].Val + fn23;

                    // part 2:
                    Double u1ij_inv2 = 1.0 / (u1ij * u1ij);
                    Double u1ji_inv2 = 1.0 / (u1ji * u1ji);

                    Double Cf1Aij = 0.5 * fn3 * (u1ij_inv2 + u1ji_inv2);
                    Double Cf1Bij = -0.5 * ((u1ij - fn3) * u1ij_inv2 + (u1ji - fn3) * u1ji_inv2);

                    // part 3:
                    Double exp_delt22 = exppboc2i + exppboc2j;
                    Double Cf1ij = (-Cf1Aij * bond[inxn].pboc1 * exppboc1i) + (Cf1Bij * exppboc2i) / (exp_delt22);
                    Double Cf1ji = (-Cf1Aij * bond[inxn].pboc1 * exppboc1j) + (Cf1Bij * exppboc2j) / (exp_delt22);

                    // part 4:
                    Double pboc34 = bond[inxn].pboc3 * bond[inxn].pboc4;
                    Double BOpij_2 = BOpsqr;

                    Double u45ij = bond[inxn].pboc5 + bond[inxn].pboc3 * atomData[i].deltap[1] - pboc34 * BOpij_2;
                    Double u45ji = bond[inxn].pboc5 + bond[inxn].pboc3 * atomData[j].deltap[1] - pboc34 * BOpij_2;

                    // part 5:
                    Double exph_45ij = exp(u45ij);
                    Double exph_45ji = exp(u45ji);
                    Double exp1 = 1.0 / (1.0 + exph_45ij);
                    Double exp2 = 1.0 / (1.0 + exph_45ji);
                    Double exp12 = exp1 * exp2;

                    Double Cf45ij = -exph_45ij * exp12 * exp1;
                    Double Cf45ji = -exph_45ji * exp12 * exp2;

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
                    Double fn45_inv = 1.0 / fn45;
                    Double Cf1ij_div1 = Cf1ij / fn1;
                    Double Cf1ji_div1 = Cf1ji / fn1;

                    arrayAppend(atomData[i].A0, InPlaceInit, fn145);
                    arrayAppend(atomData[i].A1, InPlaceInit, -2.0 * pboc34 * BOp0 * (Cf45ij + Cf45ji) * fn45_inv);
                    arrayAppend(atomData[i].A2, InPlaceInit, Cf1ij_div1 + (bond[inxn].pboc3 * Cf45ij * fn45_inv));
                    arrayAppend(atomData[i].A3, InPlaceInit, atomData[i].A2[atomData[i].A2.size() - 1] + Cf1ij_div1);

                    arrayAppend(atomData[j].A0, InPlaceInit, atomData[i].A0[atomData[i].A0.size() - 1]);
                    arrayAppend(atomData[j].A1, InPlaceInit, atomData[i].A1[atomData[i].A1.size() - 1]);
                    arrayAppend(atomData[j].A2, InPlaceInit, Cf1ji_div1 + (bond[inxn].pboc3 * Cf45ji * fn45_inv));
                    arrayAppend(atomData[j].A3, InPlaceInit, atomData[i].A2[atomData[i].A2.size() - 1] + Cf1ji_div1);
                }

                // atomData[i].delta = -atom[atomData[i].type].Val + sum(all BO?) tired... sleep....
            }
        }
    }
}