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
        running = true;
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
            //atomData[i].position.y() *= 0.5;
            atomData[i].velocity = (tempVelocity * 200.0 - Vector3d{100.0}).resized(randomVelocity);
            atomData[i].type = tempType;
            atomFloatPositions[i] = Vector3(atomData[i].position);

            atomInstanceData[i].transformationMatrix = Matrix4::translation(atomFloatPositions[i]) * Matrix4::scaling(Vector3{atomRadius});
            atomInstanceData[i].normalMatrix = atomInstanceData[i].transformationMatrix.normalMatrix();
            atomInstanceData[i].color = atom[atomData[i].type].color;
        }

        atomMesh.setInstanceCount(atomInstanceData.size());
        octree->setPoints(atomFloatPositions);
        octree->build();
        Debug{} << " [Octree] Allocated nodes:" << octree->numAllocatedNodes();
        Debug{} << " [Octree] Max number of points per node:" << octree->maxNumPointInNodes();

        UPDATE_OCTREE();

        QEq();
        FORCE();
        
    }

    void Simulation::GETPARAMS()
    {
        // arrayResize(TBL_Eclmb_QEq, 4); // +1 from last inxn2
        for (std::size_t i = 0; i < NTABLE + 1; ++i)
        {
            arrayResize(TBL_Eclmb_QEq[i], nso);
            arrayResize(TBL_Eclmb_p[i], nso);
            arrayResize(TBL_Eclmb_d[i], nso);
            arrayResize(TBL_Evdw_p[i], nso);
            arrayResize(TBL_Evdw_d[i], nso);
        }
        for (std::size_t i = 0; i < nso; ++i)
        {
            arrayResize(atom[i].gamW, nso);
            arrayResize(atom[i].gamij, nso);
            arrayResize(atom[i].inxn2, nso);
            arrayResize(atom[i].r0s, nso);
            arrayResize(atom[i].r0p, nso);
            arrayResize(atom[i].r0pp, nso);
            arrayResize(atom[i].Dij, nso);
            arrayResize(atom[i].alpij, nso);
            arrayResize(atom[i].rvdW, nso);
        }

        atom[0].inxn2[0] = 1;
        atom[1].inxn2[1] = 2;
        atom[0].inxn2[1] = 3;
        atom[1].inxn2[0] = 3;

        // atom
        Containers::StaticArray<3, std::string> _name{"H", "O", "X"};
        Containers::StaticArray<3, Color3> _color{Color3::fromSrgbInt(0xffffff), Color3::fromSrgbInt(0xff0d0d), Color3::fromSrgbInt(0x000000)};
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
        Containers::StaticArray<3, Double> _eps{0.093, 0.1038, 0.0};
        Containers::StaticArray<3, Double> _alf{8.218, 9.7942, 10.0};
        Containers::StaticArray<3, Double> _rvdw1{1.355, 2.3808, 2.0};

        // bond
        Containers::StaticArray<3, Double> _pbo1{-0.079, -0.1225, -0.0924};
        Containers::StaticArray<3, Double> _pbo2{6.0552, 5.5, 4.2778};
        Containers::StaticArray<3, Double> _pbo3{1.0, -0.1055, 1.0};
        Containers::StaticArray<3, Double> _pbo4{0.0, 9.0, 0.0};
        Containers::StaticArray<3, Double> _pbo5{1.0, 1.0, 1.0};
        Containers::StaticArray<3, Double> _pbo6{6.0, 29.7503, 6.0};
        Containers::StaticArray<3, Double> _ovc{0.0, 1.0, 0.0};
        Containers::StaticArray<3, Double> _v13cor{1.0, 1.0, 1.0};
        Containers::StaticArray<3, Double> _pbe1{-0.46, 0.2506, -0.577};
        Containers::StaticArray<3, Double> _pbe2{6.25, 0.3451, 1.1413};
        Containers::StaticArray<3, Double> _Desig{153.3934, 142.2858, 167.2086};
        Containers::StaticArray<3, Double> _Depi{0.0, 145.0, 0.0};
        Containers::StaticArray<3, Double> _Depipi{0.0, 50.8293, 0.0};

        // ? change atom type values
        for (std::size_t i = 0; i < nso; ++i)
        {
            atom[i].name = _name[i];
            atom[i].color = _color[i];
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
            atom[i].eps = _eps[i];
            atom[i].alf = _alf[i];
            atom[i].rvdw1 = _rvdw1[i];
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
            bond[i].pbe1 = _pbe1[i];
            bond[i].pbe2 = _pbe2[i];
            bond[i].Desig = _Desig[i];
            bond[i].Depi = _Depi[i];
            bond[i].Depipi = _Depipi[i];
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
                atom[i].rvdW[j] = sqrt(4.0 * atom[i].rvdw1 * atom[j].rvdw1);
                atom[i].Dij[j] = sqrt(atom[i].eps * atom[j].eps);
                atom[i].alpij[j] = sqrt(atom[i].alf * atom[j].alf);
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
    }

    void Simulation::INITSYSTEM()
    {
        // TODO declare variables outside loop for performance (not that important since this runs only once)

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

        // van der Waals
        pvdW1h = 0.5 * pvdW1;
        pvdW1inv = 1.0 / pvdW1;

        for (std::size_t i = 0; i < nso; ++i)
        {
            for (std::size_t j = 0; j < nso; ++j)
            {
                std::size_t inxn = atom[i].inxn2[j];
                // * inxn2 doesnt use 0, thats the "null" value, any table using it should have +1 value
                if (inxn == 0)
                    continue;
                inxn -= 1;

                const Double gamWij = atom[i].gamW[j];
                const Double gamij = atom[i].gamij[j];
                const Double Dij0 = atom[i].Dij[j];
                const Double alphaij = atom[i].alpij[j];
                const Double rvdW0 = atom[i].rvdW[j];

                const Double gamwinvp = pow(1.0 / gamWij, pvdW1);

                for (std::size_t k = 0; k < NTABLE + 1; ++k)
                {
                    Double dr2 = UDR * k;
                    Double dr1 = sqrt(dr2);

                    Double dr3 = dr1 * dr2;
                    Double dr4 = dr2 * dr2;
                    Double dr5 = dr1 * dr2 * dr2;
                    Double dr6 = dr2 * dr2 * dr2;
                    Double dr7 = dr1 * dr2 * dr2 * dr2;

                    Double Tap = CTap[7] * dr7 + CTap[6] * dr6 + CTap[5] * dr5 + CTap[4] * dr4 + CTap[0];

                    Double rij_vd1 = pow(dr2, pvdW1h);
                    Double fn13 = pow(rij_vd1 + gamwinvp, pvdW1inv);
                    Double exp1 = exp(alphaij * (1.0 - fn13 / rvdW0));
                    Double exp2 = sqrt(exp1);

                    Double dr3gamij = pow(dr3 + gamij, -1.0 / 3.0);

                    TBL_Evdw_p[k][inxn] = Tap * Dij0 * (exp1 - 2.0 * exp2);
                    TBL_Eclmb_p[k][inxn] = Tap * Cclmb0 * dr3gamij;
                    TBL_Eclmb_QEq[k][inxn] = Tap * Cclmb0_qeq * dr3gamij;
                    
                    //Force calculation:
                    Double dTap = 7.0 * CTap[7] * dr5 + 6.0 * CTap[6] * dr4 + 5.0 * CTap[5] * dr3 + 4.0 * CTap[4] * dr2;
                    Double dfn13 = pow(rij_vd1 + gamwinvp, pvdW1inv - 1.0) * pow(dr2, pvdW1h - 1.0);

                    TBL_Evdw_d[k][inxn] = Dij0 * (dTap * (exp1 - 2.0 * exp2) - Tap * (alphaij / rvdW0) * (exp1 - exp2) * dfn13);
                    TBL_Eclmb_d[k][inxn] = Cclmb0 * dr3gamij * (dTap - pow(dr3gamij, 3) * Tap * dr1);

                    //TODO if(isLG) then
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
        FORCE();
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
            // TODO should just allocate array for every atom, more memory but better performance
            atomInstanceData[i].color = atom[atomData[i].type].color;
            arrayResize(atomData[i].hessian, 0);
            arrayResize(atomData[i].dpq2, 0);
            arrayResize(atomData[i].neighbors, 0);
            arrayResize(atomData[i].bonds, 0);
            arrayResize(atomData[i].bo, 0);
            arrayResize(atomData[i].BO, 0);
            arrayResize(atomData[i].dln_BOp, 0);
            arrayResize(atomData[i].dBOp, 0);
            arrayResize(atomData[i].bo_sum, 0);
            arrayResize(atomData[i].BO_sum, 0);
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

        // TODO declare variables outside loop for performance
        for (const OctreePoint *const point : node.pointList())
        {
            const std::size_t j = point->idx();
            if (j > i)
            {
                const Vector3d pospq = ppos - atomData[j].position; // qpos
                const Float dpq2 = (pospq * 0.5).dot();

                if (dpq2 < rctap2)
                {
                    arrayAppend(atomData[i].neighbors, InPlaceInit, j);
                    arrayAppend(atomData[j].neighbors, InPlaceInit, i);

                    arrayAppend(atomData[i].dpq2, InPlaceInit, dpq2);
                    arrayAppend(atomData[j].dpq2, InPlaceInit, dpq2);

                    atomInstanceData[i].color.g() = atomInstanceData[i].color.g() + 0.4f;
                    atomInstanceData[j].color.g() = atomInstanceData[j].color.g() + 0.4f;

                    const std::size_t inxn = atom[atomData[i].type].inxn2[atomData[j].type] - 1;
                    if (dpq2 < Float(bond[inxn].rc2))
                    {
                        arrayAppend(atomData[i].bonds, InPlaceInit, j);
                        arrayAppend(atomData[j].bonds, InPlaceInit, i);
                    }
                }
            }
        }
    }
}