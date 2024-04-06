/* A B C D E F G H I J K L M N O P Q R S T U V W X Y Z */

#include <Magnum/Math/Functions.h>
#include <Magnum/Math/Constants.h>

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
            // atomData[i].position.y() *= 0.5;
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
            arrayResize(atom[i].inxn3hb, nso);
            arrayResize(atom[i].inxn3, nso);
            arrayResize(atom[i].inxn4, nso);

            for (std::size_t j = 0; j < nso; ++j)
            {
                arrayResize(atom[i].inxn3hb[j], nso);
                arrayResize(atom[i].inxn3[j], nso);
                arrayResize(atom[i].inxn4[j], nso);

                for (std::size_t k = 0; k < nso; ++k)
                {
                    arrayResize(atom[i].inxn4[j][k], nso);
                }
            }
        }

        // TODO change to loop when there is a input system
        atom[0].inxn2[0] = 1;
        atom[1].inxn2[1] = 2;

        atom[0].inxn2[1] = 3;
        atom[1].inxn2[0] = 3;

        atom[1].inxn3hb[0][1] = 1;

        atom[0].inxn3[0][0] = 1;
        atom[1].inxn3[1][1] = 2;

        atom[0].inxn3[1][1] = 3;
        atom[1].inxn3[1][0] = 3;

        atom[0].inxn3[1][0] = 4;
        atom[1].inxn3[0][1] = 5;

        atom[0].inxn3[0][1] = 6;
        atom[1].inxn3[0][0] = 6;

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
        Containers::StaticArray<3, Double> _Vale{1.0, 6.0, 6.0};
        Containers::StaticArray<3, Double> _mass{1.008, 15.999, 1.0};
        Containers::StaticArray<3, Double> _plp2{0.0, 0.1, 0.0};
        Containers::StaticArray<3, Double> _povun5{0.0, 37.5, 0.0};
        Containers::StaticArray<3, Double> _Valboc{1.0, 4.0, 4.0};
        Containers::StaticArray<3, Double> _pval3{4.2733, 2.7952, 2.7466};

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
        Containers::StaticArray<3, Double> _povun1{0.73, 0.6051, 0.6019};
        Containers::StaticArray<3, Double> _povun2{-19.4571, -3.6039, -11.0};

        // h_bond
        Containers::StaticArray<1, Double> _phb1{-3.6983};
        Containers::StaticArray<1, Double> _phb2{1.7831};
        Containers::StaticArray<1, Double> _phb3{17.0964};
        Containers::StaticArray<1, Double> _r0hb{2.1653};

        // angle
        Containers::StaticArray<6, Double> _theta00{0.0, 80.7324, 75.6935, 85.1864, 0.0, 0.0};
        Containers::StaticArray<6, Double> _pval1{27.9213, 30.4554, 50.0, 8.5843, 11.8475, 6.4269};
        Containers::StaticArray<6, Double> _pval2{5.8635, 0.9953, 2.0, 2.2985, 2.7571, 2.85};
        Containers::StaticArray<6, Double> _pcoa1{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        Containers::StaticArray<6, Double> _pval7{0.0, 1.631, 1.0, 2.9142, 0.0, 0.0};
        Containers::StaticArray<6, Double> _ppen1{0.0, 50.0, 0.0, 0.0, 0.0, 0.0};
        Containers::StaticArray<6, Double> _pval4{1.04, 1.0783, 1.168, 2.0521, 2.9, 1.0772};

        // torsions
        Containers::StaticArray<6, Double> _i1{1, 1, 2, 0, 0, 0};
        Containers::StaticArray<6, Double> _i2{2, 2, 2, 1, 1, 2};
        Containers::StaticArray<6, Double> _i3{2, 2, 2, 1, 2, 2};
        Containers::StaticArray<6, Double> _i4{1, 2, 2, 0, 0, 0};
        Containers::StaticArray<6, Double> _V1{2.5, 0.8302, -2.5, 0.0, 0.0, 0.5511};
        Containers::StaticArray<6, Double> _V2{-4.0, -4.0, -4.0, 0.0, 0.1, 25.415};
        Containers::StaticArray<6, Double> _V3{0.9, -0.7763, 1.0, 0.0, 0.02, 1.133};
        Containers::StaticArray<6, Double> _ptor1{-2.5, -2.5, -2.5, 0.0, -2.5415, -5.1903};
        Containers::StaticArray<6, Double> _pcot1{-1.0, -1.0, -1.0, 0.0, 0.0, -1.0};

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
            atom[i].Vale = _Vale[i];
            atom[i].plp1 = 6.0891; // vpar(16)
            atom[i].plp2 = _plp2[i];
            atom[i].nlpopt = 0.5 * (_Vale[i] - _Val[i]);
            atom[i].mass = _mass[i];
            atom[i].povun2 = _povun2[i];
            atom[i].povun3 = 50.0;   // vpar(33)
            atom[i].povun4 = 0.6991; // vpar(32)
            atom[i].povun5 = _povun5[i];
            atom[i].povun6 = 1.0588;  // vpar(7)
            atom[i].povun7 = 12.1176; // vpar(9)
            atom[i].povun8 = 13.3056; // vpar(10)
            atom[i].Valboc = _Valboc[i];

            if (atom[i].mass < 21.0 && atom[i].Valboc != atom[i].Valval)
                atom[i].Valboc = atom[i].Valval;
            atom[i].Valangle = atom[i].Valboc;

            atom[i].pval3 = _pval3[i];
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
            bond[i].povun1 = _povun1[i];
        }

        for (std::size_t i = 0; i < nhbty; ++i)
        {
            h_bond[i].phb1 = _phb1[i];
            h_bond[i].phb2 = _phb2[i];
            h_bond[i].phb3 = _phb3[i];
            h_bond[i].r0hb = _r0hb[i];
        }

        for (std::size_t i = 0; i < nvaty; ++i)
        {
            angle[i].theta00 = _theta00[i];
            angle[i].pval1 = _pval1[i];
            angle[i].pval2 = _pval2[i];
            angle[i].pcoa1 = _pcoa1[i];
            angle[i].pval7 = _pval7[i];
            angle[i].ppen1 = _ppen1[i];
            angle[i].pval4 = _pval4[i];
            // Terms that don't depend on atom types
            angle[i].pval6 = 33.8667; // vpar(15)
            angle[i].pval8 = 1.8512;  // vpar(34)
            angle[i].pval9 = 1.0563;  // vpar(17)
            angle[i].pval10 = 2.0384; // vpar(18)
            angle[i].ppen2 = 6.929;   // vpar(20)
            angle[i].ppen3 = 0.3989;  // vpar(21)
            angle[i].ppen4 = 3.9954;  // vpar(22)
            angle[i].pcoa2 = 26.5405; // vpar(3)
            angle[i].pcoa3 = 2.6962;  // vpar(39)
            angle[i].pcoa4 = 2.1365;  // vpar(31)
            // convert from degrees to radians:
            angle[i].theta00 = (Constantsd::pi() / 180.0) * angle[i].theta00;
        }

        for (std::size_t i = 0; i < ntoty; ++i)
        {
            std::size_t i1 = _i1[i];
            std::size_t i2 = _i2[i];
            std::size_t i3 = _i3[i];
            std::size_t i4 = _i4[i];
            if (i1 == 0)
            {
                for (i1 = 0; i1 < nso; ++i1)
                {
                    for (i4 = 0; i4 < nso; ++i4)
                    {
                        if (atom[i1].inxn4[i2][i3][i4] == 0 && atom[i1].inxn4[i3][i2][i4] == 0)
                        {
                            atom[i1].inxn4[i2][i3][i4] = i + 1;
                            atom[i4].inxn4[i2][i3][i1] = i + 1;
                            atom[i1].inxn4[i3][i2][i4] = i + 1;
                            atom[i4].inxn4[i3][i2][i1] = i + 1;
                        }
                    }
                }
            }
            else
            {
                atom[i1].inxn4[i2][i3][i4] = i + 1;
                atom[i4].inxn4[i2][i3][i1] = i + 1;
                atom[i1].inxn4[i3][i2][i4] = i + 1;
                atom[i4].inxn4[i3][i2][i1] = i + 1;
            }

            torsion[i].ptor1 = _ptor1[i];
            torsion[i].ptor2 = 5.7796; // vpar(24)
            torsion[i].ptor3 = 10.0;   // vpar(25)
            torsion[i].ptor4 = 1.9487; // vpar(26)
            torsion[i].V1 = _V1[i];
            torsion[i].V2 = _V2[i];
            torsion[i].V3 = _V3[i];
            torsion[i].pcot1 = _pcot1[i];
            torsion[i].pcot2 = 2.1645; // vpar(28)
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

        // *time unit conversion from [fs] -> time unit
        dt = dt / UTIME;

        // *square the spring const in the extended Lagrangian method
        Lex_w2 = 2.0 * Lex_k / dt / dt;

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

        std::size_t inxn;

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            const std::size_t itype = atomData[i].type;
            atom[itype].natoms_per_type = atom[itype].natoms_per_type + 1;
        }

        for (std::size_t i = 0; i < nso; ++i)
        {
            atom[i].dthm = dt * 0.5 / atom[i].mass;
            atom[i].hmas = 0.5 * atom[i].mass;

            for (std::size_t j = 0; j < nso; ++j)
            {
                inxn = atom[i].inxn2[j];
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

                    // Force calculation:
                    Double dTap = 7.0 * CTap[7] * dr5 + 6.0 * CTap[6] * dr4 + 5.0 * CTap[5] * dr3 + 4.0 * CTap[4] * dr2;
                    Double dfn13 = pow(rij_vd1 + gamwinvp, pvdW1inv - 1.0) * pow(dr2, pvdW1h - 1.0);

                    TBL_Evdw_d[k][inxn] = Dij0 * (dTap * (exp1 - 2.0 * exp2) - Tap * (alphaij / rvdW0) * (exp1 - exp2) * dfn13);
                    TBL_Eclmb_d[k][inxn] = Cclmb0 * dr3gamij * (dTap - pow(dr3gamij, 3) * Tap * dr1);

                    // TODO if(isLG) then
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

        for (std::size_t itype = 0; itype < nso; ++itype)
        {
            if (atom[itype].natoms_per_type == 0)
            {
                for (std::size_t jtype = 0; jtype < nso; ++jtype)
                {
                    inxn = atom[itype].inxn2[jtype];
                    if (inxn != 0)
                        bond[inxn - 1].rc = 0.0;
                    inxn = atom[jtype].inxn2[itype];
                    if (inxn != 0)
                        bond[inxn - 1].rc = 0.0;
                }
            }
        }
    }

    void Simulation::UPDATE_ATOMS()
    {
        std::size_t itype;

        vkick(1.0);
        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            atomData[i].qsfv = atomData[i].qsfv + 0.5 * dt * Lex_w2 * atomData[i].q - atomData[i].qsfp;
            atomData[i].qsfp = atomData[i].qsfp + dt * atomData[i].qsfv;
        }
        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            Vector3d pos = atomData[i].position + atomData[i].velocity * dt;
            for (std::size_t j = 0; j < 3; ++j)
            {
                if (pos[j] < -BORDER)
                {
                    pos[j] = pos[j] + BORDER2;
                }
                else if (pos[j] > BORDER)
                {
                    pos[j] = pos[j] - BORDER2;
                }
            }

            atomData[i].position = pos;
            atomFloatPositions[i] = Vector3(pos);
            atomInstanceData[i].transformationMatrix.translation() = atomFloatPositions[i];
        }
        atomInstanceBuffer.setData(atomInstanceData, GL::BufferUsage::DynamicDraw);

        if (nstep % qstep == 0)
            QEq();
        FORCE();

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            itype = atomData[i].type;
            astr[0] = astr[0] + atomData[i].velocity[0] * atomData[i].velocity[0] * atom[itype].mass;
            astr[1] = astr[1] + atomData[i].velocity[1] * atomData[i].velocity[1] * atom[itype].mass;
            astr[2] = astr[2] + atomData[i].velocity[2] * atomData[i].velocity[2] * atom[itype].mass;
            astr[3] = astr[3] + atomData[i].velocity[1] * atomData[i].velocity[2] * atom[itype].mass;
            astr[4] = astr[4] + atomData[i].velocity[2] * atomData[i].velocity[0] * atom[itype].mass;
            astr[5] = astr[5] + atomData[i].velocity[0] * atomData[i].velocity[1] * atom[itype].mass;
        }
        //vkick(1.0);

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            atomData[i].qsfv = atomData[i].qsfv = 0.5 * dt * Lex_w2 * (atomData[i].q - atomData[i].qsfp);
        }
        ++nstep;
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
                const std::size_t inxn = atom[atomData[i].type].inxn2[atomData[j].type] - 1;

                if (dpq2 < rctap2)
                {
                    arrayAppend(atomData[i].neighbors, InPlaceInit, j);
                    arrayAppend(atomData[j].neighbors, InPlaceInit, i);

                    arrayAppend(atomData[i].dpq2, InPlaceInit, dpq2);
                    arrayAppend(atomData[j].dpq2, InPlaceInit, dpq2);

                    atomInstanceData[i].color.g() = atomInstanceData[i].color.g() + 0.4f;
                    atomInstanceData[j].color.g() = atomInstanceData[j].color.g() + 0.4f;
                }
                if (dpq2 < Float(bond[inxn].rc2))
                {
                    arrayAppend(atomData[i].bonds, InPlaceInit, j);
                    arrayAppend(atomData[j].bonds, InPlaceInit, i);
                }
            }
        }
    }
}