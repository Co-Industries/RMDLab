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
        INITSYSTEM();
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

        for (std::size_t i = 0; i < NATOMS; ++i)
        {
            const Vector3d tempPosition = Vector3d(std::rand(), std::rand(), std::rand()) / Double(RAND_MAX);
            const Vector3d tempVelocity = Vector3d(std::rand(), std::rand(), std::rand()) / Double(RAND_MAX);

            atomData[i].position = tempPosition * 200.0 - Vector3d(100.0);
            atomData[i].position.y() *= 0.5;
            atomData[i].velocity = (tempVelocity * 200.0 - Vector3d{100.0}).resized(randomVelocity);
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
    }

    void Simulation::INITSYSTEM()
    {
        // for QEq (simulation should be Qeq instead of PQeq)
        rctap = rctap0;
        rctap2 = pow(rctap, 2);

        // unit distance in r^2 scale
        UDR = rctap2 / NTABLE;
        UDRi = 1.0 / UDR;
        Debug{} << "UDR:" << UDR << UDRi;

        // CTap
        CTap[0] = 1.0;
        CTap[1] = 0.0;
        CTap[2] = 0.0;
        CTap[3] = 0.0;
        CTap[4] = -35.0 / pow(rctap, 4);
        CTap[5] = 84.0 / pow(rctap, 5);
        CTap[6] = -70.0 / pow(rctap, 6);
        CTap[7] = 20.0 / pow(rctap, 7);
    }

    void Simulation::UPDATE_ATOMS()
    {
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
        //
        for (std::size_t i = 0; i < NATOMS; i++)
        {
            atomInstanceData[i].color = Color3(1.0, 0.0, 0.0);
            arrayResize(atomData[i].hessian, 0);
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
                if (dpq2 < rctap2)
                {
                    atomInstanceData[i].color = Color3(1.0f, 0.5f + atomInstanceData[i].color.y(), 1.0f);
                    atomInstanceData[j].color = Color3(1.0f, 0.5f + atomInstanceData[j].color.y(), 1.0f);

                    // TODO
                    // array with atype[atom] + inxn(atype, atype);
                    // ity (i type), jty (j type)

                    // get table value from distance
                    std::size_t itb = Int(dpq2 * UDRi);
                    Float drtb = dpq2 - itb * UDR;
                    drtb = drtb * UDRi;

                    // inxn = inxn2[atomData[i].type]
                    // Debug{} << drtb;

                    // arrayAppend(atomData[i].hessian, InPlaceInit, (1.0 - drtb));
                }
            }
        }
    }
}