/* A B C D E F G H I J K L M N O P Q R S T U V W X Y Z */
#include <Magnum/MeshTools/Compile.h>

#include <Magnum/Primitives/Cube.h>
#include <Magnum/Primitives/Icosphere.h>

#include <Magnum/Trade/MeshData.h>

#include "Simulation_old.h"
#include "../components/scene/Scene.h"

namespace Magnum
{
    using namespace Math::Literals;

    SimulationOld::SimulationOld(Scene3D &scene,
                                 SceneGraph::DrawableGroup3D &drawables,
                                 UnsignedInt atomCount, bool &drawOctreeBounds) : _scene(scene),
                                                                                  _drawables(drawables),
                                                                                  _atomCount(atomCount),
                                                                                  _drawOctreeBounds(drawOctreeBounds)
    {
        Debug{} << "HI!";
        _atomRadius = 0.02f;
        _atomVelocity = 0.5f;

        _atomPositions = Containers::Array<Vector3>{NoInit, _atomCount};
        _atomVelocities = Containers::Array<Vector3>{NoInit, _atomCount};
        _atomInstanceData = Containers::Array<AtomInstanceDataOld>{NoInit, _atomCount};

        for (std::size_t i = 0; i < _atomCount; ++i)
        {
            const Vector3 tmpPos = Vector3(std::rand(), std::rand(), std::rand()) / Float(RAND_MAX);
            const Vector3 tmpVel = Vector3(std::rand(), std::rand(), std::rand()) / Float(RAND_MAX);

            _atomPositions[i] = tmpPos * 2.0f - Vector3{1.0f};
            _atomPositions[i].y() *= 0.5f;
            _atomVelocities[i] = (tmpVel * 2.0f - Vector3{1.0f}).resized(_atomVelocity);

            _atomInstanceData[i].transformationMatrix =
                Matrix4::translation(_atomPositions[i]) *
                Matrix4::scaling(Vector3{_atomRadius});
            _atomInstanceData[i].normalMatrix =
                _atomInstanceData[i].transformationMatrix.normalMatrix();
            _atomInstanceData[i].color = Color3(1.0f, 0.0f, 0.0f);
        }

        _atomShader = Shaders::PhongGL{Shaders::PhongGL::Configuration{}
                                           .setFlags(Shaders::PhongGL::Flag::VertexColor |
                                                     Shaders::PhongGL::Flag::InstancedTransformation)};
        _atomShader.setShininess(10.0f);
        _atomShader.setSpecularColor(Color3({0.3}));
        _atomInstanceBuffer = GL::Buffer{};
        _atomMesh = MeshTools::compile(Primitives::icosphereSolid(2));
        _atomMesh.addVertexBufferInstanced(_atomInstanceBuffer, 1, 0,
                                           Shaders::PhongGL::TransformationMatrix{},
                                           Shaders::PhongGL::NormalMatrix{},
                                           Shaders::PhongGL::Color3{});
        _atomInstanceBuffer.setData(_atomInstanceData, GL::BufferUsage::DynamicDraw); /*prevent crash*/
        _atomMesh.setInstanceCount(_atomInstanceData.size());
        auto atomObject = new Object3D{&_scene};
        new AtomDrawable{*atomObject, _atomShader, _atomMesh, _drawables};
        /*Setup octree*/
        _octree.emplace(Vector3{0}, 1.0f, Math::max(_atomRadius, 0.1f));

        _octree->setPoints(_atomPositions);
        _octree->build();
        Debug{} << " [Octree] Allocated nodes:" << _octree->numAllocatedNodes();
        Debug{} << " [Octree] Max number of points per node:" << _octree->maxNumPointInNodes();

        _octreeShader = Shaders::FlatGL3D{Shaders::FlatGL3D::Configuration{}
                                              .setFlags(Shaders::FlatGL3D::Flag::VertexColor |
                                                        Shaders::FlatGL3D::Flag::InstancedTransformation)};
        _octreeInstanceBuffer = GL::Buffer{};
        _octreeMesh = MeshTools::compile(Primitives::cubeWireframe());
        _octreeMesh.addVertexBufferInstanced(_octreeInstanceBuffer, 1, 0,
                                             Shaders::FlatGL3D::TransformationMatrix{},
                                             Shaders::FlatGL3D::Color3{});
        /*prevent crash*/
        _octreeInstanceBuffer.setData(_octreeInstanceData, GL::BufferUsage::DynamicDraw);
        _octreeMesh.setInstanceCount(_octreeInstanceData.size());
        auto octreeObject = new Object3D{&_scene};
        new FlatGLDrawable{*octreeObject, _octreeShader, _octreeMesh, _drawables};
    }

    void SimulationOld::octreeCollisionDetection()
    {
        const OctreeNode &rootNode = _octree->rootNode();
        for (std::size_t i = 0; i < _atomPositions.size(); ++i)
        {
            checkCollisionWithSubTree(rootNode, i,
                                      _atomPositions[i], _atomVelocities[i],
                                      Range3D::fromCenter(_atomPositions[i], Vector3{_atomRadius}));
        }
    }

    bool SimulationOld::findCollision(const Containers::Array<Int> &collisions, Int value)
    {
        for (Int collision : collisions)
        {
            if (collision == value)
            {
                return true;
            }
        }
        return false;
    }

    void SimulationOld::removeCollision(Containers::Array<Int> &collisions, Int value)
    {
        arrayResize(collisions, 0);
        for (Int collision : collisions)
        {
            if (collision == value)
            {
                continue;
            }
            arrayAppend(collisions, InPlaceInit, collision);
        }
    }

    void SimulationOld::checkCollisionWithSubTree(const OctreeNode &node,
                                                  std::size_t i, const Vector3 &ppos, const Vector3 &pvel, const Range3D &bounds)
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
                const Vector3 qpos = _atomPositions[j];
                const Vector3 qvel = _atomVelocities[j];
                const Vector3 velpq = pvel - qvel;
                const Vector3 pospq = ppos - qpos;
                const Float vp = Math::dot(velpq, pospq);
                /* INFO Velocity vector, to stop attraction*/
                if (vp < 0.0f)
                {
                    const Float dpq = pospq.length();
                    if (dpq < 2.0f * _atomRadius)
                    {

                        /* INFO - Collisions*/
                        const Vector3 vNormal = vp * pospq / (dpq * dpq);
                        _atomVelocities[i] = (_atomVelocities[i] - vNormal).resized(_atomVelocity);
                        _atomVelocities[j] = (_atomVelocities[j] + vNormal).resized(_atomVelocity);
                    }
                }
            }
        }
    }

    void SimulationOld::updateAtoms()
    {
        constexpr Float dt = 1.0f / 120.0f;

        for (std::size_t i = 0; i < _atomPositions.size(); ++i)
        {
            Vector3 pos = _atomPositions[i] + _atomVelocities[i] * dt;
            for (std::size_t j = 0; j < 3; ++j)
            {
                if (pos[j] < -1.0f || pos[j] > 1.0f)
                    _atomVelocities[i][j] = -_atomVelocities[i][j];
                pos[j] = Math::clamp(pos[j], -1.0f, 1.0f);
            }

            _atomPositions[i] = pos;
            _atomInstanceData[i].transformationMatrix.translation() = pos;
        }
        _atomInstanceBuffer.setData(_atomInstanceData, GL::BufferUsage::DynamicDraw);
    }

    void SimulationOld::updateOctree()
    {
        octreeCollisionDetection();
        _octree->update();
        arrayResize(_octreeInstanceData, 0);
        arrayAppend(_octreeInstanceData, InPlaceInit,
                    Matrix4::translation(_octree->center()) *
                        Matrix4::scaling(Vector3{_octree->halfWidth()}),
                    0x00ffff_rgbf);

        if (_drawOctreeBounds)
        {

            const auto &activeTreeNodeBlocks = _octree->activeTreeNodeBlocks();
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
                        arrayAppend(_octreeInstanceData, InPlaceInit, t, 0x197f99_rgbf);
                    }
                }
            }
        }
        _octreeInstanceBuffer.setData(_octreeInstanceData, GL::BufferUsage::DynamicDraw);
        _octreeMesh.setInstanceCount(_octreeInstanceData.size());
    }
}