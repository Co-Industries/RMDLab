#ifndef RMD_SimulationOld_h
#define RMD_SimulationOld_h

/* A B C D E F G H I J K L M N O P Q R S T U V W X Y Z */
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/Pointer.h>
#include <unordered_set>

#include <Magnum/GL/Mesh.h>

#include <Magnum/Math/Color.h>
#include <Magnum/Math/Matrix4.h>

#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/SceneGraph/Scene.h>

#include <Magnum/Shaders/FlatGL.h>
#include <Magnum/Shaders/PhongGL.h>

#include "../components/octree/LooseOctree.h"

namespace Magnum
{
    struct AtomInstanceDataOld
    {
        Matrix4 transformationMatrix;
        Matrix3x3 normalMatrix;
        Color3 color;
    };

    // struct AtomData
    //{
    //     Containers::Array<int> collisions;
    // };

    struct OctreeInstanceDataOld
    {
        Matrix4 transformationMatrix;
        Color3 color;
    };

    using Object3D = SceneGraph::Object<SceneGraph::MatrixTransformation3D>;
    using Scene3D = SceneGraph::Scene<SceneGraph::MatrixTransformation3D>;

    using namespace Math::Literals;

    class SimulationOld
    {
    public:
        explicit SimulationOld(Scene3D &scene, SceneGraph::DrawableGroup3D &drawables, UnsignedInt atomCount, bool &drawOctreeBounds);
        void updateAtoms();
        void updateOctree();
        void updateColor(Color3 color);

    private:
        void octreeCollisionDetection();
        bool findCollision(const Containers::Array<Int> &collisions, Int value);
        void removeCollision(Containers::Array<Int> &collisions, Int value);
        void checkCollisionWithSubTree(const OctreeNode &node, std::size_t i,
                                       const Vector3 &ppos, const Vector3 &pvel, const Range3D &bounds);

    protected:
        Scene3D &_scene;
        SceneGraph::DrawableGroup3D &_drawables;

        Containers::Array<AtomInstanceDataOld> _atomInstanceData;
        // Containers::Array<AtomData> _atomData;
        Containers::Array<Vector3> _atomPositions;
        Containers::Array<Vector3> _atomVelocities;
        GL::Buffer _atomInstanceBuffer;
        Float _atomRadius, _atomVelocity, _atomRange;
        UnsignedInt _atomCount;

        Containers::Pointer<LooseOctree> _octree;
        Containers::Array<OctreeInstanceDataOld> _octreeInstanceData;
        GL::Buffer _octreeInstanceBuffer;
        bool &_drawOctreeBounds;

        Shaders::PhongGL _atomShader;
        Shaders::FlatGL3D _octreeShader;
        GL::Mesh _atomMesh;
        GL::Mesh _octreeMesh;
    };
}

#endif