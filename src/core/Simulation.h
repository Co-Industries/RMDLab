#ifndef RMD_Simulation_h
#define RMD_Simulation_h

#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/Pointer.h>
#include <Corrade/Containers/StaticArray.h>
#include <Corrade/Containers/String.h>

#include <Magnum/GL/Mesh.h>

#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Math/Constants.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/Range.h>
#include <Magnum/Math/Vector.h>
#include <Magnum/Math/Vector3.h>

#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/SceneGraph/Scene.h>

#include <Magnum/Shaders/FlatGL.h>
#include <Magnum/Shaders/PhongGL.h>

#include "../components/octree/LooseOctree.h"
#include "Data.h"

namespace Magnum
{
    using namespace Math::Literals;

    using Object3D = SceneGraph::Object<SceneGraph::MatrixTransformation3D>;
    using Scene3D = SceneGraph::Scene<SceneGraph::MatrixTransformation3D>;

    struct SimulationParameters
    {
        std::size_t atomCount;
        Float atomRadius;
        Float randomVelocity;
    };

    struct AtomInstanceData
    {
        Matrix4 transformationMatrix;
        Matrix3x3 normalMatrix;
        Color3 color;
    };

    struct OctreeInstanceData
    {
        Matrix4 transformationMatrix;
        Color3 color;
    };

    class Simulation
    {
    public:
        explicit Simulation(Scene3D &_scene, SceneGraph::DrawableGroup3D &_drawables);
        void UPDATE_ATOMS();
        void UPDATE_OCTREE();
        void RUN(const SimulationParameters &parameters);
        
        bool running = false;

    private:
        void octreeCollisionDetection();
        void checkCollisionWithSubTree(const OctreeNode &node, std::size_t i, const Vector3d &ppos, const Vector3d &pvel, const Range3D &bounds);
        void GETPARAMS();
        void INITSYSTEM();
        void ALLOCATE();

    protected:
        // ? ATOMS
        Containers::Array<AtomInstanceData> atomInstanceData;
        Containers::Array<Vector3> atomFloatPositions;
        GL::Buffer atomInstanceBuffer;
        Shaders::PhongGL atomShader;
        GL::Mesh atomMesh;

        // ? OCTREE
        Containers::Pointer<LooseOctree> octree;
        Containers::Array<OctreeInstanceData> octreeInstanceData;
        GL::Buffer octreeInstanceBuffer;
        Shaders::FlatGL3D octreeShader;
        GL::Mesh octreeMesh;
        
    };
}

#endif