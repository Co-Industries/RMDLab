#include <Magnum/Math/Matrix4.h>

#include <Magnum/MeshTools/Compile.h>
#include <Magnum/Meshtools/Copy.h>
#include <Magnum/MeshTools/Transform.h>

#include <Magnum/Primitives/Grid.h>
#include <Magnum/Trade/MeshData.h>

#include "Grid.h"
#include "../components/scene/Scene.h"

namespace Magnum
{
    using namespace Math::Literals;

    Grid::Grid(
        Scene3D &scene,
        SceneGraph::DrawableGroup3D &drawables,
        float size,
        Vector2i subdivisions,
        Color3 color) : _scene(scene), _drawables(drawables)
    {
        const Trade::MeshData _meshData = Primitives::grid3DWireframe(subdivisions);
        Trade::MeshData _mutableMeshData = MeshTools::copy(_meshData);
        MeshTools::transformPointsInPlace(Matrix4::translation(Vector3(0.0f, -150.0f, 0.0f)) *
                                              Matrix4::scaling(Vector3{size}) *
                                              Matrix4::rotationX(90.0_degf),
                                          _mutableMeshData.mutableAttribute<Vector3>(Trade::MeshAttribute::Position));
        _mesh = MeshTools::compile(_mutableMeshData);
        _shader = Shaders::FlatGL3D{};
        _shader.setColor(color);
        auto object = new Object3D{&_scene};
        new FlatGLDrawable{*object, _shader, _mesh, _drawables};
    }
}