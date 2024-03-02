#include <Magnum/Trade/MeshData.h>

#include <Magnum/MeshTools/Compile.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/MeshTools/FlipNormals.h>
#include <Magnum/MeshTools/Copy.h>
#include <Magnum/MeshTools/Transform.h>

#include <Magnum/Shaders/FlatGL.h>
#include <Magnum/Primitives/Icosphere.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Corrade/Containers/GrowableArray.h>

#include <Magnum/Math/Color.h>
#include <Magnum/Math/Vector.h>
#include <Magnum/Math/Matrix4.h>

#include "Skybox.h"
#include "../scene/Scene.h"

namespace Magnum
{
        using namespace Math::Literals;

        Skybox::Skybox(
            Scene3D &scene,
            SceneGraph::DrawableGroup3D &drawables,
            float size) : _scene(scene), _drawables(drawables)
        {
                const Trade::MeshData _meshData = Primitives::icosphereSolid(4);
                Trade::MeshData _mutableMeshData = MeshTools::copy(_meshData);
                MeshTools::flipFaceWindingInPlace(_mutableMeshData.mutableIndices());
                MeshTools::transformPointsInPlace(Matrix4::scaling(Vector3{size}),
                                                  _mutableMeshData.mutableAttribute<Vector3>(Trade::MeshAttribute::Position));
                _mesh = MeshTools::compile(_mutableMeshData);

                _shader = Shaders::FlatGL3D{Shaders::FlatGL3D::Configuration{}.setFlags(Shaders::FlatGL3D::Flag::VertexColor)};
                Containers::Array<Color3> colorData;
                GL::Buffer vertices;
                for (Vector3 vertex : _meshData.positions3DAsArray())
                {
                        float yValue = (vertex.y() + 0.2f) / 8.0f + 0.3f;
                        arrayAppend(colorData, InPlaceInit, Color3({yValue}));
                }

                vertices.setData(MeshTools::interleave(_mutableMeshData.positions3DAsArray(), colorData));
                _mesh.addVertexBuffer(std::move(vertices), 0, Shaders::FlatGL3D::Position{}, Shaders::FlatGL3D::Color3{});
                auto object = new Object3D{&_scene};
                new FlatGLDrawable{*object, _shader, _mesh, _drawables};
        }
}