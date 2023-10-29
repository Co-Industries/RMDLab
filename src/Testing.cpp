/*
    This file is part of Magnum.

    Original authors — credit is appreciated but not required:

        2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019,
        2020, 2021, 2022, 2023 — Vladimír Vondruš <mosra@centrum.cz>

    This is free and unencumbered software released into the public domain.

    Anyone is free to copy, modify, publish, use, compile, sell, or distribute
    this software, either in source code form or as a compiled binary, for any
    purpose, commercial or non-commercial, and by any means.

    In jurisdictions that recognize copyright laws, the author or authors of
    this software dedicate any and all copyright interest in the software to
    the public domain. We make this dedication for the benefit of the public
    at large and to the detriment of our heirs and successors. We intend this
    dedication to be an overt act of relinquishment in perpetuity of all
    present and future rights to this software under copyright law.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
    IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/Renderer.h>

#include <Magnum/Math/Color.h>
#include <Magnum/Math/Matrix4.h>

#include <Magnum/Platform/Sdl2Application.h>
#include <Corrade/Containers/GrowableArray.h>

#include <Magnum/Primitives/Icosphere.h>
#include <Magnum/Primitives/Cube.h>
#include <Magnum/Primitives/Grid.h>

#include <Magnum/Shaders/FlatGL.h>
#include <Magnum/Shaders/VertexColorGL.h>

#include <Magnum/Trade/MeshData.h>
#include <Magnum/Trade/Data.h>

#include <Magnum/MeshTools/Compile.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/MeshTools/FlipNormals.h>
#include <Magnum/MeshTools/Copy.h>
#include <Magnum/MeshTools/Transform.h>

namespace Magnum
{

    using namespace Magnum::Math::Literals;

    class Testing : public Platform::Application
    {
    public:
        explicit Testing(const Arguments &arguments);

    private:
        void drawEvent() override;
        void mousePressEvent(MouseEvent &event) override;
        void mouseReleaseEvent(MouseEvent &event) override;
        void mouseMoveEvent(MouseMoveEvent &event) override;

        GL::Mesh _backgroundMesh;
        GL::Mesh _gridMesh;
        Shaders::VertexColorGL3D _backgroundShader;
        Shaders::FlatGL3D _gridShader;
        Containers::Array<Color3> _vertexColors;

        Matrix4 _transformation, _projection;
        Color3 _color;
    };

    Testing::Testing(const Arguments &arguments) : Platform::Application{arguments, Configuration{}
                                                                                        .setTitle("Testing")}
    {
        GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
        GL::Renderer::enable(GL::Renderer::Feature::DepthTest);

        Trade::MeshData _backgroundMeshData = Primitives::icosphereSolid(2);

        for (Vector3 _vertex : _backgroundMeshData.positions3DAsArray())
        {
            float yValue = ((_vertex.y() + 1.0f) / 8.0f) + 0.1f;
            Color3 _vertexColor = Color3({yValue});
            arrayAppend(_vertexColors, InPlaceInit, _vertexColor);
        }

        Trade::MeshData _backgroundMeshDataMutable = MeshTools::copy(_backgroundMeshData);
        MeshTools::flipFaceWindingInPlace(_backgroundMeshDataMutable.mutableIndices());
        MeshTools::transformPointsInPlace(Matrix4::scaling(Vector3{20.0f}),
                                          _backgroundMeshDataMutable.mutableAttribute<Vector3>(Trade::MeshAttribute::Position));
        /* INFO Flip Normals
        MeshTools::flipNormalsInPlace(_backgroundMeshDataMutable.mutableAttribute<Vector3>(Trade::MeshAttribute::Normal));
        */

        _backgroundMesh = MeshTools::compile(_backgroundMeshDataMutable);
        _backgroundMesh.addVertexBuffer(GL::Buffer{_vertexColors}, 0, Shaders::VertexColorGL3D::Color3{});

        _transformation = Matrix4::rotationX(0.0_degf) * Matrix4::rotationY(0.0_degf);

        _projection = Matrix4::perspectiveProjection(35.0_degf, Vector2{windowSize()}.aspectRatio(), 0.01f, 100.0f) *
                      Matrix4::translation(Vector3::zAxis(-10.0f));

        Trade::MeshData _gridMeshData = Primitives::grid3DWireframe({32, 32});
        Trade::MeshData _gridMeshDataMutable = MeshTools::copy(_gridMeshData);
        MeshTools::transformPointsInPlace(Matrix4::scaling(Vector3{32.0f}) * Matrix4::rotationX(90.0_degf),
                                          _gridMeshDataMutable.mutableAttribute<Vector3>(Trade::MeshAttribute::Position));
        _gridMesh = MeshTools::compile(_gridMeshDataMutable);
    }

    void Testing::drawEvent()
    {
        GL::defaultFramebuffer.clear(
            GL::FramebufferClear::Color | GL::FramebufferClear::Depth);
        _backgroundShader.setTransformationProjectionMatrix(_projection * _transformation)
            .draw(_backgroundMesh);
        _gridShader.setTransformationProjectionMatrix(_projection * _transformation)
            .draw(_gridMesh);
        swapBuffers();
    }

    void Testing::mousePressEvent(MouseEvent &event)
    {
        if (event.button() != MouseEvent::Button::Left)
            return;

        event.setAccepted();
    }

    void Testing::mouseReleaseEvent(MouseEvent &event)
    {
        _color = Color3::fromHsv({_color.hue() + 50.0_degf, 1.0f, 1.0f});

        event.setAccepted();
        redraw();
    }

    void Testing::mouseMoveEvent(MouseMoveEvent &event)
    {
        if (!(event.buttons() & MouseMoveEvent::Button::Left))
            return;

        Vector2 delta = 3.0f * Vector2{event.relativePosition()} / Vector2{windowSize()};

        _transformation =
            Matrix4::rotationX(Rad{delta.y()}) *
            _transformation *
            Matrix4::rotationY(Rad{delta.x()});

        event.setAccepted();
        redraw();
    }
}

MAGNUM_APPLICATION_MAIN(Magnum::Testing)