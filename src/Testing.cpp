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

/*
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

    Testing::Testing(const Arguments &arguments) : Platform::Application{arguments, NoCreate}
    {
        const Vector2 dpiScaling = this->dpiScaling({});
        Configuration conf;
        conf.setTitle("Testing")
            .setSize(conf.size(), dpiScaling)
            .setWindowFlags(Configuration::WindowFlag::Resizable);
        GLConfiguration glConf;
        glConf.setSampleCount(dpiScaling.max() < 2.0f ? 8 : 2);
        if (!tryCreate(conf, glConf))
        {
            create(conf, glConf.setSampleCount(0));
        }
        GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
        GL::Renderer::enable(GL::Renderer::Feature::DepthTest);

        Trade::MeshData _backgroundMeshData = Primitives::icosphereSolid(2);

        for (Vector3 _vertex : _backgroundMeshData.positions3DAsArray())
        {
        float yValue = ((_vertex.y() + 1.0f) / 6.0f) + 0.0f;
        Color3 _vertexColor = Color3({yValue});
        arrayAppend(_vertexColors, InPlaceInit, _vertexColor);
        }

        Trade::MeshData _backgroundMeshDataMutable = MeshTools::copy(_backgroundMeshData);
        MeshTools::flipFaceWindingInPlace(_backgroundMeshDataMutable.mutableIndices());
        MeshTools::transformPointsInPlace(Matrix4::scaling(Vector3{100.0f}),
                                          _backgroundMeshDataMutable.mutableAttribute<Vector3>(Trade::MeshAttribute::Position));
        INFO Flip Normals
        MeshTools::flipNormalsInPlace(_backgroundMeshDataMutable.mutableAttribute<Vector3>(Trade::MeshAttribute::Normal));


        _backgroundMesh = MeshTools::compile(_backgroundMeshDataMutable);
        _backgroundMesh.addVertexBuffer(GL::Buffer{_vertexColors}, 0, Shaders::VertexColorGL3D::Color3{});

        _transformation = Matrix4::rotationX(0.0_degf) * Matrix4::rotationY(0.0_degf);

        _projection = Matrix4::perspectiveProjection(35.0_degf, Vector2{windowSize()}.aspectRatio(), 0.01f, 1000.0f) *
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

*/

#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/Optional.h>
#include <Corrade/Containers/Pointer.h>
#include <Corrade/Utility/Arguments.h>

#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/DebugTools/FrameProfiler.h>

#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/Renderer.h>

#include <Magnum/Math/Color.h>
#include <Magnum/Math/Matrix4.h>

#include <Magnum/Primitives/Icosphere.h>
#include <Magnum/Primitives/Cube.h>
#include <Magnum/Primitives/Grid.h>

#include <Magnum/Shaders/FlatGL.h>
#include <Magnum/Shaders/VertexColorGL.h>
#include <Magnum/Shaders/PhongGL.h>

#include <Magnum/Trade/MeshData.h>
#include <Magnum/Trade/Data.h>

#include <Magnum/MeshTools/Compile.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/MeshTools/FlipNormals.h>
#include <Magnum/MeshTools/Copy.h>
#include <Magnum/MeshTools/Transform.h>

#include "./arcball/ArcBall.h"

namespace Magnum
{
  struct BoxInstanceData
  {
    Matrix4 transformationMatrix;
    Color3 color;
  };

  struct GridData
  {
    Matrix4 transformationMatrix;
  };

  class Testing : public Platform::Application
  {
  public:
    explicit Testing(const Arguments &arguments);

  protected:
    void viewportEvent(ViewportEvent &event) override;
    void keyPressEvent(KeyEvent &event) override;
    void drawEvent() override;
    void mousePressEvent(MouseEvent &event) override;
    void mouseReleaseEvent(MouseEvent &event) override;
    void mouseMoveEvent(MouseMoveEvent &event) override;
    void mouseScrollEvent(MouseScrollEvent &event) override;

    void drawBackground();
    void drawBox();

    Containers::Optional<ArcBall> _arcballCamera;
    Matrix4 _projectionMatrix;

    /* Profiling */
    DebugTools::FrameProfilerGL _profiler{
        DebugTools::FrameProfilerGL::Value::FrameTime |
            DebugTools::FrameProfilerGL::Value::CpuDuration,
        180};

    GL::Mesh _boxMesh{NoCreate};
    GL::Buffer _boxInstanceBuffer{NoCreate};
    Shaders::FlatGL3D _boxShader{NoCreate};
    Containers::Array<BoxInstanceData> _boxInstanceData;

    /* Background rendering */

    GL::Mesh _gridMesh{NoCreate};
    GL::Buffer _gridBuffer{NoCreate};
    Shaders::FlatGL3D _gridShader{NoCreate};
    Containers::Array<GridData> _gridData;

    GL::Mesh _backgroundMesh{NoCreate};
    GL::Buffer _backgroundBuffer{NoCreate};
    Shaders::VertexColorGL3D _backgroundShader{NoCreate};
  };

  using namespace Math::Literals;

  Testing::Testing(const Arguments &arguments) : Platform::Application{arguments, NoCreate}
  {
    /* Setup window and parameters */
    {
      const Vector2 dpiScaling = this->dpiScaling({});
      Configuration conf;
      conf.setTitle("Testing")
          .setSize(conf.size(), dpiScaling)
          .setWindowFlags(Configuration::WindowFlag::Resizable);
      GLConfiguration glConf;
      glConf.setSampleCount(dpiScaling.max() < 2.0f ? 8 : 2);
      if (!tryCreate(conf, glConf))
      {
        create(conf, glConf.setSampleCount(0));
      }

      GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
      GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);

      /* Loop at 60 Hz max */
      setSwapInterval(1);
      /* 16 */
      setMinimalLoopPeriod(16);
    }

    /* Setup camera */
    {
      const Vector3 eye = Vector3::zAxis(5.0f);
      const Vector3 viewCenter;
      const Vector3 up = Vector3::yAxis();
      const Deg fov = 45.0_degf;
      _arcballCamera.emplace(eye, viewCenter, up, fov, windowSize());
      _arcballCamera->setLagging(0.85f);

      _projectionMatrix = Matrix4::perspectiveProjection(fov,
                                                         Vector2{framebufferSize()}.aspectRatio(), 0.01f, 100.0f);
    }

    /* Setup background */
    {
      Trade::MeshData _backgroundMeshData = Primitives::icosphereSolid(2);
      Containers::Array<Color3> _vertexColors;

      for (Vector3 vertex : _backgroundMeshData.positions3DAsArray())
      {
        float yValue = ((vertex.y() + 1.0f) / 8.0f) + 0.1f;
        Color3 color = Color3({yValue});
        arrayAppend(_vertexColors, InPlaceInit, color);
      }

      Trade::MeshData _backgroundMutableData = MeshTools::copy(_backgroundMeshData);
      MeshTools::flipFaceWindingInPlace(_backgroundMutableData.mutableIndices());
      MeshTools::transformPointsInPlace(Matrix4::scaling(Vector3{20.0f}) * Matrix4::translation(Vector3{0.0f}),
                                        _backgroundMutableData.mutableAttribute<Vector3>(Trade::MeshAttribute::Position));
      /* INFO Flip Normals
      MeshTools::flipNormalsInPlace(_backgroundMeshDataMutable.mutableAttribute<Vector3>(Trade::MeshAttribute::Normal));
      */

      _backgroundMesh = MeshTools::compile(_backgroundMutableData);
      /* TODO -  Buffer with Color3 and Position*/
      _backgroundMesh.addVertexBuffer(GL::Buffer{_vertexColors}, 0, Shaders::VertexColorGL3D::Color3{});
    }

    /* Setup grid */
    {
      _gridBuffer = GL::Buffer{};
      /*
      Trade::MeshData _gridMeshData = Primitives::grid3DWireframe({32, 32});

      Trade::MeshData _gridMutableData = MeshTools::copy(_gridMeshData);
      MeshTools::transformPointsInPlace(Matrix4::scaling(Vector3{32.0f}) * Matrix4::rotationX(85.0_degf),
                                        _gridMutableData.mutableAttribute<Vector3>(Trade::MeshAttribute::Position));
      */
      _gridMesh = MeshTools::compile(Primitives::grid3DWireframe({16, 16}));
      _gridShader = Shaders::FlatGL3D{Shaders::FlatGL3D::Configuration{}.setFlags(Shaders::FlatGL3D::Flag::InstancedTransformation)};
      _gridMesh.addVertexBufferInstanced(_gridBuffer, 1, 0, Shaders::FlatGL3D::TransformationMatrix{});
    }

    /* Treenode bounding boxes render variables */
    {
      _boxShader = Shaders::FlatGL3D{Shaders::FlatGL3D::Configuration{}
                                         .setFlags(Shaders::FlatGL3D::Flag::VertexColor |
                                                   Shaders::FlatGL3D::Flag::InstancedTransformation)};
      _boxInstanceBuffer = GL::Buffer{};
      _boxMesh = MeshTools::compile(Primitives::cubeWireframe());
      _boxMesh.addVertexBufferInstanced(_boxInstanceBuffer, 1, 0,
                                        Shaders::FlatGL3D::TransformationMatrix{},
                                        Shaders::FlatGL3D::Color3{});
    }
  }

  void Testing::drawEvent()
  {
    GL::defaultFramebuffer.clear(GL::FramebufferClear::Color | GL::FramebufferClear::Depth);
    /* Update camera before drawing instances */
    const bool moving = _arcballCamera->updateTransformation();

    drawBackground();
    drawBox();
    swapBuffers();

    /* If the camera is moving or the animation is running, redraw immediately */
    if (moving)
      redraw();
  }

  void Testing::drawBackground()
  {
    _backgroundShader
        .setTransformationProjectionMatrix(_arcballCamera->viewMatrix())
        .draw(_backgroundMesh);

    arrayResize(_gridData, 0);
    arrayAppend(_gridData, InPlaceInit,
                _arcballCamera->viewMatrix() *
                    Matrix4::translation(Vector3{0.0f}) *
                    Matrix4::scaling(Vector3{10.f}) *
                    Matrix4::rotationX(90.0_degf));
    _gridBuffer.setData(_gridData, GL::BufferUsage::DynamicDraw);
    _gridMesh.setInstanceCount(_gridData.size());
    _gridShader
        .setTransformationProjectionMatrix(_projectionMatrix)
        .draw(_gridMesh);
  }

  void Testing::drawBox()
  {
    arrayResize(_boxInstanceData, 0);

    /* Always draw the root node */
    arrayAppend(_boxInstanceData, InPlaceInit,
                _arcballCamera->viewMatrix() *
                    Matrix4::translation(Vector3{0.0f}) *
                    Matrix4::scaling(Vector3{5.0f}),
                0x00ffff_rgbf);
    _boxInstanceBuffer.setData(_boxInstanceData, GL::BufferUsage::DynamicDraw);
    _boxMesh.setInstanceCount(_boxInstanceData.size());
    _boxShader.setTransformationProjectionMatrix(_projectionMatrix)
        .draw(_boxMesh);
  }

  void Testing::viewportEvent(ViewportEvent &event)
  {
    GL::defaultFramebuffer.setViewport({{}, event.framebufferSize()});
    _arcballCamera->reshape(event.windowSize());

    _projectionMatrix = Matrix4::perspectiveProjection(_arcballCamera->fov(),
                                                       Vector2{event.framebufferSize()}.aspectRatio(), 0.01f, 100.0f);
  }

  void Testing::keyPressEvent(KeyEvent &event)
  {
    if (event.key() == KeyEvent::Key::P)
    {
      if (_profiler.isEnabled())
        _profiler.disable();
      else
        _profiler.enable();
    }
    else if (event.key() == KeyEvent::Key::R)
    {
      _arcballCamera->reset();
    }
    else
      return;

    event.setAccepted();
    redraw();
  }

  void Testing::mousePressEvent(MouseEvent &event)
  {
    /* Enable mouse capture so the mouse can drag outside of the window */
    /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
    SDL_CaptureMouse(SDL_TRUE);
    _arcballCamera->initTransformation(event.position());
    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
  }

  void Testing::mouseReleaseEvent(MouseEvent &)
  {
    /* Disable mouse capture again */
    /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
    SDL_CaptureMouse(SDL_FALSE);
  }

  void Testing::mouseMoveEvent(MouseMoveEvent &event)
  {
    if (!event.buttons())
      return;

    if (event.modifiers() & MouseMoveEvent::Modifier::Shift)
      _arcballCamera->translate(event.position());
    else
      _arcballCamera->rotate(event.position());

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
  }

  void Testing::mouseScrollEvent(MouseScrollEvent &event)
  {
    const Float delta = event.offset().y();
    if (Math::abs(delta) < 1.0e-2f)
      return;

    _arcballCamera->zoom(delta);

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
  }
}

MAGNUM_APPLICATION_MAIN(Magnum::Testing)