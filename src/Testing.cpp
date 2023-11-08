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

#include "./octree/LooseOctree.h"
#include "./arcball/ArcBall.h"

namespace Magnum
{
  struct AtomData
  {
    bool isColliding;
    std::size_t lastCollision;
    float hue;
  };

  struct AtomInstanceData
  {
    Matrix4 transformationMatrix;
    Matrix3x3 normalMatrix;
    Color3 color;
    /* CHANGE - MISSING ATOM INFO*/
  };

  struct BoxInstanceData
  {
    Matrix4 transformationMatrix;
    Color3 color;
  };

  struct GridData
  {
    Matrix4 transformationMatrix;
    Color3 color;
  };

  struct BackgroundData
  {
    Matrix4 transformationMatrix;
  };

  class RMD : public Platform::Application
  {
  public:
    explicit RMD(const Arguments &arguments);

  protected:
    void viewportEvent(ViewportEvent &event) override;
    void keyPressEvent(KeyEvent &event) override;
    void drawEvent() override;
    void mousePressEvent(MouseEvent &event) override;
    void mouseReleaseEvent(MouseEvent &event) override;
    void mouseMoveEvent(MouseMoveEvent &event) override;
    void mouseScrollEvent(MouseScrollEvent &event) override;

    void movePoints();
    void collisionDetectionAndHandlingBruteForce();
    void collisionDetectionAndHandlingUsingOctree();
    void checkCollisionWithSubTree(const OctreeNode &node, std::size_t i,
                                   const Vector3 &ppos, const Vector3 &pvel, const Range3D &bounds);
    void drawBackground();
    void drawAtoms();
    void drawTreeNodeBoundingBoxes();

    Containers::Optional<ArcBall> _arcballCamera;
    Matrix4 _projectionMatrix;

    Containers::Array<AtomData> _atomData;
    bool _paused = false;
    bool _skipFrame = false;

    /* Points data as spheres with size */
    Containers::Array<Vector3> _atomPositions;
    Containers::Array<Vector3> _atomVelocities;
    Float _atomRadius, _atomVelocity;
    bool _animation = true;
    bool _collisionDetectionByOctree = true;

    /* RMD and boundary boxes */
    Containers::Pointer<LooseOctree> _octree;

    /* Profiling */
    DebugTools::FrameProfilerGL _profiler{
        DebugTools::FrameProfilerGL::Value::FrameTime |
            DebugTools::FrameProfilerGL::Value::CpuDuration,
        180};

    /* Atom rendering */
    GL::Mesh _sphereMesh{NoCreate};
    GL::Buffer _atomInstanceBuffer{NoCreate};
    Shaders::PhongGL _atomShader{NoCreate};
    Containers::Array<AtomInstanceData> _atomInstanceData;

    /* Treenode bounding boxes rendering */
    GL::Mesh _boxMesh{NoCreate};
    GL::Buffer _boxInstanceBuffer{NoCreate};
    Shaders::FlatGL3D _boxShader{NoCreate};
    Containers::Array<BoxInstanceData> _boxInstanceData;
    bool _drawBoundingBoxes = true;

    /* Background rendering */
    GL::Mesh _backgroundMesh{NoCreate};
    GL::Buffer _backgroundBuffer{NoCreate};
    GL::Buffer _backgroundVertexColorBuffer{NoCreate};
    Shaders::FlatGL3D _backgroundShader{NoCreate};
    Containers::Array<BackgroundData> _backgroundData;
    Trade::MeshData _backgroundMeshData = Primitives::icosphereSolid(4);
    Containers::Array<Color3> _backgroundVertexColors;
    GL::Mesh _gridMesh{NoCreate};
    GL::Buffer _gridBuffer{NoCreate};
    Shaders::FlatGL3D _gridShader{NoCreate};
    Containers::Array<GridData> _gridData;
  };

  using namespace Math::Literals;

  RMD::RMD(const Arguments &arguments) : Platform::Application{arguments, NoCreate}
  {
    Utility::Arguments args;
    args.addOption('s', "spheres", "100")
        .setHelp("spheres", "number of spheres to simulate", "N")
        .addOption('r', "sphere-radius", "0.1")
        .setHelp("sphere-radius", "sphere radius", "R")
        .addOption('v', "sphere-velocity", "1.0")
        .setHelp("sphere-velocity", "sphere velocity", "V")
        .addSkippedPrefix("magnum")
        .parse(arguments.argc, arguments.argv);

    _atomRadius = args.value<Float>("sphere-radius");
    _atomVelocity = args.value<Float>("sphere-velocity");

    /* Setup window and parameters */
    {
      const Vector2 dpiScaling = this->dpiScaling({});
      Configuration conf;
      conf.setTitle("RMD v0.1.2.1")
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
      /*
      for (Vector3 vertex : _backgroundMeshData.positions3DAsArray())
      {
        float yValue = (vertex.y() + 1.0f) / 2.0f;
        Color3 color = Color3::fromHsv(ColorHsv((360.0_degf) * yValue, 1.0f, 1.0f));
        arrayAppend(_backgroundVertexColors, InPlaceInit, color);
      }
      */
      Trade::MeshData _backgroundMutableData = MeshTools::copy(_backgroundMeshData);
      MeshTools::flipFaceWindingInPlace(_backgroundMutableData.mutableIndices());
      /* INFO Flip Normals
      MeshTools::flipNormalsInPlace(_backgroundMeshDataMutable.mutableAttribute<Vector3>(Trade::MeshAttribute::Normal));
      */

      _backgroundMesh = MeshTools::compile(_backgroundMutableData);

      _backgroundShader = Shaders::FlatGL3D{Shaders::FlatGL3D::Configuration{}.setFlags(Shaders::FlatGL3D::Flag::VertexColor | Shaders::FlatGL3D::Flag::InstancedTransformation)};
      _backgroundVertexColorBuffer = GL::Buffer{};
      for (Vector3 vertex : _backgroundMeshData.positions3DAsArray())
      {
        float yValue = (vertex.y() + 0.2f) / 8.0f + 0.3f;
        Color3 color = Color3({yValue});
        arrayAppend(_backgroundVertexColors, InPlaceInit, color);
      }
      _backgroundVertexColorBuffer.setData(_backgroundVertexColors, GL::BufferUsage::DynamicDraw);

      _backgroundMesh.addVertexBuffer(_backgroundVertexColorBuffer, 0, Shaders::FlatGL3D::Color3{});
      _backgroundBuffer = GL::Buffer{};
      _backgroundMesh.addVertexBufferInstanced(_backgroundBuffer, 1, 0, Shaders::FlatGL3D::TransformationMatrix{});

      _gridMesh = MeshTools::compile(Primitives::grid3DWireframe({16, 16}));
      _gridShader = Shaders::FlatGL3D{Shaders::FlatGL3D::Configuration{}.setFlags(Shaders::FlatGL3D::Flag::VertexColor | Shaders::FlatGL3D::Flag::InstancedTransformation)};
      _gridBuffer = GL::Buffer{};
      _gridMesh.addVertexBufferInstanced(_gridBuffer, 1, 0, Shaders::FlatGL3D::TransformationMatrix{}, Shaders::FlatGL3D::Color3{});
    }

    /* Setup atoms (render as spheres) */
    {
      const UnsignedInt numSpheres = args.value<UnsignedInt>("spheres");
      _atomPositions = Containers::Array<Vector3>{NoInit, numSpheres};
      _atomVelocities = Containers::Array<Vector3>{NoInit, numSpheres};
      _atomInstanceData = Containers::Array<AtomInstanceData>{NoInit, numSpheres};
      _atomData = Containers::Array<AtomData>{NoInit, numSpheres};

      for (std::size_t i = 0; i < numSpheres; ++i)
      {
        const Vector3 tmpPos = Vector3(std::rand(), std::rand(), std::rand()) /
                               Float(RAND_MAX);
        const Vector3 tmpVel = Vector3(std::rand(), std::rand(), std::rand()) /
                               Float(RAND_MAX);
        _atomPositions[i] = tmpPos * 2.0f - Vector3{1.0f};
        _atomPositions[i].y() *= 0.5f;
        _atomVelocities[i] = (tmpVel * 2.0f - Vector3{1.0f}).resized(_atomVelocity);

        /* Fill in the instance data. Most of this stays the same, except
           for the translation */
        _atomInstanceData[i].transformationMatrix =
            Matrix4::translation(_atomPositions[i]) *
            Matrix4::scaling(Vector3{_atomRadius});
        _atomInstanceData[i].normalMatrix =
            _atomInstanceData[i].transformationMatrix.normalMatrix();
        _atomInstanceData[i].color = Color3(1.0, 0.0, 0.0);
      }

      _atomShader = Shaders::PhongGL{Shaders::PhongGL::Configuration{}
                                         .setFlags(Shaders::PhongGL::Flag::VertexColor |
                                                   Shaders::PhongGL::Flag::InstancedTransformation)};
      _atomInstanceBuffer = GL::Buffer{};
      _sphereMesh = MeshTools::compile(Primitives::icosphereSolid(2));
      _sphereMesh.addVertexBufferInstanced(_atomInstanceBuffer, 1, 0,
                                           Shaders::PhongGL::TransformationMatrix{},
                                           Shaders::PhongGL::NormalMatrix{},
                                           Shaders::PhongGL::Color3{});
      _sphereMesh.setInstanceCount(_atomInstanceData.size());
    }

    /* Setup octree */
    {
      /* RMD nodes should have half width no smaller than the sphere
         radius */
      _octree.emplace(Vector3{0}, 1.0f, Math::max(_atomRadius, 0.1f));

      _octree->setPoints(_atomPositions);
      _octree->build();
      Debug{} << "  Allocated nodes:" << _octree->numAllocatedNodes();
      Debug{} << "  Max number of points per node:" << _octree->maxNumPointInNodes();

      /* Disable profiler by default */
      _profiler.disable();
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

    Debug{} << "Collision detection using octree";
  }

  void RMD::drawEvent()
  {
    GL::defaultFramebuffer.clear(GL::FramebufferClear::Color | GL::FramebufferClear::Depth);

    if (!_paused || _skipFrame)
    {
      _skipFrame = false;
      collisionDetectionAndHandlingUsingOctree();
      movePoints();
      _octree->update();
    }
    /* Update camera before drawing instances */
    const bool moving = _arcballCamera->updateTransformation();

    drawBackground();
    drawAtoms();
    drawTreeNodeBoundingBoxes();

    swapBuffers();

    /* If the camera is moving or the animation is running, redraw immediately */
    if (moving || _animation)
      redraw();
  }

  void RMD::collisionDetectionAndHandlingBruteForce()
  {
    for (std::size_t i = 0; i < _atomPositions.size(); ++i)
    {
      const Vector3 ppos = _atomPositions[i];
      const Vector3 pvel = _atomVelocities[i];
      for (std::size_t j = i + 1; j < _atomPositions.size(); ++j)
      {
        const Vector3 qpos = _atomPositions[j];
        const Vector3 qvel = _atomVelocities[j];
        const Vector3 velpq = pvel - qvel;
        const Vector3 pospq = ppos - qpos;
        const Float vp = Math::dot(velpq, pospq);
        if (vp < 0.0f)
        {
          const Float dpq = pospq.length();
          if (dpq < 2.0f * _atomRadius)
          {
            const Vector3 vNormal = vp * pospq / (dpq * dpq);
            _atomVelocities[i] = (_atomVelocities[i] - vNormal).resized(_atomVelocity);
            _atomVelocities[j] = (_atomVelocities[j] + vNormal).resized(_atomVelocity);
          }
        }
      }
    }
  }

  void RMD::collisionDetectionAndHandlingUsingOctree()
  {
    const OctreeNode &rootNode = _octree->rootNode();
    for (std::size_t i = 0; i < _atomPositions.size(); ++i)
    {
      checkCollisionWithSubTree(rootNode, i,
                                _atomPositions[i], _atomVelocities[i],
                                Range3D::fromCenter(_atomPositions[i], Vector3{_atomRadius}));
    }
  }

  void RMD::checkCollisionWithSubTree(const OctreeNode &node,
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
        const Vector3 pospq = ppos - qpos;
        const Vector3 qvel = _atomVelocities[j];
        const Vector3 velpq = pvel - qvel;
        const Float vp = Math::dot(velpq, pospq);
        /*
        if (vp < 0.0f)
        {
          */
        const Float dpq = pospq.length();
        if (dpq < 2.0f * _atomRadius)
        {
          const Vector3 vNormal = vp * pospq / (dpq * dpq);
          _atomVelocities[i] = (_atomVelocities[i] - vNormal).resized(_atomVelocity);
          _atomVelocities[j] = (_atomVelocities[j] + vNormal).resized(_atomVelocity);
          _atomData[i].isColliding = true;
          _atomData[j].isColliding = true;
          _atomData[i].lastCollision = j;
          _atomData[j].lastCollision = i;
          _atomInstanceData[i].color = Color3::fromHsv(ColorHsv(_atomData[i].hue * 360.0_degf, 1.0, 2.0));
          _atomInstanceData[j].color = Color3::fromHsv(ColorHsv(_atomData[i].hue * 360.0_degf, 1.0, 2.0));
          _atomData[i].hue += 0.0001f;
          _atomData[j].hue += 0.0001f;
          continue;
        }

        if (_atomData[i].lastCollision == j)
        {
          _atomData[i].isColliding = false;
          _atomInstanceData[i].color = Color3(1.0, 0.0, 0.0);
        }
        if (_atomData[j].lastCollision == i)
        {
          _atomData[j].isColliding = false;
          _atomInstanceData[j].color = Color3(1.0, 0.0, 0.0);
        }
        if (_atomData[i].isColliding)
        {
          _atomInstanceData[i].color = Color3::fromHsv(ColorHsv(_atomData[i].hue * 360.0_degf, 1.0, 2.0));
          _atomData[i].hue += 0.0001f;
        }
        if (_atomData[j].isColliding)
        {
          _atomInstanceData[j].color = Color3::fromHsv(ColorHsv(_atomData[j].hue * 360.0_degf, 1.0, 2.0));
          _atomData[j].hue += 0.0001f;
        }
      }
    }
  }

  void RMD::movePoints()
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
    }
  }

  auto _backgroundHueOffset = 0.0_degf;
  double _backgroundLightnessOffset = -1.0f;
  bool _backgroundDirection = true;
  void RMD::drawBackground()
  {
    arrayResize(_backgroundData, 0);
    arrayAppend(_backgroundData, InPlaceInit,
                _arcballCamera->viewMatrix() *
                    Matrix4::translation({_octree->center()}) *
                    Matrix4::scaling(Vector3{20.f}));
    _backgroundBuffer.setData(_backgroundData, GL::BufferUsage::DynamicDraw);
    _backgroundMesh.setInstanceCount(_backgroundData.size());
    _backgroundShader
        .setTransformationProjectionMatrix(_projectionMatrix)
        .draw(_backgroundMesh);

    arrayResize(_gridData, 0);
    arrayAppend(_gridData, InPlaceInit,
                _arcballCamera->viewMatrix() *
                    Matrix4::translation(Vector3(0.0f, -_octree->halfWidth(), 0.0f)) *
                    Matrix4::translation({_octree->center()}) *
                    Matrix4::scaling(Vector3{2 * _octree->halfWidth()}) *
                    Matrix4::rotationX(90.0_degf),
                Color3(0.7f, 0.7f, 0.7f));
    _gridBuffer.setData(_gridData, GL::BufferUsage::DynamicDraw);
    _gridMesh.setInstanceCount(_gridData.size());
    _gridShader
        .setTransformationProjectionMatrix(_projectionMatrix)
        .draw(_gridMesh);
  }

  void RMD::drawAtoms()
  {
    for (std::size_t i = 0; i != _atomPositions.size(); ++i)
    {
      _atomInstanceData[i].transformationMatrix.translation() =
          _atomPositions[i];
    }

    _atomInstanceBuffer.setData(_atomInstanceData, GL::BufferUsage::DynamicDraw);
    _atomShader
        .setProjectionMatrix(_projectionMatrix)
        .setTransformationMatrix(_arcballCamera->viewMatrix())
        .setNormalMatrix(_arcballCamera->viewMatrix().normalMatrix())
        .draw(_sphereMesh);
  }

  void RMD::drawTreeNodeBoundingBoxes()
  {
    arrayResize(_boxInstanceData, 0);

    /* Always draw the root node */
    arrayAppend(_boxInstanceData, InPlaceInit,
                _arcballCamera->viewMatrix() *
                    Matrix4::translation(_octree->center()) *
                    Matrix4::scaling(Vector3{_octree->halfWidth()}),
                0x00ffff_rgbf);

    /* Draw the remaining non-empty nodes */
    if (_drawBoundingBoxes)
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
            const Matrix4 t = _arcballCamera->viewMatrix() *
                              Matrix4::translation(pNode.center()) *
                              Matrix4::scaling(Vector3{pNode.halfWidth()});
            arrayAppend(_boxInstanceData, InPlaceInit, t, 0x197f99_rgbf);
          }
        }
      }
    }

    _boxInstanceBuffer.setData(_boxInstanceData, GL::BufferUsage::DynamicDraw);
    _boxMesh.setInstanceCount(_boxInstanceData.size());
    _boxShader.setTransformationProjectionMatrix(_projectionMatrix)
        .draw(_boxMesh);
  }

  void RMD::viewportEvent(ViewportEvent &event)
  {
    GL::defaultFramebuffer.setViewport({{}, event.framebufferSize()});
    _arcballCamera->reshape(event.windowSize());

    _projectionMatrix = Matrix4::perspectiveProjection(_arcballCamera->fov(),
                                                       Vector2{event.framebufferSize()}.aspectRatio(), 0.01f, 100.0f);
  }

  void RMD::keyPressEvent(KeyEvent &event)
  {
    if (event.key() == KeyEvent::Key::B)
    {
      _drawBoundingBoxes ^= true;
    }
    else if (event.key() == KeyEvent::Key::O)
    {
      if ((_collisionDetectionByOctree ^= true))
        Debug{} << "Collision detection using octree";
      else
        Debug{} << "Collision detection using brute force";
      /* Reset the profiler to avoid measurements of the two methods mixed
         together */
      if (_profiler.isEnabled())
        _profiler.enable();
    }
    else if (event.key() == KeyEvent::Key::P)
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
    else if (event.key() == KeyEvent::Key::Space)
    {
      _paused ^= true;
    }
    else if (event.key() == KeyEvent::Key::Right)
    {
      _skipFrame = true;
    }
    else
      return;

    event.setAccepted();
    redraw();
  }

  void RMD::mousePressEvent(MouseEvent &event)
  {
    /* Enable mouse capture so the mouse can drag outside of the window */
    /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
    SDL_CaptureMouse(SDL_TRUE);
    _arcballCamera->initTransformation(event.position());
    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
  }

  void RMD::mouseReleaseEvent(MouseEvent &)
  {
    /* Disable mouse capture again */
    /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
    SDL_CaptureMouse(SDL_FALSE);
  }

  void RMD::mouseMoveEvent(MouseMoveEvent &event)
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

  void RMD::mouseScrollEvent(MouseScrollEvent &event)
  {
    const Float delta = event.offset().y();
    if (Math::abs(delta) < 1.0e-2f)
      return;

    _arcballCamera->zoom(delta);

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
  }
}

MAGNUM_APPLICATION_MAIN(Magnum::RMD)