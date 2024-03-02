#ifndef RMD_Grid_h
#define RMD_Grid_h

#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Shaders/FlatGL.h>

namespace Magnum
{
    using Object3D = SceneGraph::Object<SceneGraph::MatrixTransformation3D>;
    using Scene3D = SceneGraph::Scene<SceneGraph::MatrixTransformation3D>;

    class Grid
    {
    public:
        explicit Grid(
            Scene3D &scene,
            SceneGraph::DrawableGroup3D &drawables,
            float size,
            Vector2i subdivisions,
            Color3 color);

    protected:
        Scene3D &_scene;
        SceneGraph::DrawableGroup3D &_drawables;
        GL::Mesh _mesh;
        Shaders::FlatGL3D _shader;
    };
}

#endif