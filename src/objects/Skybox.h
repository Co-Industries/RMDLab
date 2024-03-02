#ifndef RMD_Skybox_h
#define RMD_Skybox_h

#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Drawable.h>

namespace Magnum
{
    using Scene3D = SceneGraph::Scene<SceneGraph::MatrixTransformation3D>;
    class Skybox
    {
    public:
        explicit Skybox(Scene3D &scene, SceneGraph::DrawableGroup3D &drawables);

    protected:
        Scene3D &_scene;
        SceneGraph::DrawableGroup3D &_drawables;
        GL::Mesh _mesh;
        Shaders::FlatGL3D _shader;
    };
}

#endif