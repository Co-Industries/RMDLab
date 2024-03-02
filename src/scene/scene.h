#ifndef RMD_Scene_h
#define RMD_Scene_h

#include <Magnum/Math/Matrix4.h>

#include <Magnum/Shaders/FlatGL.h>

#include <Magnum/SceneGraph/Object.h>
#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/Camera.h>

namespace Magnum
{
    using Object3D = SceneGraph::Object<SceneGraph::MatrixTransformation3D>;

    class FlatGLDrawable : public SceneGraph::Drawable3D
    {
    public:
        explicit FlatGLDrawable(Object3D &object, Shaders::FlatGL3D &shader, GL::Mesh &mesh, SceneGraph::DrawableGroup3D &drawables);

        void draw(const Matrix4 &transformation, SceneGraph::Camera3D &camera);

    protected:
        Shaders::FlatGL3D &_shader;
        GL::Mesh &_mesh;
    };
}
#endif