#include "Scene.h"

namespace Magnum
{
    using Object3D = SceneGraph::Object<SceneGraph::MatrixTransformation3D>;

    FlatGLDrawable::FlatGLDrawable(Object3D &object,
                                   Shaders::FlatGL3D &shader,
                                   GL::Mesh &mesh,
                                   SceneGraph::DrawableGroup3D &drawables) : SceneGraph::Drawable3D{object, &drawables},
                                                                             _shader(shader),
                                                                             _mesh(mesh) {}

    void FlatGLDrawable::draw(const Matrix4 &transformation, SceneGraph::Camera3D &camera)
    {
        _shader
            .setTransformationProjectionMatrix(camera.projectionMatrix() * transformation)
            .draw(_mesh);
    };

    AtomDrawable::AtomDrawable(Object3D &object,
                               Shaders::PhongGL &shader,
                               GL::Mesh &mesh,
                               SceneGraph::DrawableGroup3D &drawables) : SceneGraph::Drawable3D{object, &drawables},
                                                                         _shader(shader),
                                                                         _mesh(mesh) {}

    void AtomDrawable::draw(const Matrix4 &transformation, SceneGraph::Camera3D &camera)
    {
        _shader
            .setTransformationMatrix(transformation)
            .setProjectionMatrix(camera.projectionMatrix())
            .setNormalMatrix(transformation.normalMatrix())
            .draw(_mesh);
    }
}