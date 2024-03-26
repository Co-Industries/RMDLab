
/* A B C D E F G H I J K L M N O P Q R S T U V W X Y Z */
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/Optional.h>
#include <Corrade/Containers/Pointer.h>
#include <Corrade/Utility/Arguments.h>

#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Renderer.h>

#include <Magnum/Math/Color.h>
#include <Magnum/Math/Vector.h>
#include <Magnum/Math/Vector3.h>

#include <Magnum/ImGuiIntegration/Context.hpp>
#include <Magnum/Platform/Sdl2Application.h>

#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/SceneGraph/Scene.h>

#include "./components/arcball/ArcBall.h"
#include "./components/arcball/ArcBallCamera.h"
#include "./content/Grid.h"
#include "./content/Skybox.h"
#include "./core/Simulation.h"
#include "./core/Simulation_old.h"
#include "Version.h"

namespace Magnum

{
    using namespace Math::Literals;

    class RMD : public Platform::Application
    {
    public:
        explicit RMD(const Arguments &arguments);

    private:
        void drawEvent() override;
        void viewportEvent(ViewportEvent &event) override;
        void keyPressEvent(KeyEvent &event) override;
        void keyReleaseEvent(KeyEvent &event) override;
        void mousePressEvent(MouseEvent &event) override;
        void mouseReleaseEvent(MouseEvent &event) override;
        void mouseMoveEvent(MouseMoveEvent &event) override;
        void mouseScrollEvent(MouseScrollEvent &event) override;
        void textInputEvent(TextInputEvent &event) override;

        ImGuiIntegration::Context _imgui{NoCreate};
        bool _showDemoWindow = false;
        Color3 _clearColor = Color3(1.0f, 1.0f, 0.0f);
        // Simulation
        Double _pvdW1 = 1.5591; // Double _pvdW1min = -5.0; Double _pvdW1max = 5.0;
        Double _cutoff_vpar30 = 0.0100;
        UnsignedInt _nso = 7, _nboty = 18;
        Double _plp1param = 3.5895;
        Double _povun3param = 38.5241,
               _povun4param = 3.4021,
               _povun6param = 1.0701,
               _povun7param = 11.9083,
               _povun8param = 13.3822;
        UnsignedLong _NATOMS = 100;

        Scene3D _scene;
        SceneGraph::DrawableGroup3D _drawables;
        Containers::Pointer<SimulationOld> _simulation;
        Containers::Pointer<Simulation> _newSimulation;
        Containers::Optional<ArcBallCamera> _arcballCamera;

        bool _drawOctreeBounds = true;
        bool _paused = true;
        bool _skipFrame = false;
    };
    RMD::RMD(const Arguments &arguments) : Platform::Application{arguments, NoCreate}
    {
        /* INFO Settings */
        Utility::Arguments args;
        args.addSkippedPrefix("magnum")
            .parse(arguments.argc, arguments.argv);
        /* INFO Window and parameters */
        {
            const Vector2 dpiScaling = this->dpiScaling({});
            Configuration conf;
            conf.setTitle(std::string("RMD ") + std::string(RMD_VERSION))
                .setSize(conf.size(), dpiScaling)
                .setWindowFlags(Configuration::WindowFlag::Resizable);
            GLConfiguration glConf;
            glConf.setSampleCount(dpiScaling.max() < 2.0f ? 8 : 2);
            if (!tryCreate(conf, glConf))
            {
                create(conf, glConf.setSampleCount(0));
            }
        }
        _imgui = ImGuiIntegration::Context(Vector2{windowSize()} / dpiScaling(),
                                           windowSize(), framebufferSize());

        GL::Renderer::setBlendEquation(GL::Renderer::BlendEquation::Add, GL::Renderer::BlendEquation::Add);
        GL::Renderer::setBlendFunction(GL::Renderer::BlendFunction::SourceAlpha, GL::Renderer::BlendFunction::OneMinusSourceAlpha);
        GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
        GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);

        /* INFO Initialize scene objects */
        {
            new Skybox(_scene, _drawables, 500.0f);
            new Grid(_scene, _drawables, 200.0f, Vector2i{40}, Color3{0.7f});
            _simulation.emplace(_scene, _drawables, UnsignedInt(100), _drawOctreeBounds);
            _newSimulation.emplace();
        }

        /* INFO Camera */
        {
            const Vector3 eye = Vector3::zAxis(500.0f);
            const Vector3 center{};
            const Vector3 up = Vector3::yAxis();
            const Deg fov = 45.0_degf;
            _arcballCamera.emplace(_scene, eye, center, up, fov, windowSize(), framebufferSize());
            _arcballCamera->setLagging(0.85f);
        }
        /* Loop at 60 Hz max (16)*/
        setSwapInterval(1);
        setMinimalLoopPeriod(16);
    };

    void RMD::drawEvent()
    {
        GL::defaultFramebuffer.clear(GL::FramebufferClear::Color | GL::FramebufferClear::Depth);

        if (!_paused || _skipFrame)
        {
            _skipFrame = false;

            _simulation->updateOctree();
            _simulation->updateAtoms();
        }

        /* Update camera */
        // ! unused variable [bool camChanged = _arcballCamera->update();]
        _arcballCamera->update();
        _arcballCamera->draw(_drawables);

        _imgui.newFrame();
        /* Enable text input, if needed */
        if (ImGui::GetIO().WantTextInput && !isTextInputActive())
        {
            startTextInput();
        }
        else if (!ImGui::GetIO().WantTextInput && isTextInputActive())
        {
            stopTextInput();
        }
        /* 1. Show a simple window. Tip: if we don't call ImGui::Begin()/ImGui::End() the widgets appear in
        a window called "Debug" automatically */
        {
            ImGui::Text("Auksts laiks, bet silta mana sirds.");
            if (ImGui::Button("Demo window"))
                _showDemoWindow ^= true;
            ImGui::Text("|| SIMULATION ||");
            // ImGui::SliderScalar("vdWaals shielding", ImGuiDataType_Double, &_pvdW1, &_pvdW1min, &_pvdW1max);
            // ImGui::InputScalar("Atom types", ImGuiDataType_U32, &_nso);
            // ImGui::InputScalar("Bonds", ImGuiDataType_U32, &_nboty);
            // ImGui::InputDouble("vdWaals shielding", &_pvdW1);
            // ImGui::InputDouble("Cutoff for bond order (*100)", &_cutoff_vpar30);
            // ImGui::InputDouble("Valency angle/lone pair parameter", &_plp1param);
            // ImGui::InputDouble("Overcoordination <povun3>", &_povun3param);
            // ImGui::InputDouble("Overcoordination <povun4>", &_povun4param);
            // ImGui::InputDouble("Undercoordination <povun6>", &_povun6param);
            // ImGui::InputDouble("Undercoordination <povun7>", &_povun7param);
            // ImGui::InputDouble("Undercoordination <povun8>", &_povun8param);
            if (ImGui::ColorEdit3("Atom color", _clearColor.data()))
                _simulation->updateColor(_clearColor);
            if (ImGui::Button("Run Simulation"))
                _newSimulation->RUN();
            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                        1000.0 / Double(ImGui::GetIO().Framerate), Double(ImGui::GetIO().Framerate));
        }

        if (_showDemoWindow)
        {
            ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiCond_FirstUseEver);
            ImGui::ShowDemoWindow();
        }

        /* Update application cursor */
        _imgui.updateApplicationCursor(*this);

        /* Set appropriate states. If you only draw ImGui, it is sufficient to
           just enable blending and scissor test in the constructor. */
        GL::Renderer::enable(GL::Renderer::Feature::Blending);
        GL::Renderer::enable(GL::Renderer::Feature::ScissorTest);
        GL::Renderer::disable(GL::Renderer::Feature::FaceCulling);
        GL::Renderer::disable(GL::Renderer::Feature::DepthTest);

        _imgui.drawFrame();

        GL::Renderer::disable(GL::Renderer::Feature::Blending);
        GL::Renderer::disable(GL::Renderer::Feature::ScissorTest);
        GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
        GL::Renderer::enable(GL::Renderer::Feature::DepthTest);

        swapBuffers();

        redraw();
    }

    void RMD::viewportEvent(ViewportEvent &event)
    {
        GL::defaultFramebuffer.setViewport({{}, event.framebufferSize()});
        _imgui.relayout(Vector2{event.windowSize()} / event.dpiScaling(),
                        event.windowSize(), event.framebufferSize());
        _arcballCamera->reshape(event.windowSize(), event.framebufferSize());
    }

    void RMD::keyPressEvent(KeyEvent &event)
    {
        switch (event.key())
        {
        case KeyEvent::Key::B:
            _drawOctreeBounds ^= true;
            break;
        case KeyEvent::Key::R:
            _arcballCamera->reset();
            break;
        case KeyEvent::Key::Space:
            _paused ^= true;
            break;
        case KeyEvent::Key::Right:
            _skipFrame = true;
            break;
        default:
            if (_imgui.handleKeyPressEvent(event))
            {
                event.setAccepted();
            }
        }
    }

    void RMD::keyReleaseEvent(KeyEvent &event)
    {
        if (_imgui.handleKeyReleaseEvent(event))
            return;
    }

    void RMD::mousePressEvent(MouseEvent &event)
    {
        if (_imgui.handleMousePressEvent(event))
            return;

        /* Enable mouse capture so the mouse can drag outside of the window */
        /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
        SDL_CaptureMouse(SDL_TRUE);

        if (event.modifiers() & MouseMoveEvent::Modifier::Ctrl)
            _arcballCamera->initTransformation(event.position(), 1);
        else if (event.modifiers() & MouseMoveEvent::Modifier::Alt)
            _arcballCamera->initTransformation(event.position(), 2);
        else
            _arcballCamera->initTransformation(event.position(), 0);
        event.setAccepted();
    }

    void RMD::mouseReleaseEvent(MouseEvent &event)
    {
        if (_imgui.handleMouseReleaseEvent(event))
            return;
        /* Disable mouse capture again */
        /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
        SDL_CaptureMouse(SDL_FALSE);
    }

    void RMD::mouseMoveEvent(MouseMoveEvent &event)
    {
        if (_imgui.handleMouseMoveEvent(event))
            return;

        if (!event.buttons())
            return;

        if (event.modifiers() & MouseMoveEvent::Modifier::Shift)
            _arcballCamera->translate(event.position());
        else if (event.modifiers() & MouseMoveEvent::Modifier::Ctrl)
            _arcballCamera->rotate(event.position(), 1);
        else if (event.modifiers() & MouseMoveEvent::Modifier::Alt)
            _arcballCamera->rotate(event.position(), 2);
        else
            _arcballCamera->rotate(event.position(), 0);
        event.setAccepted();
    }

    void RMD::mouseScrollEvent(MouseScrollEvent &event)
    {
        if (_imgui.handleMouseScrollEvent(event))
            return;

        const Float delta = event.offset().y();
        if (Math::abs(delta) < 1.0e-2f)
            return;

        _arcballCamera->zoom(20 * delta);

        event.setAccepted();
    }

    void RMD::textInputEvent(TextInputEvent &event)
    {
        if (_imgui.handleTextInputEvent(event))
            event.setAccepted();
        return;
    }
}

MAGNUM_APPLICATION_MAIN(Magnum::RMD)