
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
#include "./core/Data.h"
#include "./core/Simulation.h"
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

        ImGuiIntegration::Context imgui{NoCreate};
        bool showDemoWindow = false;

        Scene3D scene;
        SceneGraph::DrawableGroup3D drawables;
        Containers::Pointer<Simulation> simulation;

        //? Parameters
        std::size_t atomCount = 200;
        Float atomRadius = 0.3;
        Float randomVelocity = 2.0;
        Double timestep = 0.25;
        Double border = 20.0;

        Containers::Optional<ArcBallCamera> arcballCamera;

        bool paused = true;
        bool skipFrame = false;
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
        imgui = ImGuiIntegration::Context(Vector2{windowSize()} / dpiScaling(),
                                          windowSize(), framebufferSize());

        GL::Renderer::setBlendEquation(GL::Renderer::BlendEquation::Add, GL::Renderer::BlendEquation::Add);
        GL::Renderer::setBlendFunction(GL::Renderer::BlendFunction::SourceAlpha, GL::Renderer::BlendFunction::OneMinusSourceAlpha);
        GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
        GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);

        /* INFO Initialize scene objects */
        {
            new Skybox(scene, drawables, 5000.0f);
            new Grid(scene, drawables, 400.0f, Vector2i{31}, Color3{0.7f});
            simulation.emplace(scene, drawables);
        }

        /* INFO Camera */
        {
            const Vector3 eye = Vector3::zAxis(20.0f);
            const Vector3 center{};
            const Vector3 up = Vector3::yAxis();
            const Deg fov = 45.0_degf;
            arcballCamera.emplace(scene, eye, center, up, fov, windowSize(), framebufferSize());
            arcballCamera->setLagging(0.85f);
        }
        /* Loop at 60 Hz max (16)*/
        setSwapInterval(1);
        setMinimalLoopPeriod(16);
    };

    void RMD::drawEvent()
    {
        GL::defaultFramebuffer.clear(GL::FramebufferClear::Color | GL::FramebufferClear::Depth);

        if (!paused || skipFrame)
        {
            skipFrame = false;
            simulation->UPDATE_ATOMS();
        }

        /* Update camera */
        // ! unused variable [bool camChanged = arcballCamera->update();]
        arcballCamera->update();
        arcballCamera->draw(drawables);

        imgui.newFrame();
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
                showDemoWindow ^= true;
            ImGui::Text("SIMULATION");
            // if (ImGui::ColorEdit3("Atom color", _clearColor.data()))
            // oldSimulation->updateColor(_clearColor);
            ImGui::InputScalar("<NATOMS>", ImGuiDataType_U64, &atomCount);
            ImGui::InputFloat("<RADIUS>", &atomRadius);
            ImGui::InputFloat("<VELOCITY>", &randomVelocity);
            ImGui::InputDouble("<TIMESTEP>", &timestep);
            ImGui::InputDouble("<BORDER>", &border);
            if (ImGui::Button("Run Simulation"))
            {    
                simulation->RUN(SimulationParameters{atomCount, atomRadius, randomVelocity, timestep, border});
                paused = false;
            }
            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                        1000.0 / Double(ImGui::GetIO().Framerate), Double(ImGui::GetIO().Framerate));
        }

        if (showDemoWindow)
        {
            ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiCond_FirstUseEver);
            ImGui::ShowDemoWindow();
        }

        /* Update application cursor */
        imgui.updateApplicationCursor(*this);

        /* Set appropriate states. If you only draw ImGui, it is sufficient to
           just enable blending and scissor test in the constructor. */
        GL::Renderer::enable(GL::Renderer::Feature::Blending);
        GL::Renderer::enable(GL::Renderer::Feature::ScissorTest);
        GL::Renderer::disable(GL::Renderer::Feature::FaceCulling);
        GL::Renderer::disable(GL::Renderer::Feature::DepthTest);

        imgui.drawFrame();

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
        imgui.relayout(Vector2{event.windowSize()} / event.dpiScaling(),
                       event.windowSize(), event.framebufferSize());
        arcballCamera->reshape(event.windowSize(), event.framebufferSize());
    }

    void RMD::keyPressEvent(KeyEvent &event)
    {
        switch (event.key())
        {
        case KeyEvent::Key::B:
            drawOctreeBounds ^= true;
            break;
        case KeyEvent::Key::R:
            arcballCamera->reset();
            break;
        case KeyEvent::Key::Space:
            if (simulation->running)
                paused ^= true;
            break;
        case KeyEvent::Key::Right:
            if (simulation->running)
                skipFrame = true;
            break;
        default:
            if (imgui.handleKeyPressEvent(event))
            {
                event.setAccepted();
            }
        }
    }

    void RMD::keyReleaseEvent(KeyEvent &event)
    {
        if (imgui.handleKeyReleaseEvent(event))
            return;
    }

    void RMD::mousePressEvent(MouseEvent &event)
    {
        if (imgui.handleMousePressEvent(event))
            return;

        /* Enable mouse capture so the mouse can drag outside of the window */
        /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
        SDL_CaptureMouse(SDL_TRUE);

        if (event.modifiers() & MouseMoveEvent::Modifier::Ctrl)
            arcballCamera->initTransformation(event.position(), 1);
        else if (event.modifiers() & MouseMoveEvent::Modifier::Alt)
            arcballCamera->initTransformation(event.position(), 2);
        else
            arcballCamera->initTransformation(event.position(), 0);
        event.setAccepted();
    }

    void RMD::mouseReleaseEvent(MouseEvent &event)
    {
        if (imgui.handleMouseReleaseEvent(event))
            return;
        /* Disable mouse capture again */
        /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
        SDL_CaptureMouse(SDL_FALSE);
    }

    void RMD::mouseMoveEvent(MouseMoveEvent &event)
    {
        if (imgui.handleMouseMoveEvent(event))
            return;

        if (!event.buttons())
            return;

        if (event.modifiers() & MouseMoveEvent::Modifier::Shift)
            arcballCamera->translate(event.position());
        else if (event.modifiers() & MouseMoveEvent::Modifier::Ctrl)
            arcballCamera->rotate(event.position(), 1);
        else if (event.modifiers() & MouseMoveEvent::Modifier::Alt)
            arcballCamera->rotate(event.position(), 2);
        else
            arcballCamera->rotate(event.position(), 0);
        event.setAccepted();
    }

    void RMD::mouseScrollEvent(MouseScrollEvent &event)
    {
        if (imgui.handleMouseScrollEvent(event))
            return;

        const Float delta = event.offset().y();
        if (Math::abs(delta) < 1.0e-2f)
            return;

        arcballCamera->zoom(50 * delta);

        event.setAccepted();
    }

    void RMD::textInputEvent(TextInputEvent &event)
    {
        if (imgui.handleTextInputEvent(event))
            event.setAccepted();
        return;
    }
}

MAGNUM_APPLICATION_MAIN(Magnum::RMD)