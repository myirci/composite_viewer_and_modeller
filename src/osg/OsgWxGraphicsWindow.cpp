#include "OsgWxGraphicsWindow.hpp"

OsgWxGraphicsWindow::OsgWxGraphicsWindow(OsgWxGLCanvas* canvas) : m_canvas(canvas) {

    _traits = new GraphicsContext::Traits;
    wxPoint pos = m_canvas->GetPosition();
    wxSize  size = m_canvas->GetSize();
    _traits->x = pos.x;
    _traits->y = pos.y;
    _traits->width = size.x;
    _traits->height = size.y;
    init();
}

void OsgWxGraphicsWindow::init() {

    if (valid()) {
        setState(new osg::State );
        getState()->setGraphicsContext(this);
        if (_traits.valid() && _traits->sharedContext.valid()) {
            getState()->setContextID( _traits->sharedContext->getState()->getContextID() );
            incrementContextIDUsageCount( getState()->getContextID() );
        }
        else {
            getState()->setContextID( osg::GraphicsContext::createNewContextID() );
        }
    }
}

void OsgWxGraphicsWindow::grabFocus() {
    // focus the canvas
    m_canvas->SetFocus();
}

void OsgWxGraphicsWindow::grabFocusIfPointerInWindow() {

    // focus this window, if the pointer is in the window
    wxPoint pos = wxGetMousePosition();
    if (wxFindWindowAtPoint(pos) == m_canvas)
        m_canvas->SetFocus();
}

void OsgWxGraphicsWindow::useCursor(bool cursorOn) {

    m_canvas->UsrUseCursor(cursorOn);
}

bool OsgWxGraphicsWindow::makeCurrentImplementation() {

    m_canvas->UsrMakeContextCurrent();
    return true;
}

void OsgWxGraphicsWindow::swapBuffersImplementation() {

    m_canvas->SwapBuffers();
}

//bool OsgWxGraphicsWindow::resizedImplementation(int x, int y, int width, int height) {
//    std::cerr << "x: " << x << std::endl;
//    std::cerr << "y: " << y << std::endl;
//    std::cerr << "w: " << width << std::endl;
//    std::cerr << "h: " << height << std::endl;
//    return true;
//}
