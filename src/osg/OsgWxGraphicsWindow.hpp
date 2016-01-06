#ifndef _OSG_WX_GRRAPHICS_WINDOW_HPP
#define _OSG_WX_GRRAPHICS_WINDOW_HPP

#include "OsgWxGLCanvas.hpp"
#include <osgViewer/Viewer>

class OsgWxGraphicsWindow : public osgViewer::GraphicsWindow {
public:
    OsgWxGraphicsWindow(OsgWxGLCanvas *canvas);
    virtual ~OsgWxGraphicsWindow() { }
    void init();

    // GraphicsWindow interface
    void grabFocus();
    void grabFocusIfPointerInWindow();
    void useCursor(bool cursorOn);
    bool makeCurrentImplementation();
    void swapBuffersImplementation();

    // not implemented yet...just use dummy implementation to get working.
    virtual bool valid() const { return true; }
    virtual bool realizeImplementation() { return true; }
    virtual bool isRealizedImplementation() const  { return m_canvas->IsShownOnScreen(); }
    virtual void closeImplementation() { }
    virtual bool releaseContextImplementation() { return true; }
private:
    // XXX need to set _canvas to NULL when the canvas is deleted by
    // its parent. for this, need to add event handler in OsgWxGLCanvas
    OsgWxGLCanvas*  m_canvas;
};

#endif
