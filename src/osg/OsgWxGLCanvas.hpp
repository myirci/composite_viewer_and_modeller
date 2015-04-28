#ifndef _OSG_WX_CANVAS_HPP
#define _OSG_WX_CANVAS_HPP

#include "../modeller/components/GeneralizedCylinderGeometry.hpp"

#include <wx/glcanvas.h>
#include <osgViewer/Viewer>
#include <memory>

class OsgWxFrame;
class ImageModeller;
class CoordinateTransformations;
class OsgSelectionHandler;

template <typename T> class Point2D;

class OsgWxGLCanvas : public wxGLCanvas {
private:
    osg::ref_ptr<osgViewer::GraphicsWindow> m_graphics_window;
    wxGLContext* m_context;
    wxCursor m_oldCursor;
    OsgWxFrame* m_parent;
    ImageModeller* m_modeller;
    std::unique_ptr<OsgSelectionHandler> m_selection_handler;
public:
    OsgWxGLCanvas(wxWindow* parent, wxWindowID id = wxID_ANY, int* attributes = 0, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = 0, const wxString& name = wxT("wxOsgGLCanvas"));
    virtual ~OsgWxGLCanvas();
    void UsrSetGraphicsWindow(osgViewer::GraphicsWindow *gw) { m_graphics_window = gw; }
    void UsrUseCursor(bool value);
    void UsrMakeContextCurrent();
    void UsrInitializeModeller(const std::shared_ptr<CoordinateTransformations>& ppp, const wxString& fpath);
    void UsrAddSelectableNodeToDisplay(osg::Node* node, unsigned int component_id);
    void UsrAddToBackgroundDisplay(osg::Geometry* geom);
    void UsrTransformCoordinates(Point2D<int>& pt) const;
    void UsrLogicalToDevice(wxPoint& p) const;
    void UsrSaveModel(const wxString& path) const;
    void UsrLogErrorMessage(const std::string& str) const;
    void UsrSetRenderingType(rendering_type rtype);
    ImageModeller* UsrGetModeller();
    osg::Switch* UsrGetSelectionBoxes();
    const osg::Camera* const UsrGetMainCamera() const;

private:
    // Private member functions
    inline wxPoint usrDeviceToLogical(const wxPoint& p) const;

    // Event handlers
    void OnPaint(wxPaintEvent& event);
    void OnSize(wxSizeEvent& event);
    void OnEraseBackground(wxEraseEvent& event);
    void OnChar(wxKeyEvent &event);
    void OnKeyUp(wxKeyEvent &event);
    void OnMouseEnter(wxMouseEvent &event);
    void OnMouseDown(wxMouseEvent &event);
    void OnMouseUp(wxMouseEvent &event);
    void OnMouseMotion(wxMouseEvent &event);
    void OnMouseWheel(wxMouseEvent &event);
    DECLARE_EVENT_TABLE()
};

#endif
